'''
Workflow tools generate
Invoked: $ dane_wf wf: example <params/options/io>
'''
import logging
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from datetime import datetime
from pathlib import Path

from bioinformatics_tools.file_classes.base_classes import command
from bioinformatics_tools.workflow_tools.bapptainer import (
    CacheSifError, cache_sif_files, locate_local_sif_files)
from bioinformatics_tools.workflow_tools.models import WorkflowKey
from bioinformatics_tools.workflow_tools.output_cache import cached_tools, log_workflow_run, restore, restore_all, store, store_all
from bioinformatics_tools.workflow_tools.load_to_db import is_already_processed, PIPELINE_VERSION, compute_fasta_hash
from bioinformatics_tools.workflow_tools.programs import ProgramBase
from bioinformatics_tools.workflow_tools.workflow_helpers import (
    discover_genomes, get_workflow_prefix_for, WORKFLOW_PATH_DEFAULTS)
from bioinformatics_tools.workflow_tools.workflow_registry import (
    MARGIE_SB_PHASED_TOOLS,
    WORKFLOWS,
    margie_sb_sif_files,
)

LOGGER = logging.getLogger(__name__)

# Must mirror margie_sb.smk's INTERPRO_ANALYSIS_TO_BASENAME.values() exactly --
# duplicated here because margie_sb.smk isn't an importable Python module.
INTERPRO_DB_BASENAMES = [
    "antifam", "cdd", "coils", "funfam", "gene3d", "hamap", "mobidb", "ncbifam",
    "panther", "pfam", "pirsf", "pirsr", "prints", "prosite_patterns",
    "prosite_profiles", "sfld", "smart", "superfamily",
]
WORKFLOW_DIR = Path(__file__).parent

# --cores for SLURM-mode runs. NOT 'all', because Snakemake clamps every rule's
# threads: to --cores when it builds the DAG -- including rules it will hand to
# SLURM, which run on a different node than the driver and so have nothing to do
# with the driver's own core count. The driver now runs inside a small batch
# allocation (see ssh_slurm.DRIVER_CPUS), where 'all' resolves to that
# allocation's handful of cores; gtdbtk's threads: 64 silently became 4, and
# sbatch rejected the child job because Negishi's highmem partition requires at
# least 64 cores. Any ceiling at or above the largest rule's threads: restores
# the declared value; 256 is the biggest node on this cluster, so no rule can
# usefully ask for more. This does not oversubscribe the driver: in cluster mode
# it is --local-cores, defaulting to the host's real core count, that bounds the
# local rules actually running here, and --jobs that bounds the submitted ones.
SLURM_MODE_CORES = 256


def _subprocess_env() -> dict:
    '''Env for every snakemake subprocess we launch. Sets BASH_ENV to
    Negishi's Lmod init script so signalp4/signalp6 (the only two rules
    using envmodules: rather than container:) can find their "module" shell
    function -- Snakemake's --use-envmodules runs each rule via a non-login
    bash -c that never sources /etc/profile.d/*, so bash needs BASH_ENV to
    pick it up instead. Harmless no-op for every other workflow/rule.'''
    env = os.environ.copy()
    env['BASH_ENV'] = '/etc/profile.d/modules.sh'
    return env


def _snakemake_executable() -> str:
    '''Resolve snakemake next to the running interpreter rather than relying
    on subprocess PATH lookup of the bare name 'snakemake'. dane_wf's
    console-script entry point can be invoked by its absolute venv path
    (e.g. ~/bioinformatics-tools/.venv/bin/dane_wf) without that venv ever
    being activated -- PATH then has no reason to include the venv's bin/,
    even though snakemake is installed right alongside python there.
    '''
    candidate = Path(sys.executable).parent / 'snakemake'
    return str(candidate) if candidate.exists() else 'snakemake'


class _BackgroundProcess:
    '''Handle returned by WorkflowBase._start_background_subprocess(): wraps
    a Popen whose stdout/stderr are already being drained on background
    threads, so the caller can check .poll()/.wait() without itself
    blocking on log lines the way _run_subprocess()'s caller does.'''

    def __init__(self, proc: subprocess.Popen, stdout_lines: list[str], stderr_lines: list[str]):
        self.proc = proc
        self.stdout_lines = stdout_lines
        self.stderr_lines = stderr_lines

    def poll(self):
        '''None while still running, else the process's exit code.'''
        return self.proc.poll()

    def wait(self) -> int:
        return self.proc.wait()


class WorkflowBase(ProgramBase):
    '''Snakemake workflow execution. Inherits single-program commands from ProgramBase.
    '''

    def __init__(self, workflow_id=None):
        LOGGER.debug('Starting __init__ of WorkflowBase')
        self.workflow_id = workflow_id
        self.timestamp = datetime.now().strftime("%d%m%y-%H%M")

        LOGGER.debug('Using the workflow id of %s', self.workflow_id)

        super().__init__()

    def build_executable(self, key: WorkflowKey, config_overrides: dict = None, mode='notdev', compute_config: dict = None,
                          extra_resources: dict = None, rerun_triggers: str = None, target: str = None,
                          max_jobs_override: int = None) -> list[str]:
        '''
        Build snakemake command from workflow key and config.

        Args:
            key: WorkflowKey defining the workflow
            config_overrides: Only workflow-specific overrides (input_fasta, output_dir, main_database)
            mode: Execution mode ('dev' or other for slurm)
            compute_config: Compute cluster config (account, partition, resources)
            extra_resources: Optional named Snakemake resources (--resources k=v ...) to cap
                concurrency of a specific group of rules independently of --jobs, e.g.
                {'margie_sb_phase4_slot': 4} to let only 4 phase4 tools run at once.
            rerun_triggers: Optional value for Snakemake's --rerun-triggers (e.g. 'mtime').
                None (default) leaves the combined input+code+params triggers untouched;
                resume_job's relaunch asks for 'mtime' since a resumed run's output_dir
                always has a new timestamp, which combined triggers would otherwise
                treat as a params change and rerun everything.
            target: Optional explicit Snakemake target rule name (e.g. 'rasttk_all' or
                'phase4_12_one_genome' -- see margie_sb.smk) instead of the implicit
                rule all. Must come before --config (which greedily consumes every
                following token as key=value). Also passes --nolock, since
                _run_pipeline_batch_sequential runs 'rasttk_all' concurrently with a
                sequence of 'phase4_12_one_genome' targets against the same working
                directory -- safe because the two touch disjoint rule sets and outputs.
            max_jobs_override: Caller-specified --jobs value that wins over
                compute_config's max_jobs. Used by Stage 1's rasttk_all invocation
                (see _run_pipeline_batch_sequential) to force real BV-BRC submissions
                one at a time -- a named Snakemake resource pool was tried for this
                instead (margie_sb_phase3_slot) and deadlocked a real run (pool stuck
                at 0 with jobs never selected), so concurrency is capped here via
                Snakemake's own --jobs scheduler instead, which every cluster executor
                already has to get right. Safe because Stage 1's whole DAG is just
                gtdbtk+rasttk -- nothing else needs the higher Stage 2 concurrency.
        '''
        smk_path = WORKFLOW_DIR / key.snakemake_file

        # Use compute config to determine max_jobs (default to 5)
        max_jobs = 5
        if compute_config:
            max_jobs = compute_config.get('max_jobs', 5)
        if max_jobs_override is not None:
            max_jobs = max_jobs_override

        # 'all' only in dev mode, where every rule really does run on this
        # machine and should be sized to it. See SLURM_MODE_CORES.
        cores = 'all' if mode == 'dev' else SLURM_MODE_CORES

        core_command = [
            _snakemake_executable(),
            '-s', str(smk_path),
            f'--cores={cores}',
            '--keep-going',
            '--use-apptainer',
            '--sdm=apptainer',
            '--apptainer-args', '-B /home/ddeemer -B /depot/lindems/data/Databases/',  #TODO: HARDCODED!
            # signalp4/signalp6 are the only two rules using envmodules: instead
            # of container: (the cluster already provides them, see margie_sb.smk).
            # Without this flag those rules run with neither signalp4 nor
            # signalp6 on PATH.
            '--use-envmodules',
            # Prints each job's rule/wildcards/jobid block to its own live log --
            # the only place run_ssh_task can reliably read a job's genome
            # wildcard from, instead of racing the remote log file's cleanup.
            '--verbose',
            f'--jobs={max_jobs}',
            '--latency-wait=60',
            '--scheduler=greedy',
        ]

        if target:
            core_command.append('--nolock')
            core_command.append(target)

        # Check for dry-run/test mode (from config or command line)
        if self.conf.get('dry_run', False) or self.conf.get('test_only', False):
            core_command.append('--dry-run')
            LOGGER.info("Running in DRY-RUN mode - no actual execution")

        if mode != 'dev':
            core_command.append('--executor=slurm')

        # Add default SLURM resources from compute config
        if mode != 'dev' and compute_config:
            default_resources = ['--default-resources']

            # Required: account
            account = compute_config.get('account', '').strip()
            if account:
                default_resources.append(f'slurm_account={account}')

            # Optional: partition
            partition = compute_config.get('partition', '').strip()
            if partition:
                default_resources.append(f'slurm_partition={partition}')

            # Optional: default runtime and memory
            if 'default_runtime' in compute_config:
                default_resources.append(f'runtime={compute_config["default_runtime"]}')

            if 'default_mem_mb' in compute_config:
                default_resources.append(f'mem_mb={compute_config["default_mem_mb"]}')

            # Only add if we have at least the account
            if len(default_resources) > 1:
                core_command.extend(default_resources)

        # Pass original config file(s) to Snakemake to preserve types
        # This loads the full user config with proper int/bool/nested dict types
        if hasattr(self, 'config_paths') and self.config_paths:
            for config_path in self.config_paths:
                core_command.extend(['--configfile', str(config_path)])

        # Override only workflow-specific values (all strings, so no type issues)
        if config_overrides:
            config_pairs = [f'{k}={v}' for k, v in config_overrides.items()]
            core_command.append('--config')
            core_command.extend(config_pairs)

        # Named resources cap concurrency of a specific group of rules,
        # independent of --jobs (which governs everything else).
        if extra_resources:
            core_command.append('--resources')
            core_command.extend(f'{k}={v}' for k, v in extra_resources.items())

        if rerun_triggers:
            core_command.extend(['--rerun-triggers', rerun_triggers])

        return core_command

    @staticmethod
    def _parse_snakemake_output(stderr: str) -> dict:
        '''Best-effort parse of snakemake stderr for structured reporting.'''
        result = {'total': 0, 'completed': 0, 'failed': 0, 'failed_rules': []}

        # Extract "X of Y steps (Z%) done"
        steps_match = re.search(r'(\d+) of (\d+) steps \(\d+%\) done', stderr)
        if steps_match:
            result['completed'] = int(steps_match.group(1))
            result['total'] = int(steps_match.group(2))

        # Extract failed rule names from "Error in rule <name>:"
        failed_rules = re.findall(r'Error in rule (\w+):', stderr)
        result['failed_rules'] = failed_rules
        result['failed'] = len(failed_rules)

        # If we found failed rules but no total, estimate total from completed + failed
        if result['failed'] and not result['total']:
            result['total'] = result['completed'] + result['failed']

        return result

    def _run_subprocess(self, wf_command):
        '''Wrapper for subprocess.run(). Returns CompletedProcess on any exit
        code (even non-zero), or None on launch failure (e.g. snakemake not installed).'''
        LOGGER.debug('Received command and running: %s', wf_command)

        # Pin snakemake's working directory to output_dir so that .snakemake/
        # and any relative rule paths resolve there, regardless of the SSH
        # session's CWD on the cluster.
        output_dir = self.conf.get('output_dir', '')
        cwd = output_dir or None
        if cwd:
            Path(cwd).mkdir(parents=True, exist_ok=True)

        try:
            proc = subprocess.Popen(
                wf_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, cwd=cwd, env=_subprocess_env(),
            )

            # Collect stderr on a background thread so it doesn't block stdout reads.
            stderr_lines: list[str] = []

            def _read_stderr():
                for line in proc.stderr:
                    line = line.rstrip()
                    LOGGER.info('[snakemake] %s', line)
                    stderr_lines.append(line)

            stderr_thread = threading.Thread(target=_read_stderr, daemon=True)
            stderr_thread.start()

            stdout_lines: list[str] = []
            for line in proc.stdout:
                line = line.rstrip()
                LOGGER.info('[snakemake] %s', line)
                stdout_lines.append(line)

            stderr_thread.join()
            proc.wait()

            return subprocess.CompletedProcess(
                args=wf_command,
                returncode=proc.returncode,
                stdout='\n'.join(stdout_lines),
                stderr='\n'.join(stderr_lines),
            )
        except Exception as e:
            LOGGER.error('Failed to launch subprocess %s: %s', wf_command[0], e)
            self.failed(f'Failed to launch subprocess: {e}')
            return None

    def _start_background_subprocess(self, wf_command) -> '_BackgroundProcess | None':
        '''Non-blocking sibling of _run_subprocess() -- starts wf_command the
        same way (cwd pinned to output_dir, stdout/stderr streamed
        line-by-line to LOGGER as "[snakemake] ..." so job_runner.py's SSH
        log parsing sees it exactly the same way it already does for every
        other snakemake call) but returns immediately instead of joining/
        waiting, since the caller (_run_pipeline_batch_sequential's Stage 1)
        needs to keep polling for newly-ready genomes while this keeps
        running in the background. Returns None on launch failure, same
        convention as _run_subprocess().'''
        LOGGER.debug('Received command and running in background: %s', wf_command)

        output_dir = self.conf.get('output_dir', '')
        cwd = output_dir or None
        if cwd:
            Path(cwd).mkdir(parents=True, exist_ok=True)

        try:
            proc = subprocess.Popen(
                wf_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, cwd=cwd, env=_subprocess_env(),
            )
        except Exception as e:
            LOGGER.error('Failed to launch background subprocess %s: %s', wf_command[0], e)
            return None

        stdout_lines: list[str] = []
        stderr_lines: list[str] = []

        def _read_stdout():
            for line in proc.stdout:
                line = line.rstrip()
                LOGGER.info('[snakemake] %s', line)
                stdout_lines.append(line)

        def _read_stderr():
            for line in proc.stderr:
                line = line.rstrip()
                LOGGER.info('[snakemake] %s', line)
                stderr_lines.append(line)

        threading.Thread(target=_read_stdout, daemon=True).start()
        threading.Thread(target=_read_stderr, daemon=True).start()

        return _BackgroundProcess(proc, stdout_lines, stderr_lines)

    def _build_result(self, key_name, proc):
        '''Build a structured result dict from a completed snakemake process.'''
        rules_summary = self._parse_snakemake_output(proc.stderr)
        return {
            'workflow': key_name,
            'returncode': proc.returncode,
            'rules_summary': rules_summary,
            'stdout_tail': proc.stdout[-2000:] if proc.stdout else '',
            'stderr_tail': proc.stderr[-2000:] if proc.stderr else '',
        }

    def _output_prefix(self) -> str:
        """Return the filesystem prefix to prepend to all output paths for this run.

        Reads ``output_dir`` from the caragols config (set via CLI arg or passed
        from the API). If present, returns ``'{output_dir}/'``; otherwise returns
        ``''`` so that output paths remain relative to the SSH working directory.

        Every ``do_*`` workflow method should call this instead of reading
        ``output_dir`` directly, so the logic stays in one place.
        """
        output_dir = self.conf.get('output_dir', '')
        return f"{output_dir.rstrip('/')}/" if output_dir else ''

    def _run_pipeline(self, key_name: str, smk_config: dict, cache_map: dict = None, mode='dev', compute_config: dict = None):
        '''Shared pipeline execution: cache containers, restore outputs, run snakemake, store outputs.'''
        run_id = str(uuid.uuid4())
        LOGGER.info('Starting workflow "%s" run_id=%s', key_name, run_id)

        selected_wf = WORKFLOWS.get(key_name)
        if not selected_wf:
            self.failed(f'No workflow key found for "{key_name}"')
            return 1

        # Download / ensure .sif files are cached (skip if none needed, e.g. selftest)
        if selected_wf.sif_files:
            if selected_wf.local_sif_only:
                # Never contact the registry — just report what's already on disk.
                locate_local_sif_files(selected_wf.sif_files, local_sif_dir=self.conf.get('sif_path', None))
            else:
                try:
                    cache_sif_files(selected_wf.sif_files, local_sif_dir=self.conf.get('sif_path', None))
                except CacheSifError as e:
                    LOGGER.critical('Error with cache_sif_files: %s', e)
                    self.failed(f'Error with cache_sif_files: {e}')
                    return 1

        # Restore cached outputs from DB so snakemake skips completed rules
        db_path = smk_config.get('main_database')
        input_file = smk_config.get('input_fasta') or smk_config.get('input_file')
        restored = {}
        if cache_map and db_path and input_file:
            restored = restore_all(db_path, input_file, cache_map)
            LOGGER.info('Cache restore results: %s', restored)

        # Build and run snakemake
        wf_command = self.build_executable(selected_wf, config_overrides=smk_config, mode=mode, compute_config=compute_config)
        LOGGER.info('Running snakemake command: %s', ' '.join(wf_command))
        proc = self._run_subprocess(wf_command)

        # Launch failure (e.g. snakemake not installed) — already called self.failed()
        if proc is None:
            return 1

        result = self._build_result(key_name, proc)

        if proc.returncode != 0:
            LOGGER.error('Snakemake failed (rc=%d): %s', proc.returncode, result['rules_summary'])
            if cache_map and db_path and input_file:
                log_workflow_run(db_path, run_id, input_file, key_name,
                                 result['rules_summary'].get('completed', 0), status='failed')
            self.failed(msg=f'Workflow "{key_name}" failed', dex=result)
            return proc.returncode

        # Success — store outputs and log the run
        if cache_map and db_path and input_file:
            # Only store outputs that were cache misses (newly computed)
            tools_to_store = {tool: paths for tool, paths in cache_map.items()
                            if not restored.get(tool, False)}
            if tools_to_store:
                store_all(db_path, input_file, tools_to_store)
            else:
                LOGGER.info('All outputs were cache hits — skipping redundant storage')
            log_workflow_run(db_path, run_id, input_file, key_name,
                             result['rules_summary'].get('completed', 0), status='success')

        self.succeeded(msg=f'Workflow "{key_name}" completed successfully', dex=result)

    def _run_pipeline_batch(self, key_name: str, smk_config: dict, genome_files: dict[str, str],
                             cache_tool: str, cache_paths_fn, mode='dev', compute_config: dict = None,
                             extra_resources: dict = None, sif_files_override: list[tuple] = None,
                             rerun_triggers: str = None):
        '''Batch sibling of _run_pipeline() for workflows that accept a folder of
        genomes (WorkflowKey.supports_batch_input=True). Launches exactly ONE
        snakemake subprocess covering every genome via Snakemake's own {genome}
        wildcard DAG — parallelism across genomes is handled by Snakemake's
        --jobs scheduler, not by looping subprocess calls here. The per-genome
        SQLite output cache (output_cache.py) is still checked/updated once per
        genome, since each genome has its own content hash.

        cache_paths_fn(stem) -> list[str] must return that genome's output
        files to cache, mirroring whatever path shape the .smk file writes.

        sif_files_override, when given, replaces selected_wf.sif_files for
        this call only -- lets a caller validate only the SIFs needed for a
        partial-phase run instead of the workflow's full, static list.
        '''
        run_id = str(uuid.uuid4())
        LOGGER.info('Starting batch workflow "%s" run_id=%s genomes=%d', key_name, run_id, len(genome_files))

        selected_wf = WORKFLOWS.get(key_name)
        if not selected_wf:
            self.failed(f'No workflow key found for "{key_name}"')
            return 1

        sif_files = selected_wf.sif_files if sif_files_override is None else sif_files_override
        if sif_files:
            local_sif_dir = (self.conf.get(key_name, {}).get('sif_path', None)
                            or WORKFLOW_PATH_DEFAULTS.get(key_name, {}).get('sif_path'))
            if selected_wf.local_sif_only:
                locate_local_sif_files(sif_files, local_sif_dir=local_sif_dir)
            else:
                try:
                    cache_sif_files(sif_files, local_sif_dir=local_sif_dir)
                except CacheSifError as e:
                    LOGGER.critical('Error with cache_sif_files: %s', e)
                    self.failed(f'Error with cache_sif_files: {e}')
                    return 1

        db_path = smk_config.get('main_database')

        # Per-genome cache restore (cheap SQLite check, before Snakemake runs)
        restored: dict[str, bool] = {}
        if db_path:
            for stem, genome_file in genome_files.items():
                cache_map = {cache_tool: cache_paths_fn(stem)}
                restored[stem] = restore_all(db_path, genome_file, cache_map).get(cache_tool, False)
            LOGGER.info('Cache restore results: %s', restored)

        # One snakemake subprocess for the whole batch.
        wf_command = self.build_executable(selected_wf, config_overrides=smk_config, mode=mode, compute_config=compute_config,
                                           extra_resources=extra_resources, rerun_triggers=rerun_triggers)
        LOGGER.info('Running snakemake command: %s', ' '.join(wf_command))
        proc = self._run_subprocess(wf_command)

        if proc is None:
            return 1

        result = self._build_result(key_name, proc)

        if proc.returncode != 0:
            LOGGER.error('Snakemake failed (rc=%d): %s', proc.returncode, result['rules_summary'])
            if db_path:
                for stem, genome_file in genome_files.items():
                    log_workflow_run(db_path, run_id, genome_file, key_name,
                                     result['rules_summary'].get('completed', 0), status='failed')
            self.failed(msg=f'Workflow "{key_name}" failed', dex=result)
            return proc.returncode

        # Success — store outputs for genomes that were cache misses, log every genome's run
        if db_path:
            for stem, genome_file in genome_files.items():
                if not restored.get(stem, False):
                    store_all(db_path, genome_file, {cache_tool: cache_paths_fn(stem)})
                log_workflow_run(db_path, run_id, genome_file, key_name,
                                 result['rules_summary'].get('completed', 0), status='success')

        self.succeeded(msg=f'Workflow "{key_name}" completed successfully ({len(genome_files)} genome(s))', dex=result)

    @staticmethod
    def _gtdbtk_split_paths(genome: str, smk_config: dict) -> list[str]:
        prefix = get_workflow_prefix_for(genome, smk_config)
        return [f"{prefix}gtdbtk/gtdbtk_results.tsv", f"{prefix}gtdbtk/translation_table.tsv"]

    @staticmethod
    def _rasttk_paths(genome: str, smk_config: dict) -> list[str]:
        prefix = get_workflow_prefix_for(genome, smk_config)
        return [f"{prefix}rasttk/rast.tsv", f"{prefix}rasttk/rast.faa", f"{prefix}rasttk/rast.gff",
                f"{prefix}rasttk/rasttk_db.tkn"]

    @staticmethod
    def _restore_gtdbtk_batch(db_path: str, genome_files: dict[str, str], smk_config: dict) -> bool:
        '''All-or-nothing cache check for GTDB-Tk's batch rule (run_gtdbtk_batch
        in margie_sb.smk): GTDB-Tk runs once across every genome in the set, so
        a per-genome cache hit alone isn't enough -- Snakemake's DAG still sees
        run_gtdbtk_batch's own combined outputs as missing and reruns the whole
        batch, overwriting any per-genome split files restored individually.

        If every genome's split files are individually cached, restores them
        all, then synthesizes run_gtdbtk_batch's own combined files by
        concatenating the rows just written. A single genome miss aborts the
        shortcut; the real batch rerun then overwrites the partial restore
        with a fresher mtime, so downstream rules still see the change.
        '''
        if not db_path:
            return False

        for genome, fasta_path in genome_files.items():
            if not restore(db_path, fasta_path, 'gtdbtk', WorkflowBase._gtdbtk_split_paths(genome, smk_config)):
                return False

        output_dir = (smk_config.get('output_dir') or '').rstrip('/')
        if not output_dir:
            return False
        batch_dir = f"{output_dir}/original_container_outputs/gtdbtk"
        Path(batch_dir).mkdir(parents=True, exist_ok=True)

        result_header, result_rows = None, []
        translation_header, translation_rows = None, []
        for genome in genome_files:
            results_path, translation_path = WorkflowBase._gtdbtk_split_paths(genome, smk_config)
            r_header, *r_rows = Path(results_path).read_text().splitlines()
            t_header, *t_rows = Path(translation_path).read_text().splitlines()
            result_header, translation_header = r_header, t_header
            # Force column 0 (genome name) to this loop's own genome rather
            # than trusting whatever name is embedded in the restored row --
            # output_cache keys purely on FASTA content hash, so two
            # byte-identical genomes under different names (e.g. a synteny
            # reference copy) collapse to one cached row and would otherwise
            # silently mislabel one of them, leaving the other entirely
            # absent from this combined file (split_gtdbtk_batch_per_genome
            # then fails with "Missing GTDB-Tk row" for the absent one).
            if r_rows:
                r_rows = ["\t".join([genome, *r_rows[0].split("\t")[1:]])]
            if t_rows:
                t_rows = ["\t".join([genome, *t_rows[0].split("\t")[1:]])]
            result_rows.extend(r_rows)
            translation_rows.extend(t_rows)

        Path(f"{batch_dir}/gtdbtk_results.tsv").write_text(
            "\n".join([result_header, *result_rows]) + "\n")
        Path(f"{batch_dir}/gtdbtk.translation_table_summary.tsv").write_text(
            "\n".join([translation_header, *translation_rows]) + "\n")
        Path(f"{batch_dir}/gtdbtk_batch.done").touch()

        # Bump every per-genome split file's mtime past what we just wrote
        # for the batch files above -- otherwise Snakemake's mtime rerun
        # trigger sees split_gtdbtk_batch_per_genome's input (these batch
        # files) as newer than its own already-correct output (the
        # per-genome files), reruns it for no reason, and that rewrite then
        # makes run_rasttk's gtdbtk input look newer than rasttk's own
        # cache-restored output too -- cascading into a real (and possibly
        # failing) RASTtk/BV-BRC call for a genome that was already cached.
        for genome in genome_files:
            results_path, translation_path = WorkflowBase._gtdbtk_split_paths(genome, smk_config)
            Path(results_path).touch()
            Path(translation_path).touch()

        LOGGER.info('GTDB-Tk batch cache HIT for all %d genomes — skipping the container run entirely', len(genome_files))
        # Only logged once the WHOLE batch is confirmed a hit -- logging
        # per-genome inside the loop above would have been misleading for a
        # batch that ends up a partial miss, since the real rerun overwrites
        # every genome's restored files regardless of its own individual hit.
        for genome in genome_files:
            LOGGER.info("Cache HIT for gtdbtk (genome=%s) — skipping recomputation", genome)
        return True

    def _genome_stage2_progress(self, genome: str, genome_file: str, smk_config: dict) -> tuple[int, int]:
        '''Cheap, read-only "how far along is this genome" proxy used to run the
        most-complete genomes first in Stage 2 (see _run_pipeline_batch_sequential).

        A phase is "done" for this genome if EITHER its per-tool output cache row
        exists in the DB (a prior run completed it -- the common resume case,
        where the output directory is still empty and everything shows as
        CACHED) OR its '{output_dir}/{genome}/{tool_key}/' subdir is already
        populated on disk (resume-from-disk). We union both so the scheduler
        sees real progress in either scenario -- the earlier disk-only version
        scored every cache-restored genome as 0 and fell back to FIFO. Returns
        (max_phase_reached, n_phases_done); higher means fewer phases left, i.e.
        process sooner. Never raises -- an unreadable genome scores (0, 0) and
        sorts last.'''
        done: set = set()
        # (1) phases already in the output cache (completed in a prior run)
        db_path = smk_config.get('main_database')
        if db_path:
            done |= cached_tools(db_path, genome_file)
        # (2) phases already materialised on disk this run
        try:
            prefix = Path(get_workflow_prefix_for(genome, smk_config))
            for tool in MARGIE_SB_PHASED_TOOLS:
                key = tool.get('key')
                if not key or key in done:
                    continue
                d = prefix / key
                try:
                    if d.is_dir() and any(d.iterdir()):
                        done.add(key)
                except OSError:
                    continue
        except Exception:
            pass
        phase_of = {t.get('key'): t.get('phase', 0) for t in MARGIE_SB_PHASED_TOOLS}
        max_phase = max((phase_of.get(t, 0) for t in done), default=0)
        return (max_phase, len(done))

    def _run_pipeline_batch_sequential(self, key_name: str, smk_config: dict, genome_files: dict[str, str],
                                        cache_map_fn, mode='dev', compute_config: dict = None,
                                        extra_resources: dict = None, sif_files_override: list[tuple] = None,
                                        rerun_triggers: str = None):
        '''Sequential-per-organism sibling of _run_pipeline_batch(), used by
        do_margie_sb() for margie_sb's phase-ordering requirement: RASTtk/
        GTDB-Tk (phase1-3) are bottlenecked on BV-BRC's remote service, so
        they run breadth-first across every genome via one long-running
        Snakemake invocation targeting rule rasttk_all (Stage 1). Every
        local-compute phase (4-8) plus that genome's consolidation (phase9),
        labeling (phase10), and scoring (phase11) instead processes ONE
        genome fully via rule phase4_12_one_genome (Stage 2) before starting
        the next.

        Stage 1 is launched with max_jobs_override=1 (see build_executable):
        its whole DAG is just gtdbtk+rasttk, and only one genome's real
        BV-BRC submission may ever be in flight (margie_sb.smk's run_rasttk
        also holds its own mkdir mutex around the BV-BRC call as a backstop,
        but without this --jobs=1 cap Snakemake still submits up to
        compute_config's max_jobs worth of rasttk SLURM jobs that just idle
        on that mutex, wasting cluster allocation and starving other phases
        of account quota -- observed live on 2026-06-26).

        Genomes become eligible for Stage 2 in the order their RASTtk+GTDB-Tk
        output became ready (polled via each genome's RASTtk token file), since
        that's the only order BV-BRC's queue can be observed in. Among the
        genomes ready at any moment, though, Stage 2 processes the MOST-COMPLETE
        one first (most phase outputs already on disk -- see
        _genome_stage2_progress): time-to-first-result matters, so a genome that
        is minutes from a final result must not wait behind one that needs hours,
        and the slowest (least-complete) genomes are left for last. A genome's
        Stage 2 failure halts the queue there -- deliberate: later genomes, even
        ones already RASTtk-ready, are not processed this run.

        Same setup (sif validation, per-genome output_cache restore/store)
        as _run_pipeline_batch() -- only the "how Snakemake actually runs"
        middle section differs.
        '''
        run_id = str(uuid.uuid4())
        LOGGER.info('Starting sequential batch workflow "%s" run_id=%s genomes=%d', key_name, run_id, len(genome_files))

        selected_wf = WORKFLOWS.get(key_name)
        if not selected_wf:
            self.failed(f'No workflow key found for "{key_name}"')
            return 1

        sif_files = selected_wf.sif_files if sif_files_override is None else sif_files_override
        if sif_files:
            local_sif_dir = (self.conf.get(key_name, {}).get('sif_path', None)
                            or WORKFLOW_PATH_DEFAULTS.get(key_name, {}).get('sif_path'))
            if selected_wf.local_sif_only:
                locate_local_sif_files(sif_files, local_sif_dir=local_sif_dir)
            else:
                try:
                    cache_sif_files(sif_files, local_sif_dir=local_sif_dir)
                except CacheSifError as e:
                    LOGGER.critical('Error with cache_sif_files: %s', e)
                    self.failed(f'Error with cache_sif_files: {e}')
                    return 1

        db_path = smk_config.get('main_database')

        restored: dict[str, dict[str, bool]] = {}
        rasttk_restored: dict[str, bool] = {}
        gtdbtk_batch_hit = False
        if db_path:
            # Phase4+ cache restore is deferred to just before each genome's
            # Stage 2 (see the per-genome restore below). Restoring here, before
            # Stage 1 starts, causes mtime violations: rasttk runs during Stage 1
            # and produces rast.faa with a newer mtime than the restored
            # tmbed_results.tsv / interpro_results.tsv / etc., so Snakemake's
            # mtime trigger fires and re-runs every phase4 tool despite the cache.
            LOGGER.info('Phase4+ cache restore deferred to per-genome pre-Stage-2 (after rasttk outputs are fresh)')

            # gtdbtk/rasttk (phase1-3) are excluded from cache_map_fn -- see
            # _genome_cache_map's docstring. GTDB-Tk's real rule runs once
            # across the whole genome set, so _restore_gtdbtk_batch handles
            # that batch shape directly; rasttk is a normal per-genome rule
            # once gtdbtk's split outputs are in place.
            gtdbtk_batch_hit = self._restore_gtdbtk_batch(db_path, genome_files, smk_config)
            if gtdbtk_batch_hit:
                # Only worth attempting once gtdbtk's batch is a full hit --
                # otherwise gtdbtk reruns for real and overwrites these with
                # a fresher mtime anyway, forcing rasttk to rerun too.
                for stem, genome_file in genome_files.items():
                    hit = restore(db_path, genome_file, 'rasttk', self._rasttk_paths(stem, smk_config))
                    rasttk_restored[stem] = hit
                    if hit:
                        LOGGER.info("Cache HIT for rasttk (genome=%s) — skipping recomputation", stem)
            LOGGER.info('GTDB-Tk batch cache hit: %s. RASTtk per-genome cache restore results: %s',
                        gtdbtk_batch_hit, rasttk_restored)

        # Stage 1: one long-running snakemake invocation, phase1-3 only,
        # every genome -- runs breadth-first in the background for the rest
        # of this method's lifetime. max_jobs_override=1 keeps real BV-BRC
        # submissions one genome at a time (see this method's own docstring).
        stage1_command = self.build_executable(selected_wf, config_overrides=smk_config, mode=mode,
                                                compute_config=compute_config, extra_resources=extra_resources,
                                                rerun_triggers=rerun_triggers, target='rasttk_all',
                                                max_jobs_override=1)
        LOGGER.info('Starting Stage 1 (phase1-3, all genomes): %s', ' '.join(stage1_command))
        stage1 = self._start_background_subprocess(stage1_command)
        if stage1 is None:
            self.failed(f'Workflow "{key_name}" failed to launch Stage 1 (rasttk_all)')
            return 1

        def _rasttk_token_path(genome: str) -> str:
            return f"{get_workflow_prefix_for(genome, smk_config)}rasttk/rasttk_db.tkn"

        # When LLM is enabled, Stage 2 stops before the GPU step so CPU phases
        # for genome N+1 can run while genome N's LLM waits in the SLURM queue.
        # Stage 3 (rule llm_all) then submits all LLM jobs at once; SLURM's
        # gres=gpu:1 constraint serialises them on the GPU automatically.
        _run_llm_val = smk_config.get('run_llm', False)
        run_llm_enabled = _run_llm_val not in (False, 'false', '0', 'no')
        stage2_target = 'phase4_12_one_genome_no_llm' if run_llm_enabled else 'phase4_12_one_genome'
        LOGGER.info('LLM enabled: %s  →  Stage 2 target: %s', run_llm_enabled, stage2_target)

        total = len(genome_files)
        pending = set(genome_files.keys())
        queue: list[str] = []
        skipped: list[str] = []
        processed = 0
        last_proc = None
        progress_memo: dict[str, tuple[int, int]] = {}  # genome -> completeness, for most-complete-first ordering

        while pending or queue:
            for genome in list(pending):
                if Path(_rasttk_token_path(genome)).exists():
                    pending.discard(genome)
                    queue.append(genome)
                    if db_path:
                        if not gtdbtk_batch_hit:
                            store(db_path, genome_files[genome], 'gtdbtk', self._gtdbtk_split_paths(genome, smk_config))
                        if not rasttk_restored.get(genome, False):
                            store(db_path, genome_files[genome], 'rasttk', self._rasttk_paths(genome, smk_config))

            if not queue:
                if stage1.poll() is not None and pending:
                    # Stage 1 exited and these genomes never produced a token --
                    # a real phase1-3 failure for them specifically. --keep-going
                    # already let Stage 1 continue past them; skip and move on.
                    LOGGER.warning('Stage 1 exited without producing a RASTtk token for: %s -- skipping',
                                   sorted(pending))
                    skipped.extend(pending)
                    pending.clear()
                else:
                    time.sleep(15)
                continue

            # Most-complete-first: among the RASTtk-ready genomes, process the
            # one closest to a final result next, so quick wins (genomes minutes
            # from done -- e.g. all phase4-8 tools already cached, needing only
            # consolidation/labeling/fingerprint/scoring) are not stuck behind a
            # genome that needs hours; the least-complete genomes are left for
            # last. Progress is stable for genomes still pending (nothing writes
            # to their cache until they are processed), so memoise it. Python's
            # sort is stable, so genomes at equal progress keep their
            # RASTtk-ready (FIFO) order.
            if len(queue) > 1:
                queue.sort(
                    key=lambda g: progress_memo.setdefault(
                        g, self._genome_stage2_progress(g, genome_files[g], smk_config)),
                    reverse=True,
                )
            genome = queue.pop(0)
            processed += 1
            LOGGER.info('=== SEQUENTIAL: genome %d/%d (%s) phase4-12 starting ===', processed, total, genome)

            # Restore phase4+ cached outputs NOW — after rasttk has run and
            # produced rast.faa — so the restored files have a newer mtime than
            # rast.faa and Snakemake's mtime trigger correctly skips them.
            # MUST run BEFORE the already-processed skip below: otherwise a
            # genome already at PIPELINE_VERSION would `continue` past this and
            # a fresh output_dir would never get its consolidation/phase4-8
            # files materialized on disk (only Stage 1's rasttk output would
            # land), which is exactly the "restores only up to rasttk" bug.
            if db_path:
                genome_restored = restore_all(db_path, genome_files[genome], cache_map_fn(genome))
                restored[genome] = genome_restored
                LOGGER.info('Phase4+ cache restore for %s: %s', genome, genome_restored)

            # Skip Stage 2 COMPUTE only when this genome is already at the
            # current PIPELINE_VERSION AND every SELECTED step was actually
            # restored from cache above — then there is nothing left to compute.
            # If a selected step had a cache miss (e.g. tmbed newly enabled and
            # not yet cached for this organism), fall through to Stage 2 so
            # Snakemake computes just the missing selected steps; the cache-hit
            # files restored above are skipped by their fresh mtime.
            if db_path:
                fasta_hash = compute_fasta_hash(genome_files[genome])
                missing_selected = [
                    t for t, hit in genome_restored.items()
                    if not hit and smk_config.get(f'run_{t}', True) not in (False, 'false', '0', 'no')
                ]
                # Scoring is never cached (see _genome_cache_map): the OCC operon
                # reference grows over time, so a genome must be RE-SCORED on every
                # run even when every cached step was restored. Never take the
                # already-processed fast-path while scoring is selected.
                scoring_always_recomputes = smk_config.get('run_scoring', True) not in (False, 'false', '0', 'no')
                if (is_already_processed(db_path, fasta_hash, PIPELINE_VERSION)
                        and not missing_selected and not scoring_always_recomputes):
                    LOGGER.info('Genome %s already at pipeline version %s — skipping Stage 2 compute (all selected outputs restored from cache)',
                                genome, PIPELINE_VERSION)
                    log_workflow_run(db_path, run_id, genome_files[genome], key_name, 0,
                                     status='success')
                    continue
                if missing_selected:
                    LOGGER.info('Genome %s at version %s but selected steps not in cache: %s — running Stage 2 to compute them',
                                genome, PIPELINE_VERSION, missing_selected)
                elif scoring_always_recomputes and is_already_processed(db_path, fasta_hash, PIPELINE_VERSION):
                    LOGGER.info('Genome %s already at version %s, but scoring is never cached — running Stage 2 to RE-SCORE (only scoring + fingerprint recompute; all other phases restore from cache)',
                                genome, PIPELINE_VERSION)

            stage2_config = {**smk_config, 'target_genome': genome}
            stage2_command = self.build_executable(selected_wf, config_overrides=stage2_config, mode=mode,
                                                    compute_config=compute_config, extra_resources=extra_resources,
                                                    rerun_triggers=rerun_triggers, target=stage2_target)
            LOGGER.info('Starting Stage 2 (%s) for %s: %s', stage2_target, genome, ' '.join(stage2_command))
            proc = self._run_subprocess(stage2_command)
            last_proc = proc

            if proc is None or proc.returncode != 0:
                rc = proc.returncode if proc else None
                LOGGER.error('Stage 2 (phase4-10) failed for genome "%s" (rc=%s) -- halting the sequential queue; '
                             'later genomes (even ones already RASTtk-ready) will not be processed this run.',
                             genome, rc)
                result = (self._build_result(key_name, proc) if proc else
                          {'workflow': key_name, 'returncode': rc, 'rules_summary': {}, 'stdout_tail': '', 'stderr_tail': ''})
                if db_path:
                    log_workflow_run(db_path, run_id, genome_files[genome], key_name,
                                     result['rules_summary'].get('completed', 0), status='failed')
                self.failed(msg=f'Workflow "{key_name}" failed on genome "{genome}" (phase4-10)', dex=result)
                return rc or 1

            if db_path:
                genome_cache_map = cache_map_fn(genome)
                genome_restored = restored.get(genome, {})
                tools_to_store = {tool: paths for tool, paths in genome_cache_map.items()
                                  if not genome_restored.get(tool, False)}
                if tools_to_store:
                    store_all(db_path, genome_files[genome], tools_to_store)
                else:
                    LOGGER.info('All outputs were cache hits for %s — skipping redundant storage', genome)
                log_workflow_run(db_path, run_id, genome_files[genome], key_name,
                                 self._build_result(key_name, proc)['rules_summary'].get('completed', 0), status='success')

        # Stage 3: LLM for all genomes at once (only when LLM is enabled and at
        # least one genome reached Stage 2). SLURM queues GPU jobs automatically;
        # --keep-going (already in build_executable) lets skipped genomes' LLM
        # jobs fail without aborting the rest.
        if run_llm_enabled and processed > 0:
            LOGGER.info('=== Stage 3: LLM scoring for all %d genome(s) (rule llm_all) ===', processed)
            stage3_command = self.build_executable(selected_wf, config_overrides=smk_config, mode=mode,
                                                    compute_config=compute_config, extra_resources=extra_resources,
                                                    rerun_triggers=rerun_triggers, target='llm_all')
            LOGGER.info('Starting Stage 3 (llm_all): %s', ' '.join(stage3_command))
            stage3_proc = self._run_subprocess(stage3_command)
            last_proc = stage3_proc
            if stage3_proc is None or stage3_proc.returncode != 0:
                rc_val = stage3_proc.returncode if stage3_proc else None
                LOGGER.error('Stage 3 (LLM) failed (rc=%s) -- pipeline CPU phases completed successfully; '
                             'LLM results are missing but DB is otherwise intact.', rc_val)
                result = (self._build_result(key_name, stage3_proc) if stage3_proc else
                          {'workflow': key_name, 'returncode': rc_val, 'rules_summary': {}, 'stdout_tail': '', 'stderr_tail': ''})
                self.failed(msg=f'Workflow "{key_name}" Stage 3 (LLM) failed', dex=result)
                return rc_val or 1

        stage1.wait()
        if stage1.poll() != 0:
            LOGGER.warning('Stage 1 (rasttk_all) exited with rc=%s after Stage 2/3 finished draining', stage1.poll())

        if skipped:
            LOGGER.warning('%d genome(s) skipped (RASTtk/GTDB-Tk never completed): %s', len(skipped), sorted(skipped))

        run_scoring_enabled = smk_config.get('run_scoring', True) not in (False, 'false', '0', 'no')
        queue_sqlite_backup_enabled = smk_config.get('run_sqlite_snapshot_queue', True) not in (False, 'false', '0', 'no')
        if run_scoring_enabled and queue_sqlite_backup_enabled and processed > 0 and not skipped:
            LOGGER.info('=== FINALIZE: queue background sqlite snapshot backup ===')
            finalize_command = self.build_executable(
                selected_wf,
                config_overrides=smk_config,
                mode=mode,
                compute_config=compute_config,
                extra_resources=extra_resources,
                rerun_triggers=rerun_triggers,
                target='queue_sqlite_backup_snapshot',
            )
            LOGGER.info('Starting finalize target (queue_sqlite_backup_snapshot): %s', ' '.join(finalize_command))
            finalize_proc = self._run_subprocess(finalize_command)
            if finalize_proc is None or finalize_proc.returncode != 0:
                LOGGER.warning(
                    'Could not queue sqlite backup snapshot (target queue_sqlite_backup_snapshot, rc=%s). '
                    'Pipeline outputs are complete; continuing without blocking completion status.',
                    finalize_proc.returncode if finalize_proc else None,
                )

        # Independent, downstream-only pangenome report figures. Reads finished
        # scoring outputs only; runs as an isolated finalize subprocess so a
        # figure failure degrades to a warning and can never block completion.
        report_figures_enabled = smk_config.get('run_report_figures', True) not in (False, 'false', '0', 'no')
        if run_scoring_enabled and report_figures_enabled and processed > 0 and not skipped:
            LOGGER.info('=== FINALIZE: generate pangenome report figures ===')
            figures_command = self.build_executable(
                selected_wf,
                config_overrides=smk_config,
                mode=mode,
                compute_config=compute_config,
                extra_resources=extra_resources,
                rerun_triggers=rerun_triggers,
                target='run_report_figures_global',
            )
            LOGGER.info('Starting finalize target (run_report_figures_global): %s', ' '.join(figures_command))
            figures_proc = self._run_subprocess(figures_command)
            if figures_proc is None or figures_proc.returncode != 0:
                LOGGER.warning(
                    'Could not generate pangenome report figures (target run_report_figures_global, rc=%s). '
                    'Pipeline outputs are complete; continuing without blocking completion status.',
                    figures_proc.returncode if figures_proc else None,
                )

        # FINALIZE: reorganize each organism folder into a clean 3-item layout
        # (FINAL tsv + colored FINAL xlsx + diagrams/ + per-tool-phased-output/).
        # MUST run last -- after the global report above, which reads every
        # organism's scoring/ files off disk. Plain post-step (not a Snakemake
        # rule) so moving outputs never confuses Snakemake's DAG. Non-blocking:
        # a failure degrades to a warning and never blocks completion.
        reorganize_enabled = smk_config.get('run_reorganize_outputs', True) not in (False, 'false', '0', 'no')
        output_dir = (smk_config.get('output_dir') or self.conf.get('output_dir', '') or '').rstrip('/')
        if run_scoring_enabled and reorganize_enabled and processed > 0 and not skipped and output_dir:
            LOGGER.info('=== FINALIZE: reorganize per-organism output folders ===')
            reorg_script = str(WORKFLOW_DIR / 'reorganize_outputs.py')
            excel_script = str(WORKFLOW_DIR / 'fingerprint' / 'make-final-excel.py')
            reorg_workers = str(min(4, max(1, len(genome_files))))
            reorg_command = [
                sys.executable, reorg_script,
                '--run-root', output_dir,
                '--genomes', *genome_files.keys(),
                '--figures-dirname', 'diagrams',
                '--excel-script', excel_script,
                '--python', sys.executable,
                '--workers', reorg_workers,
            ]
            LOGGER.info('Starting finalize step (reorganize_outputs): %s', ' '.join(reorg_command))
            try:
                reorg_proc = subprocess.run(reorg_command, capture_output=True, text=True)
                for line in (reorg_proc.stdout or '').splitlines():
                    LOGGER.info('[reorganize] %s', line)
                if reorg_proc.returncode != 0:
                    LOGGER.warning(
                        'Output reorganization exited rc=%s (stderr: %s). Pipeline outputs are '
                        'complete; continuing without blocking completion status.',
                        reorg_proc.returncode, (reorg_proc.stderr or '').strip()[:500],
                    )
            except Exception as reorg_exc:  # never let cleanup fail the run
                LOGGER.warning('Output reorganization raised %s; continuing.', reorg_exc)

        result = self._build_result(key_name, last_proc) if last_proc else {'workflow': key_name, 'rules_summary': {}}
        self.succeeded(
            msg=f'Workflow "{key_name}" completed ({processed}/{total} genome(s) processed'
                f'{f", {len(skipped)} skipped" if skipped else ""})',
            dex=result,
        )

    @command
    def do_example(self):
        '''example workflow to execute'''
        input_file = self.conf.get('input', None)
        if not input_file:
            LOGGER.error('No input file specified. Use: dane_wf example input: <file>')
            self.failed('No input file specified')
            return 1

        input_path = Path(input_file)
        prodigal_config = self.conf.get('prodigal', {})

        smk_config = {
            'input_fasta': input_file,
            'output_fasta': f"{input_path.stem}-output.txt",
            'prodigal_threads': prodigal_config.get('threads', 4),
        }

        self._run_pipeline('example', smk_config)

    def _selftest_config(self, stem, tmpdir, inject_failure=False):
        '''Build smk_config and cache_map for selftest workflows.'''
        td = Path(tmpdir)
        out_step_a = str(td / f"step_a/{stem}-step_a.out")
        out_step_a_extra = str(td / f"step_a/{stem}-step_a.extra")
        out_step_a_db = str(td / f"step_a/{stem}-step_a_db.tkn")
        out_step_b = str(td / f"step_b/{stem}-step_b.out")
        out_step_b_db = str(td / f"step_b/{stem}-step_b_db.tkn")
        out_step_c_primary = str(td / f"step_c/{stem}-step_c.tsv")
        out_step_c_secondary = str(td / f"step_c/{stem}-step_c_count.tsv")
        out_step_c_db = str(td / f"step_c/{stem}-step_c_db.tkn")

        # For selftest, use temp DB path (not required from config)
        selftest_db = str(td / 'selftest.db')

        smk_config = {
            'workdir': tmpdir,
            'stem': stem,
            'inject_failure': str(inject_failure).lower(),
            'out_step_a': out_step_a,
            'out_step_a_extra': out_step_a_extra,
            'out_step_a_db': out_step_a_db,
            'out_step_b': out_step_b,
            'out_step_b_db': out_step_b_db,
            'out_step_c_primary': out_step_c_primary,
            'out_step_c_secondary': out_step_c_secondary,
            'out_step_c_db': out_step_c_db,
            'main_database': selftest_db,
        }

        cache_map = {
            'step_a': [out_step_a, out_step_a_extra],
            'step_a_db': [out_step_a_db],
            'step_b': [out_step_b],
            'step_b_db': [out_step_b_db],
            'step_c': [out_step_c_primary, out_step_c_secondary],
            'step_c_db': [out_step_c_db],
        }

        return smk_config, cache_map

    @command
    def do_quick_example(self, inject_failure=False):
        '''Run selftest with real margie.db cache (deterministic input — cached on second run).'''
        stem = 'quick-example'

        with tempfile.TemporaryDirectory(prefix='dane_quick_') as tmpdir:
            # Deterministic content so the hash is stable across runs.
            # First run: cache miss → snakemake runs → store_all caches.
            # Second run: cache hit → restore_all writes files → snakemake skips.
            tmp_input = str(Path(tmpdir) / f'{stem}.txt')
            Path(tmp_input).write_text('quick-example deterministic input\n')

            smk_config, cache_map = self._selftest_config(stem, tmpdir, inject_failure)
            smk_config['input_file'] = tmp_input

            self._run_pipeline('selftest', smk_config, cache_map, mode='dev')

    @command
    def do_fresh_test(self, inject_failure=False):
        '''Run selftest with real margie.db — unique input each run so cache always misses.'''
        stem = 'fresh-test'

        with tempfile.TemporaryDirectory(prefix='dane_freshtest_') as tmpdir:
            # Unique content per run (includes timestamp) so the hash is always new.
            # restore_all will miss → snakemake runs all rules → store_all caches.
            tmp_input = str(Path(tmpdir) / f'{stem}.txt')
            Path(tmp_input).write_text(f'fresh-test {self.timestamp}\n')

            smk_config, cache_map = self._selftest_config(stem, tmpdir, inject_failure)
            smk_config['input_file'] = tmp_input

            self._run_pipeline('selftest', smk_config, cache_map, mode='dev')

    @command
    def do_margie(self, mode='slurm'):
        '''run margie workflow'''
        input_file = (self.conf.get('input', None)
                      or self.conf.get('margie', {}).get('input_path', None)
                      or WORKFLOW_PATH_DEFAULTS.get('margie', {}).get('input_path'))
        if not input_file:
            LOGGER.error('No input file specified. Use: dane_wf margie input: <file>, '
                        'or set margie.input_path in your ~/.config/bioinformatics-tools/config.yaml')
            self.failed('No input file specified')
            return 1

        # Require main_database from config - no fallback
        main_database = self.conf.get('main_database', None)
        if not main_database:
            LOGGER.error('main_database not set in config. Add main_database: <path> to your ~/.config/bioinformatics-tools/config.yaml')
            self.failed('main_database configuration is required')
            return 1

        # Expand ~ in database path (SQLite doesn't understand ~)
        main_database = str(Path(main_database).expanduser())

        # Extract and validate compute config for SLURM mode
        compute_config = None
        if mode != 'dev':
            compute_config = self.conf.get('compute', {}).get('cluster_default', {})
            slurm_account = compute_config.get('account', '').strip()
            if not slurm_account:
                LOGGER.error('compute.cluster_default.account not set in config. Add account: <your-slurm-account> to your ~/.config/bioinformatics-tools/config.yaml')
                self.failed('SLURM account configuration is required for cluster execution')
                return 1

        stem = Path(input_file).stem
        output_dir = (self.conf.get('output_dir', None)
                      or self.conf.get('margie', {}).get('output_path', None)
                      or WORKFLOW_PATH_DEFAULTS.get('margie', {}).get('output_path', ''))
        prefix = f"{output_dir.rstrip('/')}/" if output_dir else ''

        # Only pass workflow-specific overrides to Snakemake
        # The full config is loaded via --configfile from self.config_paths
        config_overrides = {
            'input_fasta': input_file,
            'output_dir': prefix.rstrip('/'),
            'main_database': main_database,
        }

        # Cache map - compute paths using same logic as workflow_helpers
        # Note: These paths must match what margie.smk generates (includes stem subdirectory)
        prefix_with_stem = f"{prefix}{stem}/"
        out_prodigal_gff = f"{prefix_with_stem}prodigal/{stem}-prodigal.gff"
        out_prodigal_faa = f"{prefix_with_stem}prodigal/{stem}-prodigal.faa"
        out_prodigal_db = f"{prefix_with_stem}prodigal/prodigal_db.tkn"
        out_pfam = f"{prefix_with_stem}pfam/pfam.tsv"
        out_pfam_db = f"{prefix_with_stem}pfam/pfam_db.tkn"
        out_cog_tkn = f"{prefix_with_stem}cog/cog.tkn"
        out_cog_classify = f"{prefix_with_stem}cog/cog_classify.tsv"
        out_cog_count = f"{prefix_with_stem}cog/cog_count.tsv"
        out_cog_db = f"{prefix_with_stem}cog/cog_db.tkn"
        out_kofam = f"{prefix_with_stem}kofam/kofam.tsv"
        out_kofam_db = f"{prefix_with_stem}kofam/kofam_db.tkn"
        out_uniop = f"{prefix_with_stem}uniop/operons.tsv"
        out_uniop_db = f"{prefix_with_stem}uniop/uniop_db.tkn"
        out_dbcan = f"{prefix_with_stem}dbcan/overview.tsv"
        out_dbcan_db = f"{prefix_with_stem}dbcan/dbcan_db.tkn"

        # Cache map - each tool includes ALL files (intermediates + token)
        # This ensures when we have a cache HIT, we restore everything Snakemake needs
        cache_map = {
            'prodigal': [out_prodigal_gff, out_prodigal_faa, out_prodigal_db],
            'pfam': [out_pfam, out_pfam_db],
            'cog': [out_cog_tkn, out_cog_classify, out_cog_count, out_cog_db],
            'kofam': [out_kofam, out_kofam_db],
            'uniop': [out_uniop, out_uniop_db],
            'dbcan': [out_dbcan, out_dbcan_db],
        }

        self._run_pipeline('margie', config_overrides, cache_map, mode=mode, compute_config=compute_config)

    @command
    def do_margie_sb(self, mode='slurm'):
        '''run margie_sb workflow — input may be a single genome file or a folder of genomes'''
        input_path_value = (self.conf.get('input', None)
                            or self.conf.get('margie_sb', {}).get('input_path', None)
                            or WORKFLOW_PATH_DEFAULTS.get('margie_sb', {}).get('input_path'))
        if not input_path_value:
            LOGGER.error('No input specified. Use: dane_wf margie_sb input: <file_or_folder>, '
                        'or set margie_sb.input_path in your ~/.config/bioinformatics-tools/config.yaml')
            self.failed('No input file or folder specified')
            return 1

        main_database = self.conf.get('main_database', None)
        if not main_database:
            LOGGER.error('main_database not set in config. Add main_database: <path> to your ~/.config/bioinformatics-tools/config.yaml')
            self.failed('main_database configuration is required')
            return 1

        main_database = str(Path(main_database).expanduser())

        compute_config = None
        if mode != 'dev':
            compute_config = self.conf.get('compute', {}).get('cluster_default', {})
            slurm_account = compute_config.get('account', '').strip()
            if not slurm_account:
                LOGGER.error('compute.cluster_default.account not set in config. Add account: <your-slurm-account> to your ~/.config/bioinformatics-tools/config.yaml')
                self.failed('SLURM account configuration is required for cluster execution')
                return 1

        # Recursive: margie_sb's synteny-input/<genome>/... reference genomes
        # (family/genus/order relatives used by the synteny tool) live nested
        # under input_fasta, and also run as primary genomes in their own
        # right, not just get read by synteny.
        genomes = discover_genomes(input_path_value, recursive=True)
        if not genomes:
            LOGGER.error('No genome files found at %s', input_path_value)
            self.failed(f'No genome files found at {input_path_value}')
            return 1

        output_dir = (self.conf.get('output_dir', None)
                      or self.conf.get('margie_sb', {}).get('output_path', None)
                      or WORKFLOW_PATH_DEFAULTS.get('margie_sb', {}).get('output_path', ''))

        config_overrides = {
            'input_fasta': input_path_value,
            'output_dir': output_dir.rstrip('/') if output_dir else '',
            'main_database': main_database,
        }

        # Production default: stop at scoring (phase11); post-scoring phases
        # remain opt-in placeholders and can be explicitly re-enabled.
        _tool_to_run_flag = {'scoring_heuristic': 'run_scoring'}
        _post_scoring_tools = {
            'fingerprint',
            'fingerprint_database',
            'operon_fingerprint',
            'annotation_gff',
            'genome_pool',
            'ani',
            'aai',
            'closest',
            'mauve',
            'synteny',
            'evidence',
            'llm',
        }
        for tool_key in _post_scoring_tools:
            config_overrides[_tool_to_run_flag.get(tool_key, f'run_{tool_key}')] = False

        # margie_sb.selected_tools: comma-joined tool keys, set by the API from
        # the caller's phase selection. Missing/empty means "use production
        # defaults" (stop at scoring). When given, selected tools are enabled
        # and unselected tools are disabled explicitly.
        selected_tools_raw = self.conf.get('margie_sb', {}).get('selected_tools', '')
        sif_files_override = None
        if selected_tools_raw:
            selected_tool_keys = {t.strip() for t in selected_tools_raw.split(',') if t.strip()}
            for tool in MARGIE_SB_PHASED_TOOLS:
                run_flag = _tool_to_run_flag.get(tool['key'], f"run_{tool['key']}")
                config_overrides[run_flag] = tool['key'] in selected_tool_keys
            # gtdbtk/rasttk's container always runs even when deselected --
            # run_gtdbtk=false only skips the DB load (rasttk can't be
            # deselected at all) -- so their SIFs stay validated regardless.
            sif_files_override = margie_sb_sif_files(selected_tool_keys | {'gtdbtk', 'rasttk'})

        # --- Licensing gate ---------------------------------------------------
        # Require accepted terms (interactive first run, or web-app acceptance
        # passed via env), then disable any license-required tool the operator
        # is not entitled to run. The web path already enforced acceptance in
        # dane-api and passes the entitlement down, so this never re-prompts it.
        # Make sure the backend secret keys exist before importing the licensing
        # gate (it pulls in the API layer, which requires them at import time).
        # On a fresh CLI clone this generates + persists them to .env once.
        from bioinformatics_tools.workflow_tools.env_keys import ensure_api_keys
        ensure_api_keys()
        from bioinformatics_tools.workflow_tools.license_gate import (
            ensure_cli_license, LicenseError,
        )
        from bioinformatics_tools.api.licensing.catalog import disabled_tool_ids
        try:
            _entitlement = ensure_cli_license()
        except LicenseError as _exc:
            LOGGER.error('%s', _exc)
            self.failed(str(_exc))
            return 1
        _gateable_keys = {tool['key'] for tool in MARGIE_SB_PHASED_TOOLS}
        _disabled = disabled_tool_ids(
            _entitlement.get('usage_type'), _entitlement.get('licensed_tools')
        ) & _gateable_keys
        if _disabled:
            # Refuse if the caller explicitly asked for a tool they can't run.
            if selected_tools_raw:
                _conflict = selected_tool_keys & _disabled
                if _conflict:
                    _msg = (
                        'These tools require a license you have not confirmed: '
                        f"{', '.join(sorted(_conflict))}. Remove them from "
                        'selected_tools, or record the license via the pipeline '
                        'setup / web app, then re-run.'
                    )
                    LOGGER.error('%s', _msg)
                    self.failed(_msg)
                    return 1
            # Otherwise disable them (and drop their containers from validation).
            for _tid in _disabled:
                config_overrides[_tool_to_run_flag.get(_tid, f'run_{_tid}')] = False
            if selected_tools_raw:
                _keep = (selected_tool_keys - _disabled) | {'gtdbtk', 'rasttk'}
            else:
                _keep = _gateable_keys - _disabled
            sif_files_override = margie_sb_sif_files(_keep)
            LOGGER.warning(
                'Licensing: disabled (not licensed for your usage type): %s',
                ', '.join(sorted(_disabled)),
            )

        def _genome_cache_map(genome: str) -> dict[str, list[str]]:
            '''Every phase4-8 tool's real output files for one genome, keyed by
            tool name, for restore_all()/store_all() -- mirrors reference-work's
            do_margie() cache_map (each tool's results + intermediates + token).

            Excludes quast/gtdbtk/rasttk (phase1-3): quast and gtdbtk each run
            as one batched rule across every genome, gated by a shared
            batch-done marker that won't exist in a fresh output_dir, so
            Snakemake must replan the full batch regardless; rasttk's input is
            gtdbtk's per-genome split output, looping back into the same
            problem.'''
            prefix = get_workflow_prefix_for(genome, config_overrides)
            simple_tools = ['cog', 'pfam', 'merops', 'tcdb', 'uniprot', 'kegg', 'eggnog',
                            'dbcan', 'pgap', 'geneprop', 'operon', 'tmbed', 'signalp6',
                            'deepsig', 'psortb', 'signalp4']
            cache_map = {t: [f'{prefix}{t}/{t}_results.tsv', f'{prefix}{t}/{t}_db.tkn'] for t in simple_tools}
            cache_map['tigrfam'] = [f'{prefix}tigrfam/tigrfam_results.tsv',
                                     f'{prefix}tigrfam/tigrfam_domtbl.out',
                                     f'{prefix}tigrfam/tigrfam_db.tkn']
            cache_map['phobius'] = [f'{prefix}phobius/phobius_results.tsv',
                                     f'{prefix}phobius/phobius_top1.tsv',
                                     f'{prefix}phobius/phobius_db.tkn']
            cache_map['envelope'] = [f'{prefix}envelope/envelope_results.tsv',
                                      f'{prefix}envelope/envelope_summary.tsv',
                                      f'{prefix}envelope/envelope_db.tkn']
            interpro_paths = [f'{prefix}interpro/interpro_results.tsv', f'{prefix}interpro/interpro_db.tkn']
            for db in INTERPRO_DB_BASENAMES:
                interpro_paths += [f'{prefix}interpro/interpro_{db}_results.tsv',
                                    f'{prefix}interpro/interpro_{db}_db.tkn']
            cache_map['interpro'] = interpro_paths
            # Phase 9 (consolidation): include the wide merged TSV so labeling
            # can re-run from cache if it ever has a cache miss while
            # consolidation does not.
            cache_map['consolidation'] = [
                f'{prefix}consolidation/detected-columns.json',
                f'{prefix}consolidation/consolidated-merged-all-columns.tsv',
                f'{prefix}consolidation/consolidated-no-stat.tsv',
                f'{prefix}consolidation/manifest.tsv',
                f'{prefix}consolidation/consolidation_compute.tkn',
                f'{prefix}consolidation/consolidation_db.tkn',
            ]
            # Phase 10 (labeling)
            cache_map['labeling'] = [
                f'{prefix}labeling/labeled-genes.tsv',
                f'{prefix}labeling/labeled-genes-ec-consensus.tsv',
                f'{prefix}labeling/labeled-genes-operon-info.tsv',
                f'{prefix}labeling/labeled-genes-cluster-agreement.tsv',
                f'{prefix}labeling/labeling_compute.tkn',
                f'{prefix}labeling/labeling_db.tkn',
            ]
            # Phase 11 (scoring) is intentionally NOT cached. The C3 Operon
            # Context Confidence factor scores each gene against a cross-organism
            # OCC operon reference that GROWS as new organisms are labeled, so the
            # same genome can score differently over time. Caching scoring would
            # restore stale scores; instead scoring recomputes on EVERY run (the
            # Stage-2 already-processed fast-path below is disabled whenever
            # run_scoring is selected) and its results overwrite the DB each run
            # (rule load_scoring_to_db passes --force). A timestamped snapshot of
            # every run's final scores is archived to depot (rule
            # archive_scoring_to_depot) so the history is never lost. Fingerprint
            # (phase 12) is likewise not cached here (it never was).
            return cache_map

        # RASTtk concurrency is enforced by a mkdir-based mutex in
        # margie_sb.smk's run_rasttk rule, not a Snakemake resource pool here
        # (the pool leaked under SLURM's jobstep executor) -- nothing to register.

        # How many phase4 tools may run at once, independent of how many
        # genomes/other phases run concurrently (margie_sb.phase4.max_parallel_tools).
        max_parallel_tools = self.conf.get('margie_sb', {}).get('phase4', {}).get('max_parallel_tools', 4)
        extra_resources = {
            'margie_sb_phase4_slot': max_parallel_tools,
        }

        # margie_sb.resume: set only by /v1/ssh/resume_job's relaunch, after it
        # has copied a failed run's output_dir forward into this one's.
        # mtime-only rerun triggers let Snakemake recognize those
        # copied-forward outputs as done despite output_dir having changed
        # (combined triggers, the default, would treat that as a params
        # change and redo everything).
        resume_raw = str(self.conf.get('margie_sb', {}).get('resume', '')).strip().lower()
        rerun_triggers = 'mtime' if resume_raw not in ('false', '0', 'no', 'off', '') else None

        self._run_pipeline_batch_sequential('margie_sb', config_overrides, genomes,
                                 cache_map_fn=_genome_cache_map, mode=mode, compute_config=compute_config,
                                 extra_resources=extra_resources, sif_files_override=sif_files_override,
                                 rerun_triggers=rerun_triggers)
