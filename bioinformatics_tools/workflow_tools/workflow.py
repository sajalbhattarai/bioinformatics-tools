'''
Workflow tools generate
Invoked: $ dane_wf wf: example <params/options/io>
'''
import hashlib
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
from bioinformatics_tools.workflow_tools.output_cache import (
    cached_tools, current_rast_genome_id, log_workflow_run, realign_rast_genome_id,
    restore, restore_all, store, store_all,
)
from bioinformatics_tools.workflow_tools.load_to_db import is_already_processed, PIPELINE_VERSION, fasta_hashes
from bioinformatics_tools.workflow_tools import protein_cache
from bioinformatics_tools.workflow_tools.programs import ProgramBase
from bioinformatics_tools.workflow_tools.workflow_helpers import (
    discover_genomes, get_workflow_prefix_for, WORKFLOW_PATH_DEFAULTS, genome_calls)
from bioinformatics_tools.workflow_tools.workflow_registry import (
    MARGIE_SB_PHASED_TOOLS,
    WORKFLOWS,
    margie_sb_sif_files,
)

LOGGER = logging.getLogger(__name__)


_SCRIPT_VERSIONS: dict[str, str] = {}


def _scripts_versioned(step: str) -> str:
    """Returns '<step>@<fingerprint of its scripts>', the cache name for steps
    that run MARGIE's own code (consolidation, labeling), so changed code is a
    cache miss while the tools' cached outputs stay valid."""
    if step not in _SCRIPT_VERSIONS:
        folder = Path(__file__).parent / step
        digest = hashlib.sha256()
        for f in sorted(folder.glob('*')):
            if f.is_file() and f.suffix in ('.py', '.sh', '.json', '.tsv'):
                digest.update(f.name.encode() + b'\0' + f.read_bytes())
        _SCRIPT_VERSIONS[step] = f'{step}@{digest.hexdigest()[:12]}'
    return _SCRIPT_VERSIONS[step]

# Mirrors margie_sb.smk's INTERPRO_ANALYSIS_TO_BASENAME.values(); the .smk is not importable.
INTERPRO_DB_BASENAMES = [
    "antifam", "cdd", "coils", "funfam", "gene3d", "hamap", "mobidb", "ncbifam",
    "panther", "pfam", "pirsf", "pirsr", "prints", "prosite_patterns",
    "prosite_profiles", "sfld", "smart", "superfamily",
]
WORKFLOW_DIR = Path(__file__).parent

# --cores for SLURM-mode runs. Not 'all': Snakemake clamps every rule's threads:
# to --cores, and inside the small driver allocation (ssh_slurm.DRIVER_CPUS)
# 'all' would shrink e.g. gtdbtk's 64 threads to 4. 256 is the largest node;
# local rules are still bounded by --local-cores.
SLURM_MODE_CORES = 256


def _subprocess_env() -> dict:
    '''Returns the env for every snakemake subprocess, with BASH_ENV set to Lmod's
    init script so the envmodules rules (signalp4/signalp6) find `module`.'''
    env = os.environ.copy()
    env['BASH_ENV'] = '/etc/profile.d/modules.sh'
    return env


def _snakemake_executable() -> str:
    '''Returns the snakemake binary beside the running interpreter, since dane_wf
    may run from a venv whose bin/ is not on PATH.
    '''
    candidate = Path(sys.executable).parent / 'snakemake'
    return str(candidate) if candidate.exists() else 'snakemake'


class _BackgroundProcess:
    '''Handle returned by WorkflowBase._start_background_subprocess(): a Popen
    whose output is drained on background threads, so .poll()/.wait() never block.'''

    def __init__(self, proc: subprocess.Popen, stdout_lines: list[str], stderr_lines: list[str]):
        self.proc = proc
        self.stdout_lines = stdout_lines
        self.stderr_lines = stderr_lines

    def poll(self):
        '''Returns None while running, else the process's exit code.'''
        return self.proc.poll()

    def wait(self) -> int:
        return self.proc.wait()


class _Stage1Waves:
    '''Stage 1 as up to three Snakemake invocations behind one
    _BackgroundProcess-shaped handle.

    GTDB-Tk runs as one batch while output_cache restores it per genome, so:

      batch    GTDB-Tk, only when some genome still needs it
      ready    RASTtk for genomes whose GTDB-Tk outputs are already on disk
      pending  RASTtk for the rest, started once the batch has landed

    "ready" runs alongside "batch"; the RASTtk waves never overlap, since
    BV-BRC must see one submission at a time.
    '''

    def __init__(self, launch, batch_targets, ready_targets, pending_targets):
        self._launch = launch          # (targets, label) -> _BackgroundProcess | None
        self._lock = threading.Lock()
        self._procs: list = []
        self._launch_failed = False
        self._done = threading.Event()
        self._thread = threading.Thread(
            target=self._drive, args=(batch_targets, ready_targets, pending_targets), daemon=True)
        self._thread.start()

    def _start(self, targets, label):
        if not targets:
            return None
        proc = self._launch(targets, label)
        if proc is None:
            self._launch_failed = True
            return None
        with self._lock:
            self._procs.append(proc)
        return proc

    def _drive(self, batch_targets, ready_targets, pending_targets):
        try:
            batch = self._start(batch_targets, 'GTDB-Tk batch')
            ready = self._start(ready_targets, 'RASTtk, GTDB-Tk already available')
            for proc in (batch, ready):
                if proc is not None:
                    proc.wait()
            # Only after the batch has landed and the first RASTtk wave has exited.
            rest = self._start(pending_targets, 'RASTtk, after GTDB-Tk batch')
            if rest is not None:
                rest.wait()
        finally:
            self._done.set()

    def poll(self):
        '''Returns None until every wave has finished, then the first non-zero exit code.'''
        if not self._done.is_set():
            return None
        if self._launch_failed:
            return 1
        with self._lock:
            return next((rc for p in self._procs if (rc := p.poll())), 0)

    def wait(self) -> int:
        self._done.wait()
        return self.poll()


class WorkflowBase(ProgramBase):
    '''Snakemake workflow execution. Inherits single-program commands from ProgramBase.
    '''

    def __init__(self, workflow_id=None):
        LOGGER.debug('Starting __init__ of WorkflowBase')
        self.workflow_id = workflow_id
        self.timestamp = datetime.now().strftime("%d%m%y-%H%M")

        LOGGER.debug('Using the workflow id of %s', self.workflow_id)

        super().__init__()

    def build_executable(self, key: WorkflowKey, config_overrides: dict = None, mode='slurm', compute_config: dict = None,
                          extra_resources: dict = None, rerun_triggers: str = None, target: str = None,
                          max_jobs_override: int = None) -> list[str]:
        '''
        Build snakemake command from workflow key and config.

        Args:
            key: WorkflowKey defining the workflow
            config_overrides: Only workflow-specific overrides (input_fasta, output_dir, main_database)
            mode: Execution mode ('dev' for local-only; anything else uses slurm executor)
            compute_config: Compute cluster config (account, partition, resources)
            extra_resources: Optional named Snakemake resources (--resources k=v ...)
                capping a group of rules independently of --jobs, e.g.
                {'margie_sb_phase4_slot': 4}.
            rerun_triggers: Optional --rerun-triggers value; resume uses 'mtime'
                because a resumed run's output_dir has a new timestamp.
            target: Optional Snakemake target rule (e.g. 'rasttk_all') or list of
                file targets, placed before --config. Adds --nolock, since Stage 1
                and Stage 2 run concurrently on disjoint rules in one directory.
            max_jobs_override: --jobs value that wins over compute_config's
                max_jobs; Stage 1 uses 1 so BV-BRC sees one submission at a time.
        '''
        smk_path = WORKFLOW_DIR / key.snakemake_file

        # Use compute config to determine max_jobs (default to 5)
        max_jobs = 5
        if compute_config:
            max_jobs = compute_config.get('max_jobs', 5)
        if max_jobs_override is not None:
            max_jobs = max_jobs_override

        # 'all' only in dev mode, where every rule runs on this machine (see SLURM_MODE_CORES).
        cores = 'all' if mode == 'dev' else SLURM_MODE_CORES

        core_command = [
            _snakemake_executable(),
            '-s', str(smk_path),
            f'--cores={cores}',
            '--keep-going',
            '--use-apptainer',
            '--sdm=apptainer',
            '--apptainer-args', '-B /home/ddeemer -B /depot/lindems/data/Databases/',  #TODO: HARDCODED!
            # Needed by signalp4/signalp6, the only rules using envmodules: (see margie_sb.smk).
            '--use-envmodules',
            # Prints each job's rule/wildcards/jobid block to the live log, where
            # run_ssh_task reads the job's genome.
            '--verbose',
            f'--jobs={max_jobs}',
            '--latency-wait=60',
            '--scheduler=greedy',
        ]

        if target:
            core_command.append('--nolock')
            # A list means explicit file targets; Stage 1 uses it for a subset of
            # genomes (see _Stage1Waves).
            core_command.extend([target] if isinstance(target, str) else list(target))

        # Check for dry-run/test mode (from config or command line)
        if self.conf.get('dry_run', False) or self.conf.get('test_only', False):
            core_command.append('--dry-run')
            LOGGER.info("Running in DRY-RUN mode - no actual execution")

        if mode != 'dev':
            core_command.append('--executor=slurm')

        # Add default SLURM resources from compute config
        if mode != 'dev' and compute_config:
            default_resources = ['--default-resources']
            min_runtime_minutes = 240

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
                runtime_value = compute_config.get('default_runtime')
                try:
                    runtime_minutes = int(runtime_value)
                except (TypeError, ValueError):
                    runtime_minutes = min_runtime_minutes
                runtime_minutes = max(min_runtime_minutes, runtime_minutes)
                default_resources.append(f'runtime={runtime_minutes}')

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

        # Named resources cap a group of rules independently of --jobs.
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
        '''Starts wf_command like _run_subprocess() (cwd output_dir, output logged
        as "[snakemake] ...") but returns at once, so Stage 1 can run in the
        background. Returns None on launch failure.'''
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

    def _run_pipeline(self, key_name: str, smk_config: dict, cache_map: dict = None, mode='slurm', compute_config: dict = None):
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
                # Never contacts the registry; reports what is already on disk.
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
                             cache_tool: str, cache_paths_fn, mode='slurm', compute_config: dict = None,
                             extra_resources: dict = None, sif_files_override: list[tuple] = None,
                             rerun_triggers: str = None):
        '''Batch sibling of _run_pipeline() for workflows that accept a folder of
        genomes: one snakemake subprocess over the {genome} wildcard, with the
        output cache checked and updated per genome.

        cache_paths_fn(stem) returns that genome's output files to cache.
        sif_files_override replaces selected_wf.sif_files for this call only.
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

        # Per-genome cache restore (SQLite check, before Snakemake runs)
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

        # Success: stores outputs for cache misses and logs every genome's run.
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
        return [f"{prefix}gtdbtk/gtdbtk_results.tsv", f"{prefix}gtdbtk/translation_table.tsv",
                f"{prefix}gtdbtk/gtdbtk_db.tkn"]

    @staticmethod
    def _rasttk_paths(genome: str, smk_config: dict) -> list[str]:
        prefix = get_workflow_prefix_for(genome, smk_config)
        return [f"{prefix}rasttk/rast.tsv", f"{prefix}rasttk/rast.faa", f"{prefix}rasttk/rast.gff",
                f"{prefix}rasttk/rasttk_db.tkn"]

    @staticmethod
    def _restore_gtdbtk_batch(db_path: str, genome_files: dict[str, str], smk_config: dict) -> bool:
        '''All-or-nothing cache restore for GTDB-Tk's batch rule (run_gtdbtk_batch).

        If every genome's split files are cached, restores them and rebuilds the
        batch's combined files from those rows; one miss aborts, and the real
        batch rerun overwrites the partial restore.
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
            results_path, translation_path, _token_path = WorkflowBase._gtdbtk_split_paths(genome, smk_config)
            r_header, *r_rows = Path(results_path).read_text().splitlines()
            t_header, *t_rows = Path(translation_path).read_text().splitlines()
            result_header, translation_header = r_header, t_header
            # Column 0 is set to this loop's genome: identical FASTAs under
            # different names share one cached row.
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

        # Bumps each split file's mtime past the batch files just written, so
        # Snakemake does not rerun the split and, in turn, RASTtk.
        for genome in genome_files:
            results_path, translation_path, _token_path = WorkflowBase._gtdbtk_split_paths(genome, smk_config)
            Path(results_path).touch()
            Path(translation_path).touch()

        LOGGER.info('GTDB-Tk batch cache HIT for all %d genomes — skipping the container run entirely', len(genome_files))
        # Logged only once the whole batch is confirmed a hit.
        for genome in genome_files:
            LOGGER.info("Cache HIT for gtdbtk (genome=%s) — skipping recomputation", genome)
        return True

    def _genome_stage2_progress(self, genome: str, genome_file: str, smk_config: dict) -> tuple[int, int]:
        '''Returns (max_phase_reached, n_phases_done) for a genome, used to run the
        most complete genomes first in Stage 2.

        A phase counts as done if its output cache row exists or its output
        folder is populated. Never raises; an unreadable genome scores (0, 0).'''
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
                                        cache_map_fn, mode='slurm', compute_config: dict = None,
                                        extra_resources: dict = None, sif_files_override: list[tuple] = None,
                                        rerun_triggers: str = None, prodigal_genomes: set[str] = frozenset()):
        '''Sequential-per-organism sibling of _run_pipeline_batch() for margie_sb.

        Stage 1 (rule rasttk_all) runs GTDB-Tk/RASTtk breadth-first in the
        background with --jobs=1, since BV-BRC takes one submission at a time.
        Stage 2 (rule phase4_12_one_genome) runs phases 4-11 for one genome at a
        time, in RASTtk-ready order, choosing the most complete ready genome
        first (see _genome_stage2_progress). A Stage 2 failure halts the queue.
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
        # Off: nothing reads GTDB-Tk (GENOME_INFO comes from margie_sb.genome_info),
        # so its batch is neither restored nor run.
        run_gtdbtk_enabled = smk_config.get('run_gtdbtk', True) not in (False, 'false', '0', 'no')

        restored: dict[str, dict[str, bool]] = {}
        rasttk_restored: dict[str, bool] = {}
        gtdbtk_batch_hit = False
        # Genomes whose GTDB-Tk outputs are committed to output_cache this run,
        # cached as soon as they exist rather than when RASTtk finishes.
        gtdbtk_stored: set[str] = set()
        if db_path:
            # Phase4+ restore is deferred to just before each genome's Stage 2,
            # so restored files are newer than the rast.faa Stage 1 writes.
            LOGGER.info('Phase4+ cache restore deferred to per-genome pre-Stage-2 (after rasttk outputs are fresh)')

            # gtdbtk/rasttk are excluded from cache_map_fn (see _genome_cache_map);
            # _restore_gtdbtk_batch handles the GTDB-Tk batch shape directly.
            gtdbtk_batch_hit = run_gtdbtk_enabled and self._restore_gtdbtk_batch(db_path, genome_files, smk_config)
            if gtdbtk_batch_hit:
                # Only worth attempting after a full gtdbtk batch hit; otherwise
                # gtdbtk reruns and forces rasttk to rerun too.
                for stem, genome_file in genome_files.items():
                    if stem in prodigal_genomes:
                        continue  # Prodigal's genes, not RASTtk's: never from the rasttk cache
                    hit = restore(db_path, genome_file, 'rasttk', self._rasttk_paths(stem, smk_config))
                    rasttk_restored[stem] = hit
                    if hit:
                        LOGGER.info("Cache HIT for rasttk (genome=%s) — skipping recomputation", stem)
            LOGGER.info('GTDB-Tk batch cache hit: %s. RASTtk per-genome cache restore results: %s',
                        gtdbtk_batch_hit, rasttk_restored)

        # Stage 1: phases 1-3 for every genome, in the background, split into
        # waves (see _Stage1Waves). Every wave uses --jobs=1 and the RASTtk waves
        # are sequenced, so BV-BRC sees one submission at a time.
        def _stage1_targets(genome: str) -> list[str]:
            '''Returns the files rule rasttk_all asks for, for one genome: the
            RASTtk DB token, plus the GTDB-Tk one when GTDB-Tk is selected.'''
            prefix = get_workflow_prefix_for(genome, smk_config)
            targets = [f"{prefix}rasttk/rasttk_db.tkn"]
            if run_gtdbtk_enabled:
                targets.append(f"{prefix}gtdbtk/gtdbtk_db.tkn")
            return targets

        # A genome is "ready" when its per-genome GTDB-Tk outputs are on disk
        # (restored above), or always when GTDB-Tk is off.
        gtdbtk_ready = [g for g in genome_files
                        if not run_gtdbtk_enabled
                        or all(Path(p).exists() for p in self._gtdbtk_split_paths(g, smk_config))]
        gtdbtk_pending = [g for g in genome_files if g not in set(gtdbtk_ready)]
        LOGGER.info('Stage 1 waves: %d genome(s) have GTDB-Tk outputs already, %d still need the batch',
                    len(gtdbtk_ready), len(gtdbtk_pending))

        output_dir_for_batch = (smk_config.get('output_dir') or '').rstrip('/')
        batch_targets = ([f"{output_dir_for_batch}/original_container_outputs/gtdbtk/gtdbtk_batch.done"]
                         if gtdbtk_pending and output_dir_for_batch else [])

        def _launch_wave(targets: list[str], label: str):
            command = self.build_executable(selected_wf, config_overrides=smk_config, mode=mode,
                                            compute_config=compute_config, extra_resources=extra_resources,
                                            rerun_triggers=rerun_triggers, target=targets,
                                            max_jobs_override=1)
            LOGGER.info('Starting Stage 1 wave [%s] over %d target(s): %s',
                        label, len(targets), ' '.join(command))
            proc = self._start_background_subprocess(command)
            if proc is None:
                LOGGER.error('Stage 1 wave [%s] failed to launch', label)
            return proc

        stage1 = _Stage1Waves(
            _launch_wave,
            batch_targets,
            [t for g in gtdbtk_ready for t in _stage1_targets(g)],
            [t for g in gtdbtk_pending for t in _stage1_targets(g)],
        )

        def _rasttk_token_path(genome: str) -> str:
            return f"{get_workflow_prefix_for(genome, smk_config)}rasttk/rasttk_db.tkn"

        # With LLM enabled, Stage 2 stops before the GPU step and Stage 3
        # (rule llm_all) submits all LLM jobs at once; gres=gpu:1 serialises them.
        _run_llm_val = smk_config.get('run_llm', False)
        run_llm_enabled = _run_llm_val not in (False, 'false', '0', 'no')
        stage2_target = 'phase4_12_one_genome_no_llm' if run_llm_enabled else 'phase4_12_one_genome'
        LOGGER.info('LLM enabled: %s  →  Stage 2 target: %s', run_llm_enabled, stage2_target)

        # What the per-protein cache reads: the account's config and this run's settings.
        pc_cfg = {**(self.conf or {}), **smk_config}
        pc_on = bool(db_path) and key_name == 'margie_sb' and protein_cache.enabled(pc_cfg)

        total = len(genome_files)
        pending = set(genome_files.keys())
        queue: list[str] = []
        skipped: list[str] = []
        processed = 0
        last_proc = None
        progress_memo: dict[str, tuple[int, int]] = {}  # genome -> completeness, for most-complete-first ordering

        while pending or queue:
            # Caches each genome's GTDB-Tk outputs as soon as its own db token
            # (gtdbtk_db.tkn) exists, so a later RASTtk failure keeps them cached.
            if db_path and run_gtdbtk_enabled and not gtdbtk_batch_hit:
                for genome in genome_files:
                    if genome in gtdbtk_stored:
                        continue
                    if all(Path(p).exists() for p in self._gtdbtk_split_paths(genome, smk_config)):
                        store(db_path, genome_files[genome], 'gtdbtk', self._gtdbtk_split_paths(genome, smk_config))
                        gtdbtk_stored.add(genome)

            for genome in list(pending):
                if Path(_rasttk_token_path(genome)).exists():
                    pending.discard(genome)
                    queue.append(genome)
                    if db_path:
                        if run_gtdbtk_enabled and not gtdbtk_batch_hit and genome not in gtdbtk_stored:
                            store(db_path, genome_files[genome], 'gtdbtk', self._gtdbtk_split_paths(genome, smk_config))
                            gtdbtk_stored.add(genome)
                        if genome not in prodigal_genomes and not rasttk_restored.get(genome, False):
                            store(db_path, genome_files[genome], 'rasttk', self._rasttk_paths(genome, smk_config))

            if not queue:
                if stage1.poll() is not None and pending:
                    # Stage 1 exited without a token for these genomes: a phase1-3
                    # failure for them, skipped here (--keep-going).
                    LOGGER.warning('Stage 1 exited without producing a RASTtk token for: %s -- skipping',
                                   sorted(pending))
                    skipped.extend(pending)
                    pending.clear()
                else:
                    time.sleep(15)
                continue

            # Most complete first among RASTtk-ready genomes; progress is memoised
            # since pending genomes do not change, and the stable sort keeps FIFO
            # order for ties.
            if len(queue) > 1:
                queue.sort(
                    key=lambda g: progress_memo.setdefault(
                        g, self._genome_stage2_progress(g, genome_files[g], smk_config)),
                    reverse=True,
                )
            genome = queue.pop(0)
            processed += 1
            LOGGER.info('=== SEQUENTIAL: genome %d/%d (%s) phase4-12 starting ===', processed, total, genome)

            # Restores phase4+ cached outputs now, after rasttk produced rast.faa,
            # so their mtime is newer and Snakemake skips them. Must precede the
            # already-processed skip below, or a fresh output_dir misses these files.
            if db_path:
                genome_cache_map_now = cache_map_fn(genome)
                genome_restored = restore_all(db_path, genome_files[genome], genome_cache_map_now)
                restored[genome] = genome_restored
                LOGGER.info('Phase4+ cache restore for %s: %s', genome, genome_restored)

                # Rewrites restored tables into this run's RASTtk namespace
                # (6666666.<job>), so they join the freshly computed rasttk output
                # (see realign_rast_genome_id).
                rast_ref = f"{get_workflow_prefix_for(genome, smk_config)}rasttk/rast.gff"
                current_id = current_rast_genome_id(rast_ref)
                if current_id:
                    hit_paths = [p for tool, hit in genome_restored.items() if hit
                                 for p in genome_cache_map_now.get(tool, [])]
                    realigned = realign_rast_genome_id(hit_paths, current_id)
                    if realigned:
                        LOGGER.warning(
                            'Realigned %d cache-restored file(s) for %s onto RASTtk genome id %s '
                            '(cached under a different submission): %s',
                            len(realigned), genome, current_id, sorted(realigned.values())[:1])
                else:
                    LOGGER.warning('No single RASTtk genome id readable from %s — skipping '
                                   'cache realignment for %s', rast_ref, genome)

            # Per-protein cache (protein_cache.py): for each selected tool not
            # restored whole, known proteins come from the cache and only the
            # rest go to the tool; every genome's proteins are stored below.
            pc_prefix = get_workflow_prefix_for(genome, smk_config)
            pc_faa_path = f"{pc_prefix}rasttk/rast.faa"
            pc_tools: list[str] = []
            pc_stored: set[str] = set()
            if pc_on and Path(pc_faa_path).exists():
                pc_tools = [t for t in protein_cache.TOOLS
                            if smk_config.get(f'run_{t}', True) not in (False, 'false', '0', 'no')]
                to_split = [t for t in pc_tools if not genome_restored.get(t, False)]
                for t in pc_tools:
                    if t not in to_split:
                        protein_cache.clear(pc_prefix, t)
                try:
                    counts = protein_cache.split(db_path, pc_faa_path, pc_prefix, genome, to_split, pc_cfg)
                    served = {t: c for t, c in counts.items() if c[0]}
                    LOGGER.info('Protein cache for %s: %s', genome,
                                ', '.join(f'{t} {c}/{c + n} cached' for t, (c, n) in counts.items()) or 'no tool to split')
                    if served:
                        LOGGER.info('Protein cache: %d tool(s) run on fewer proteins for %s', len(served), genome)
                except Exception as exc:
                    LOGGER.warning('Protein cache split failed for %s (tools run on every protein): %s', genome, exc)
                    for t in to_split:
                        protein_cache.clear(pc_prefix, t)

            # Skips Stage 2 only when the genome is at PIPELINE_VERSION and every
            # selected step was restored; otherwise Snakemake computes the missing steps.
            if db_path:
                fasta_hash = fasta_hashes(genome_files[genome])
                missing_selected = [
                    t for t, hit in genome_restored.items()
                    if not hit and smk_config.get(f"run_{t.split('@')[0]}", True) not in (False, 'false', '0', 'no')
                ]
                # Scoring is never cached (its OCC reference grows), so the
                # fast path is never taken while scoring is selected.
                scoring_always_recomputes = smk_config.get('run_scoring', True) not in (False, 'false', '0', 'no')
                # A Prodigal genome restores nothing (see _genome_cache_map), so it always computes.
                if (genome not in prodigal_genomes
                        and is_already_processed(db_path, fasta_hash, PIPELINE_VERSION)
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
            stage2_proc = self._start_background_subprocess(stage2_command)

            # Per-tool cache store, polled while Stage 2 runs: a tool is cached
            # once all its outputs and its db token (written last by
            # load_<tool>_to_db) exist, so another tool's failure cannot discard it.
            genome_cache_map = cache_map_fn(genome)
            genome_restored = restored.get(genome, {})
            genome_stored: set[str] = set()

            def _store_ready_tools():
                if not db_path:
                    return
                # Each tool's proteins, once its results are loaded (its db token).
                for t in pc_tools:
                    if t in pc_stored or not Path(f"{pc_prefix}{t}/{t}_db.tkn").exists():
                        continue
                    pc_stored.add(t)
                    try:
                        added = protein_cache.store(db_path, pc_faa_path, pc_prefix, t, pc_cfg)
                        if added:
                            LOGGER.info('Protein cache: %d new protein(s) of %s kept for %s', added, genome, t)
                    except Exception as exc:
                        LOGGER.warning('Protein cache store failed for %s / %s: %s', genome, t, exc)
                for tool, paths in genome_cache_map.items():
                    if tool in genome_stored or genome_restored.get(tool, False):
                        continue
                    if all(Path(p).exists() for p in paths):
                        store(db_path, genome_files[genome], tool, paths)
                        genome_stored.add(tool)
                        LOGGER.info('Cached %s for %s as soon as its own db token appeared '
                                    '(Stage 2 still running)', tool, genome)

            if stage2_proc is None:
                proc = None
            else:
                while stage2_proc.poll() is None:
                    _store_ready_tools()
                    time.sleep(10)
                _store_ready_tools()  # final catch-up for tokens written right before exit
                proc = subprocess.CompletedProcess(
                    args=stage2_command, returncode=stage2_proc.proc.returncode,
                    stdout='\n'.join(stage2_proc.stdout_lines), stderr='\n'.join(stage2_proc.stderr_lines))
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
                tools_to_store = {tool: paths for tool, paths in genome_cache_map.items()
                                  if tool not in genome_stored and not genome_restored.get(tool, False)}
                if tools_to_store:
                    store_all(db_path, genome_files[genome], tools_to_store)
                else:
                    LOGGER.info('All outputs were already cached incrementally or were cache hits for %s', genome)
                log_workflow_run(db_path, run_id, genome_files[genome], key_name,
                                 self._build_result(key_name, proc)['rules_summary'].get('completed', 0), status='success')

        # Stage 3: LLM for all genomes at once, when enabled and at least one
        # genome reached Stage 2; --keep-going lets one genome's failure pass.
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

        # Pangenome report figures from finished scoring outputs, as an isolated
        # subprocess; a failure only logs a warning.
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

        # FINALIZE: reorganises each organism folder (FINAL tsv, coloured xlsx,
        # diagrams/, per-tool-phased-output/). Runs last, outside Snakemake's DAG;
        # a failure only logs a warning.
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
        stage1_rc = stage1.poll()
        result['stage1_returncode'] = stage1_rc
        result['genomes_total'] = total
        result['genomes_processed'] = processed
        result['genomes_skipped'] = sorted(skipped)

        summary = (f'{processed}/{total} genome(s) processed'
                   f'{f", {len(skipped)} skipped" if skipped else ""}')

        # No genome produced output: reported as a failure, since the GUI reads
        # the exit code to decide completed vs failed.
        if total > 0 and processed == 0:
            self.failed(
                msg=f'Workflow "{key_name}" produced nothing ({summary})'
                    f'{f" -- Stage 1 (rasttk_all) exited rc={stage1_rc}" if stage1_rc else ""}',
                dex=result,
            )
            return stage1_rc or 1

        # Some genomes finished: their outputs stand, but a partial run is
        # reported as inconclusive (exit code 0, no "Success" banner).
        if skipped or stage1_rc:
            self.finished(
                msg=f'Workflow "{key_name}" completed with gaps ({summary})'
                    f'{f" -- Stage 1 (rasttk_all) exited rc={stage1_rc}" if stage1_rc else ""}',
                dex=result,
            )
            return

        self.succeeded(
            msg=f'Workflow "{key_name}" completed ({summary})',
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
        '''Runs the margie_sb workflow on a single genome file or a folder of genomes.'''
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

        # Recursive: synteny-input/<genome>/... reference genomes also run as primary genomes.
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

        # Production default: stops at scoring (phase11); later phases are opt-in.
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

        # margie_sb.selected_tools (set by the API per run) or else
        # margie_sb.default_selected_tools picks the tools; with neither, the
        # production defaults apply (stop at scoring).
        _margie_sb_conf = self.conf.get('margie_sb', {})
        selected_tools_raw = (
            _margie_sb_conf.get('selected_tools', '')
            or _margie_sb_conf.get('default_selected_tools', '')
        )
        sif_files_override = None
        if selected_tools_raw:
            selected_tool_keys = {t.strip() for t in selected_tools_raw.split(',') if t.strip()}
            for tool in MARGIE_SB_PHASED_TOOLS:
                run_flag = _tool_to_run_flag.get(tool['key'], f"run_{tool['key']}")
                config_overrides[run_flag] = tool['key'] in selected_tool_keys
            # rasttk cannot be deselected (it gates phase3); gtdbtk's SIF is needed
            # only when GTDB-Tk runs. Selecting fingerprint_database also enables the
            # operon-fingerprint databases, which grow alongside it.
            if config_overrides.get('run_fingerprint_database'):
                config_overrides['run_operon_fingerprint'] = True
            sif_files_override = margie_sb_sif_files(selected_tool_keys | {'rasttk'} | (
                {'gtdbtk'} if config_overrides.get('run_gtdbtk', True) else set()))

        # ---- Licensing gate ----
        # Requires accepted terms (interactive, or passed via env by the web app)
        # and disables tools the operator is not licensed for. ensure_api_keys()
        # creates the backend secret keys the gate's API imports need.
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
            # Refuses if the caller explicitly asked for a tool they cannot run.
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
            # Otherwise disables them and drops their containers from validation.
            for _tid in _disabled:
                config_overrides[_tool_to_run_flag.get(_tid, f'run_{_tid}')] = False
            if selected_tools_raw:
                _keep = (selected_tool_keys - _disabled) | {'rasttk'} | (
                    {'gtdbtk'} if config_overrides.get('run_gtdbtk', True) else set())
            else:
                _keep = _gateable_keys - _disabled
            sif_files_override = margie_sb_sif_files(_keep)
            LOGGER.warning(
                'Licensing: disabled (not licensed for your usage type): %s',
                ', '.join(sorted(_disabled)),
            )

        # Genomes Prodigal calls instead of RASTtk (GTDB-Tk off and no domain or
        # genetic code in margie_sb.genome_info); they never touch output_cache.
        _calls = genome_calls(genomes, {'run_gtdbtk': config_overrides.get('run_gtdbtk', True),
                                        'margie_sb': self.conf.get('margie_sb', {})})
        prodigal_genomes = {g for g, c in _calls.items() if c['gene_caller'] == 'prodigal'}
        if prodigal_genomes:
            LOGGER.info('Prodigal calls the genes of %d genome(s) with no domain or genetic code: %s',
                        len(prodigal_genomes), ', '.join(sorted(prodigal_genomes)))
            sif_files_override = [*(margie_sb_sif_files() if sif_files_override is None else sif_files_override),
                                  ('prodigal.sif', 'latest')]

        def _genome_cache_map(genome: str) -> dict[str, list[str]]:
            '''Returns each phase4-11 tool's output files for one genome, keyed by
            tool name, for restore_all()/store_all().

            Excludes quast/gtdbtk/rasttk, which run as batches the DAG replans
            anyway. Empty for Prodigal genomes, whose feature ids differ from RASTtk's.'''
            if genome in prodigal_genomes:
                return {}
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
            # Phase 9 (consolidation): includes the wide merged TSV so labeling can rerun from cache.
            cache_map[_scripts_versioned('consolidation')] = [
                f'{prefix}consolidation/detected-columns.json',
                f'{prefix}consolidation/consolidated-merged-all-columns.tsv',
                f'{prefix}consolidation/consolidated-no-stat.tsv',
                f'{prefix}consolidation/manifest.tsv',
                f'{prefix}consolidation/consolidation_compute.tkn',
                f'{prefix}consolidation/consolidation_db.tkn',
            ]
            # Phase 10 (labeling)
            cache_map[_scripts_versioned('labeling')] = [
                f'{prefix}labeling/labeled-genes.tsv',
                f'{prefix}labeling/labeled-genes-ec-consensus.tsv',
                f'{prefix}labeling/labeled-genes-operon-info.tsv',
                f'{prefix}labeling/labeled-genes-cluster-agreement.tsv',
                f'{prefix}labeling/labeling_compute.tkn',
                f'{prefix}labeling/labeling_db.tkn',
            ]
            # Phase 11 (scoring) is not cached: its OCC operon reference grows, so
            # scores are recomputed every run and archived by archive_scoring_to_depot.
            # Phase 12 (fingerprint) is not cached either.
            return cache_map

        # RASTtk concurrency is enforced by a mkdir mutex in margie_sb.smk's run_rasttk rule.

        # Number of phase4 tools running at once (margie_sb.phase4.max_parallel_tools).
        max_parallel_tools = self.conf.get('margie_sb', {}).get('phase4', {}).get('max_parallel_tools', 4)
        extra_resources = {
            'margie_sb_phase4_slot': max_parallel_tools,
        }

        # margie_sb.resume is set by resume_job's relaunch into a copied output_dir;
        # mtime-only rerun triggers then treat the copied outputs as done.
        resume_raw = str(self.conf.get('margie_sb', {}).get('resume', '')).strip().lower()
        rerun_triggers = 'mtime' if resume_raw not in ('false', '0', 'no', 'off', '') else None

        self._run_pipeline_batch_sequential('margie_sb', config_overrides, genomes,
                                 cache_map_fn=_genome_cache_map, mode=mode, compute_config=compute_config,
                                 extra_resources=extra_resources, sif_files_override=sif_files_override,
                                 rerun_triggers=rerun_triggers, prodigal_genomes=prodigal_genomes)
