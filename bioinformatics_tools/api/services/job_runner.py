"""
Job execution and monitoring.

Contains the core SSH task runner that streams remote output, parses logs
for SLURM job IDs, Snakemake progress, and container metadata, and the
SLURM status checker daemon thread.

The `connection` parameter threads a per-user SSHConnection through from the
API router all the way to the SLURM status checker daemon, so every SSH call
hits the correct cluster and account.
"""
import asyncio
import json
import logging
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor

from bioinformatics_tools.utilities import ssh_slurm
from bioinformatics_tools.utilities.ssh_connection import SSHConnection
from bioinformatics_tools.api.services.job_store import job_store

LOGGER = logging.getLogger(__name__)

# Thread pool for background SSH tasks
executor = ThreadPoolExecutor(max_workers=4)

# Regex patterns for parsing Snakemake/SLURM log output
#
# For an ungrouped rule (e.g. "rule_run_quast_batch/"), the rule name is the
# whole "run_quast_batch". For a grouped rule (e.g.
# "group_rasttk_load_rasttk_to_db_run_rasttk/" -- run_<tool> and
# load_<tool>_to_db share one SLURM submission per genome, see margie_sb.smk's
# group: directives), only the group's own short name ("rasttk") is captured
# -- the rest of that path segment is the snakemake-generated concatenation
# of every rule in the group, which read as a confusing "rule name" on its
# own (e.g. showing "load_rasttk_to_db_run_rasttk" made it look like only the
# load step ran, when the run step's actual annotation work happens in the
# very same job).
SLURM_SUBMIT_RE = re.compile(r'SLURM jobid (\d+) \(log: (.*?)\)\.?$')
RULE_NAME_FROM_LOG_PATH_RE = re.compile(r'/slurm_logs/(?:rule_(\w+)|group_([^_]+)_\w+)/')
SLURM_SUBMIT_FALLBACK_RE = re.compile(r'SLURM jobid (\d+)')
STEPS_PROGRESS_RE = re.compile(r'(\d+) of (\d+) steps \((\d+)%\) done')
CACHE_HIT_RE = re.compile(r'Cache HIT for (\w+) \(genome=([^)]+)\)')
# Genome attribution: build_executable() now passes --verbose, so Snakemake
# prints each job's "wildcards: genome=..." line to its own live log right
# before submitting it -- read here into last_genome and attached directly
# at add_slurm_job() time below. Previously genome could only come from
# lazily re-reading each SLURM job's own remote per-job log file (still done
# in _slurm_status_checker as a fallback), which raced against that file
# being cleaned up once the job finished -- intermittently losing genome
# attribution for fast-finishing rules. Safe to track as a single var (not
# per-rule) because margie_sb's sequential per-organism orchestrator (see
# SEQUENTIAL_GENOME_RE below) only ever has one genome in flight at a time,
# even when several rules for that genome are dispatched in the same batch.
WILDCARDS_GENOME_RE = re.compile(r'wildcards:.*\bgenome=([^\s,]+)')
# margie_sb's sequential per-organism orchestrator (workflow.py's
# _run_pipeline_batch_sequential) runs many short-lived Snakemake
# invocations in sequence, one per genome -- each one's own "X of Y steps"
# (STEPS_PROGRESS_RE above) resets relative to just that genome's small
# DAG, which alone would make the frontend's progress bar look like it
# keeps resetting. This is purely additive: genome_index/genome_total are
# new job_store fields, untouched by and not touching steps_done/
# steps_total/progress at all.
SEQUENTIAL_GENOME_RE = re.compile(r'SEQUENTIAL: genome (\d+)/(\d+)')


def _slurm_status_checker(job_id: str, connection: SSHConnection):
    """Daemon thread that periodically checks SLURM job statuses."""
    while job_store.get_status(job_id) == "running":
        slurm_jobs = job_store.get_slurm_jobs(job_id)
        active_ids = [sj["job_id"] for sj in slurm_jobs if sj["status"] not in ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "CACHED")]
        if active_ids:
            try:
                statuses = ssh_slurm.check_multiple_slurm_jobs(active_ids, connection=connection)
                for sj in slurm_jobs:
                    if sj["job_id"] in statuses:
                        sj["status"] = statuses[sj["job_id"]]["state"]
                        sj["time"] = statuses[sj["job_id"]]["time"]
            except Exception as e:
                LOGGER.warning("SLURM status check failed: %s", e)

        # Fallback genome backfill, for any job whose live --verbose
        # "wildcards:" line was missed (see WILDCARDS_GENOME_RE in
        # job_runner.py, the primary source now). Retried every cycle
        # (cheap, one grep each) until it succeeds; permanently empty for
        # batch rules with no genome wildcard.
        for sj in slurm_jobs:
            if not sj.get("genome") and sj.get("log_path"):
                try:
                    genome = ssh_slurm.get_job_genome(sj["log_path"], connection=connection)
                    if genome:
                        sj["genome"] = genome
                except Exception as e:
                    LOGGER.warning("Genome backfill failed for job %s: %s", sj["job_id"], e)
        # Wait 15 seconds between checks
        for _ in range(15):
            if job_store.get_status(job_id) != "running":
                break
            time.sleep(1)


def run_ssh_task(job_id: str, command: str, connection: SSHConnection):
    """Generic SSH task runner with log parsing, SLURM tracking, and progress parsing."""
    job_store.update(job_id, status="running", phase="Submitting via SSH", logs="")

    # Start SLURM status checker daemon thread (passes the same connection through)
    checker = threading.Thread(
        target=_slurm_status_checker,
        args=(job_id, connection),
        daemon=True
    )
    checker.start()

    exit_code = 0  # Track command exit code
    last_genome = ""  # Most recently seen "wildcards: genome=..." value

    try:
        for line in ssh_slurm.submit_ssh_job(cmd=command, connection=connection):
            # Detect exit code metadata from submit_ssh_job
            if line.startswith("__EXIT_CODE__:"):
                try:
                    exit_code = int(line.split(":", 1)[1])
                    LOGGER.info("Captured exit code: %d", exit_code)
                except (ValueError, IndexError):
                    LOGGER.warning("Failed to parse exit code from: %s", line)
                continue

            # Detect structured Report metadata from workflow
            if line.startswith("__REPORT__:"):
                try:
                    report_json = line.split(":", 1)[1]
                    report_data = json.loads(report_json)
                    job_store.update(job_id, report=report_data)
                    LOGGER.info("Captured structured report: status=%s", report_data.get('status', {}).get('code'))
                except (json.JSONDecodeError, IndexError) as e:
                    LOGGER.warning("Failed to parse report JSON: %s", e)
                continue

            # Detect work_dir metadata from submit_ssh_job
            if line.startswith("__WORKDIR__:"):
                job_store.update(job_id, work_dir=line.split(":", 1)[1])
                continue

            # Parse container metadata from bapptainer log lines
            if "__CONTAINER__:" in line:
                try:
                    container_json = line.split("__CONTAINER__:", 1)[1]
                    job_store.add_container(job_id, json.loads(container_json))
                except (json.JSONDecodeError, IndexError):
                    pass

            job_store.append_log(job_id, line)

            # Parse cache-restored rules (from output_cache.py "Cache HIT for <tool> (genome=...)")
            cache_match = CACHE_HIT_RE.search(line)
            if cache_match:
                rule_name, cache_genome = cache_match.groups()
                job_store.add_slurm_job(job_id, "—", rule_name, genome=cache_genome, source="from cache")
                # Immediately mark as CACHED so the checker skips it
                for sj in job_store.get_slurm_jobs(job_id):
                    if sj["rule"] == rule_name and sj["job_id"] == "—":
                        sj["status"] = "CACHED"

            # Snakemake's own --verbose "wildcards: genome=..." line, printed
            # right before it submits that same job -- see last_genome's
            # declaration above for why a single var is safe here.
            wildcards_match = WILDCARDS_GENOME_RE.search(line)
            if wildcards_match:
                last_genome = wildcards_match.group(1)

            # Parse SLURM job IDs as they appear in the log stream. log_path
            # is still captured for _slurm_status_checker's fallback backfill
            # (batch rules like quast_batch/gtdbtk_batch have no genome
            # wildcard, so last_genome is correctly empty for those).
            match = SLURM_SUBMIT_RE.search(line)
            if match:
                slurm_id, log_path = match.groups()
                rule_match = RULE_NAME_FROM_LOG_PATH_RE.search(log_path)
                ungrouped_rule_name, group_name = rule_match.groups() if rule_match else (None, None)
                job_store.add_slurm_job(job_id, slurm_id, ungrouped_rule_name or group_name or "unknown", genome=last_genome, log_path=log_path)
            elif SLURM_SUBMIT_FALLBACK_RE.search(line):
                fallback = SLURM_SUBMIT_FALLBACK_RE.search(line)
                slurm_id = fallback.group(1)
                job_store.add_slurm_job(job_id, slurm_id, "unknown", genome=last_genome)

            # Parse Snakemake step progress (e.g. "2 of 4 steps (50%) done")
            progress_match = STEPS_PROGRESS_RE.search(line)
            if progress_match:
                done, total, pct = progress_match.groups()
                job_store.update(job_id, steps_done=int(done), steps_total=int(total), progress=int(pct))

            # Parse the sequential orchestrator's own genome-transition marker
            genome_match = SEQUENTIAL_GENOME_RE.search(line)
            if genome_match:
                genome_index, genome_total = genome_match.groups()
                job_store.update(job_id, genome_index=int(genome_index), genome_total=int(genome_total))

            # Update phase based on output
            if "snakemake" in line.lower():
                job_store.update(job_id, phase="Running Snakemake")

        # Check exit code and mark as failed if non-zero
        if exit_code != 0:
            job_store.append_log(job_id, f"\n\n=== Command exited with code {exit_code} ===")
            job_store.finalize(job_id, status="failed", phase="Failed")
            LOGGER.error("Job %s failed with exit code %d", job_id, exit_code)
        else:
            job_store.finalize(job_id, status="completed", phase="Done")
    except Exception as e:
        job_store.append_log(job_id, f"\nError: {str(e)}")
        job_store.finalize(job_id, status="failed", phase="Error")


def submit_job(job_id: str, command: str, connection: SSHConnection):
    """Submit a job to the thread pool executor."""
    executor.submit(run_ssh_task, job_id, command, connection)


async def job_status_generator(job_id: str):
    """Generator that yields SSE events with job status updates."""
    last_state = None
    last_update_time = 0
    start_time = asyncio.get_event_loop().time()

    while True:
        try:
            status = ssh_slurm.check_slurm_job_status(job_id)
            current_state = status['state']
            elapsed = status['elapsed_time']

            current_time = asyncio.get_event_loop().time()
            check_duration = int(current_time - start_time)

            should_send = (
                current_state != last_state or
                (current_time - last_update_time) >= 7
            )

            if should_send:
                message = f"Job {current_state.lower()} (elapsed: {elapsed}, checking for: {check_duration}s)"
                data = {'state': current_state, 'elapsed': elapsed, 'message': message}
                yield f"data: {json.dumps(data)}\n\n"
                last_state = current_state
                last_update_time = current_time

            if current_state in ['COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'NOT_FOUND']:
                data = {'state': current_state, 'elapsed': elapsed, 'done': True}
                yield f"data: {json.dumps(data)}\n\n"
                break

            await asyncio.sleep(10)

        except Exception as e:
            LOGGER.exception("Error checking job status")
            data = {'error': f'Error checking status: {str(e)}'}
            yield f"data: {json.dumps(data)}\n\n"
            break
