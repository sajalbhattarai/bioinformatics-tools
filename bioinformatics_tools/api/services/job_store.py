"""
In-memory job state management.

Provides a JobStore class that wraps the jobs dict with structured
methods for creating, reading, and updating job state. All job state
mutations go through this module.

job_history_client.py mirrors status/phase/work_dir to the user's main_database
over SSH so history survives restarts; persistence is opt-in per job and never raises.
"""
import datetime
import logging
import time

from bioinformatics_tools.api.services import job_history_client

LOGGER = logging.getLogger(__name__)

# Fields mirrored to history on every update(); slurm_jobs/containers are
# written by checkpoint() instead.
_PERSISTED_FIELDS = ("status", "phase", "work_dir")

# Seconds between mid-run provenance writes to history. Each write is an SSH
# round-trip; this caps the loss if the API stops before finalize().
_CHECKPOINT_MIN_INTERVAL = 30.0


class JobStore:
    """In-memory job state management."""

    def __init__(self):
        self._jobs: dict[str, dict] = {}
        # Kept outside self._jobs because job dicts are returned as JSON and an
        # SSHConnection is not serialisable. Keyed by job_id.
        self._persistence: dict[str, tuple[str, object]] = {}
        # When each job's provenance was last written, for checkpoint()'s throttle.
        self._last_checkpoint: dict[str, float] = {}

    def create(self, job_id: str, genome_path: str, user_id: int | None = None,
               workflow: str | None = None, output_dir: str | None = None,
               selected_tools: str | None = None, relaunched_from: str | None = None,
               persist_owner_username: str | None = None,
               persist_owner_cluster_username: str | None = None,
               persist_db_path: str | None = None, persist_connection=None) -> dict:
        """Initialises a new job entry with all default fields.

        selected_tools (None = all tools) and relaunched_from (the job it was
        resumed or restarted from) are kept in memory too. persist_db_path /
        persist_connection mirror state changes to the user's job history; the
        self-test workflows omit them and get no history entry.
        """
        job = {
            "job_id": job_id,
            "user_id": user_id,
            "status": "pending",
            "phase": "Initializing",
            "genome_path": genome_path,
            "workflow": workflow,
            "selected_tools": selected_tools,
            "relaunched_from": relaunched_from,
            "sub_jobs": [],
            "slurm_jobs": [],
            "containers": [],
            "work_dir": None,
            "start_time": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        self._jobs[job_id] = job
        LOGGER.info("Created job %s", job_id)

        if persist_db_path and persist_connection:
            self._persistence[job_id] = (persist_db_path, persist_connection)
            try:
                job_history_client.record_job_created(
                    persist_connection, persist_db_path, job_id,
                    workflow or "unknown", genome_path, output_dir,
                    owner_username=persist_owner_username,
                    owner_cluster_username=persist_owner_cluster_username,
                    selected_tools=selected_tools, relaunched_from=relaunched_from,
                )
            except Exception as exc:
                LOGGER.warning("Could not record job %s in persistent history: %s", job_id, exc)

        return job

    def attach_persistence(self, job_id: str, db_path: str, connection) -> None:
        """Routes this job's later status/phase/work_dir changes to history without
        recording a creation; used when reattaching after a dane-api restart."""
        if job_id in self._jobs and db_path and connection is not None:
            self._persistence[job_id] = (db_path, connection)

    def get(self, job_id: str) -> dict | None:
        """Get a job by ID, or None if not found."""
        return self._jobs.get(job_id)

    def exists(self, job_id: str) -> bool:
        return job_id in self._jobs

    def update(self, job_id: str, **fields):
        """Updates fields on a job, mirroring status/phase/work_dir changes to history."""
        job = self._jobs.get(job_id)
        if job is None:
            return

        changed_persisted = {
            k: v for k, v in fields.items()
            if k in _PERSISTED_FIELDS and job.get(k) != v
        }
        job.update(fields)

        if changed_persisted and job_id in self._persistence:
            db_path, connection = self._persistence[job_id]
            try:
                job_history_client.record_job_updated(
                    connection, db_path, job_id, **changed_persisted,
                )
            except Exception as exc:
                LOGGER.warning("Could not persist update for job %s: %s", job_id, exc)

    def checkpoint(self, job_id: str, force: bool = False) -> bool:
        """Writes this job's SLURM jobs and containers to history mid-run.

        Returns True if a write happened. Throttled to one write per
        _CHECKPOINT_MIN_INTERVAL seconds unless force=True. Never raises, since
        it runs on the log-parsing hot path.
        """
        job = self._jobs.get(job_id)
        if job is None or job_id not in self._persistence:
            return False

        now = time.monotonic()
        if not force and now - self._last_checkpoint.get(job_id, 0.0) < _CHECKPOINT_MIN_INTERVAL:
            return False

        # Stamped before the write so a hanging SSH call does not let the next caller through.
        self._last_checkpoint[job_id] = now

        db_path, connection = self._persistence[job_id]
        try:
            job_history_client.record_job_updated(
                connection, db_path, job_id,
                slurm_jobs=job.get("slurm_jobs", []),
                containers=job.get("containers", []),
            )
            return True
        except Exception as exc:
            LOGGER.warning("Could not checkpoint provenance for job %s: %s", job_id, exc)
            return False

    def finalize(self, job_id: str, status: str, phase: str) -> None:
        """Marks a job completed or failed and persists a final snapshot.

        Always writes logs/slurm_jobs/containers, which update()'s change
        detection would miss (same object). The full log is written only here;
        a run that dies mid-way keeps its log in ~/.local/share/bsp/jobs/<job_id>.log.
        """
        job = self._jobs.get(job_id)
        if job is None:
            return
        job["status"] = status
        job["phase"] = phase
        job["end_time"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        LOGGER.info("Finalized job %s: status=%s phase=%s", job_id, status, phase)

        if job_id in self._persistence:
            db_path, connection = self._persistence[job_id]
            try:
                job_history_client.record_job_updated(
                    connection, db_path, job_id,
                    status=status, phase=phase,
                    logs=job.get("logs", ""),
                    slurm_jobs=job.get("slurm_jobs", []),
                    containers=job.get("containers", []),
                )
            except Exception as exc:
                LOGGER.warning("Could not persist final snapshot for job %s: %s", job_id, exc)

    def append_log(self, job_id: str, line: str):
        """Append a line to a job's log output."""
        if job_id in self._jobs:
            self._jobs[job_id]["logs"] = self._jobs[job_id].get("logs", "") + line + "\n"

    def add_slurm_job(self, job_id: str, slurm_id: str, rule: str, genome: str = "", source: str = "fresh run", log_path: str = ""):
        """Registers a newly discovered SLURM sub-job. log_path (internal) is read
        later by _slurm_status_checker to backfill the genome (see ssh_slurm.get_job_genome)."""
        if job_id not in self._jobs:
            return
        rows = self._jobs[job_id]["slurm_jobs"]
        # Ignores a job already registered: the workflow log emits every line
        # twice (two logger handlers). Cache hits share the "—" placeholder id,
        # so rule+genome identifies those rows.
        if slurm_id and slurm_id != "—":
            if any(r["job_id"] == slurm_id for r in rows):
                return
        elif any(r["job_id"] == slurm_id and r["rule"] == rule and r.get("genome") == genome
                 for r in rows):
            return
        rows.append({
            "job_id": slurm_id,
            "rule": rule,
            "status": "SUBMITTED",
            "time": "00:00:00",
            "genome": genome,
            "source": source,
            "log_path": log_path,
        })
        # Reached only for new rows, so at most one throttled call per submission.
        self.checkpoint(job_id)

    def add_container(self, job_id: str, container_info: dict):
        """Registers a container discovered from log parsing.

        Deduplicated by name + version + path, since every log line arrives
        twice and the list answers "what did this run use".
        """
        if job_id not in self._jobs:
            return
        rows = self._jobs[job_id]["containers"]
        key = (container_info.get("name"), container_info.get("version"),
               container_info.get("path"))
        if any((c.get("name"), c.get("version"), c.get("path")) == key for c in rows):
            return
        rows.append(container_info)
        self.checkpoint(job_id)

    def get_slurm_jobs(self, job_id: str) -> list[dict]:
        """Get the slurm_jobs list for a job."""
        job = self._jobs.get(job_id)
        return job.get("slurm_jobs", []) if job else []

    def get_status(self, job_id: str) -> str | None:
        """Get just the status field for a job."""
        job = self._jobs.get(job_id)
        return job.get("status") if job else None

    def cancel(self, job_id: str) -> None:
        """Mark a job as cancelled."""
        if job_id in self._jobs:
            self.update(job_id, status="cancelled", phase="Cancelled by user")
            LOGGER.info("Cancelled job %s", job_id)


# Module-level singleton
job_store = JobStore()
