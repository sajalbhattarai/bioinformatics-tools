"""
In-memory job state management.

Provides a JobStore class that wraps the jobs dict with structured
methods for creating, reading, and updating job state. All job state
mutations go through this module.

This in-memory store is wiped on every dane-api restart and has no way to
list a user's past jobs. job_history_client.py persists a small subset of
fields (status/phase/work_dir) to the user's own main_database over SSH so
that history survives restarts -- create()/update() are the single
chokepoint for all job state changes, so hooking persistence in here means
job_runner.py (which calls update() many times per job) needs no changes
at all. Persistence is opportunistic: it only fires when a job was created
with persist_db_path/persist_connection, and a failure to persist never
raises -- see job_history_client.py's module docstring.
"""
import datetime
import logging

from bioinformatics_tools.api.services import job_history_client

LOGGER = logging.getLogger(__name__)

# Fields worth round-tripping to the persistent history table; everything
# else (logs, sub_jobs, slurm_jobs, containers, report, steps_done/total,
# progress) is live-session-only detail.
_PERSISTED_FIELDS = ("status", "phase", "work_dir")


class JobStore:
    """In-memory job state management."""

    def __init__(self):
        self._jobs: dict[str, dict] = {}
        # Kept OUT of self._jobs deliberately: job_status endpoints spread a
        # job dict's fields straight into a JSON response, and an
        # SSHConnection object isn't JSON-serializable. Keyed by job_id.
        self._persistence: dict[str, tuple[str, object]] = {}

    def create(self, job_id: str, genome_path: str, user_id: int | None = None,
               workflow: str | None = None, output_dir: str | None = None,
               selected_tools: str | None = None, relaunched_from: str | None = None,
               persist_db_path: str | None = None, persist_connection=None) -> dict:
        """Initialize a new job entry with all default fields.

        selected_tools (comma-joined tool keys, None = "ran everything") and
        relaunched_from (the job_id this was resumed or restarted from, None
        for a fresh job) are stored in-memory too, not just persisted, so a
        still-live job's get_job_status response can surface them the same
        way it already does for workflow.

        persist_db_path/persist_connection are optional: when given, this
        job's status/phase/work_dir changes are mirrored to the user's
        persistent job history (see job_history_client.py). Self-test
        workflows (quick_example, fresh_test) don't pass these and simply
        get no history entry, which is correct -- they're not part of the
        user-facing workflow list to begin with.
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
                    selected_tools=selected_tools, relaunched_from=relaunched_from,
                )
            except Exception as exc:
                LOGGER.warning("Could not record job %s in persistent history: %s", job_id, exc)

        return job

    def get(self, job_id: str) -> dict | None:
        """Get a job by ID, or None if not found."""
        return self._jobs.get(job_id)

    def exists(self, job_id: str) -> bool:
        return job_id in self._jobs

    def update(self, job_id: str, **fields):
        """Update one or more fields on a job, mirroring any changed
        status/phase/work_dir to persistent history (if this job has
        persistence configured)."""
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

    def finalize(self, job_id: str, status: str, phase: str) -> None:
        """Mark a job done (completed/failed) and persist a final snapshot
        of its logs/slurm_jobs/containers alongside the status/phase change.

        Unlike update(), this always writes the snapshot fields regardless
        of whether they "changed" -- they're a point-in-time capture of
        whatever streaming-session detail has accumulated so far, not an
        incremental delta, so update()'s job.get(k) != v change-detection
        would never fire for them (the in-memory value and the value being
        "set" are the same object). This is the one point where
        streaming-only session detail (logs/slurm_jobs/containers, never
        persisted incrementally -- that would mean an SSH round-trip per
        log line) gets a chance to survive a dane-api restart and be visible
        again when a job is later resumed from history.

        A job that dies without ever reaching this method (e.g. dane-api
        itself crashes mid-run) gets no final snapshot, same limitation
        status/phase already have today -- not a new gap this introduces.
        """
        job = self._jobs.get(job_id)
        if job is None:
            return
        job["status"] = status
        job["phase"] = phase
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
        """Register a newly discovered SLURM sub-job. log_path (internal,
        not sent to the frontend table) is read later by
        _slurm_status_checker to backfill genome once that job's own log
        file exists -- see ssh_slurm.get_job_genome's docstring."""
        if job_id in self._jobs:
            self._jobs[job_id]["slurm_jobs"].append({
                "job_id": slurm_id,
                "rule": rule,
                "status": "SUBMITTED",
                "time": "00:00:00",
                "genome": genome,
                "source": source,
                "log_path": log_path,
            })

    def add_container(self, job_id: str, container_info: dict):
        """Register a container discovered from log parsing."""
        if job_id in self._jobs:
            self._jobs[job_id]["containers"].append(container_info)

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
