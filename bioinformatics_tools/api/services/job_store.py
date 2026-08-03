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
import time

from bioinformatics_tools.api.services import job_history_client

LOGGER = logging.getLogger(__name__)

# Fields worth round-tripping to the persistent history table on every
# update(); everything else (logs, sub_jobs, report, steps_done/total,
# progress) is live-session-only detail. slurm_jobs/containers are neither:
# see checkpoint() below.
_PERSISTED_FIELDS = ("status", "phase", "work_dir")

# How often a run's SLURM/container provenance is written to history while it
# is still going.
#
# It used to be written once, by finalize(). That made provenance depend on the
# API surviving to the end of the run, and it does not always: when the SSH log
# stream drops ("'NoneType' object has no attribute 'open_session'"), the job
# was finalized four minutes into a run whose in-memory list was still empty,
# so an empty list is what got snapshotted -- and the SLURM Jobs table read
# empty forever afterwards, for a run that went on to submit 458 jobs and whose
# log recorded every one of them. Three consecutive runs lost their provenance
# that way.
#
# Every history write is an SSH round-trip, so per-log-line persistence really
# would be too expensive -- but per-30-seconds is not. That is a fraction of
# what the SLURM status checker beside it already spends on its own 15s cycle,
# and it caps the loss at the last half-minute instead of everything.
_CHECKPOINT_MIN_INTERVAL = 30.0


class JobStore:
    """In-memory job state management."""

    def __init__(self):
        self._jobs: dict[str, dict] = {}
        # Kept OUT of self._jobs deliberately: job_status endpoints spread a
        # job dict's fields straight into a JSON response, and an
        # SSHConnection object isn't JSON-serializable. Keyed by job_id.
        self._persistence: dict[str, tuple[str, object]] = {}
        # When each job's provenance was last written to history, so
        # checkpoint() can throttle. Same reason as above for living outside
        # self._jobs: it is bookkeeping, not part of the API response.
        self._last_checkpoint: dict[str, float] = {}

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

    def attach_persistence(self, job_id: str, db_path: str, connection) -> None:
        """Route this job's later status/phase/work_dir changes to history
        WITHOUT recording a creation, for a job whose history row already
        exists. Used when reattaching to a run after a dane-api restart:
        create() would insert a duplicate row for a job that has been in the
        table since it was first launched."""
        if job_id in self._jobs and db_path and connection is not None:
            self._persistence[job_id] = (db_path, connection)

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

    def checkpoint(self, job_id: str, force: bool = False) -> bool:
        """Write this job's SLURM jobs and containers to history mid-run.

        Returns True if a write happened. Throttled to one write per
        _CHECKPOINT_MIN_INTERVAL seconds unless force=True, because each one
        costs an SSH round-trip.

        This exists so provenance no longer depends on the run reaching
        finalize(). A dropped log stream, a killed API, a crash -- all used to
        take the whole SLURM Jobs table with them, and the job page then showed
        an empty table permanently, because the history row it falls back to
        had never been given anything to show.

        Never raises: a job whose provenance cannot be written is still a job
        that is running fine, and this is called from the log-parsing hot path.
        """
        job = self._jobs.get(job_id)
        if job is None or job_id not in self._persistence:
            return False

        now = time.monotonic()
        if not force and now - self._last_checkpoint.get(job_id, 0.0) < _CHECKPOINT_MIN_INTERVAL:
            return False

        # Stamped BEFORE the write, not after: a slow or hanging SSH round-trip
        # would otherwise leave the previous stamp in place and let the next
        # caller straight through, stacking up round-trips on a link that is
        # already struggling.
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
        """Mark a job done (completed/failed) and persist a final snapshot
        of its logs/slurm_jobs/containers alongside the status/phase change.

        Unlike update(), this always writes the snapshot fields regardless
        of whether they "changed" -- they're a point-in-time capture of
        whatever streaming-session detail has accumulated so far, not an
        incremental delta, so update()'s job.get(k) != v change-detection
        would never fire for them (the in-memory value and the value being
        "set" are the same object).

        This is no longer the ONLY point where slurm_jobs/containers reach
        history -- checkpoint() writes them as they are discovered, so a run
        that never gets here still leaves its provenance behind. What is still
        unique to finalize() is `logs`, which stays end-of-run only: it is the
        whole log, rewritten in full on every write, and sending megabytes over
        SSH every half-minute is exactly the cost checkpoint() is shaped to
        avoid. A job that dies mid-run therefore keeps its SLURM jobs and
        containers, and loses only the log -- which is still on the cluster
        under ~/.local/share/bsp/jobs/<job_id>.log, and is what the reattach
        path replays.
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
        if job_id not in self._jobs:
            return
        rows = self._jobs[job_id]["slurm_jobs"]
        # Ignore a job already registered. A SLURM job id is unique per
        # submission, so seeing one twice never means two jobs -- and the
        # workflow's log emits every line twice (two handlers on the
        # workflow_tools logger, one timestamp-prefixed and one not), so every
        # submission and every cache hit was parsed twice: the SLURM Jobs table
        # listed a run's 759 jobs as 1518, each organism appearing under each
        # rule twice over. Deduping here rather than at the parser covers every
        # route into the table, and stays correct once the double logging is
        # itself fixed.
        #
        # Cache hits share the "—" placeholder id (nothing was submitted), so
        # for those it is the rule+genome pair that identifies the row.
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
        # Only reached when the row is genuinely new (every early return above
        # is a duplicate), so this is at most one throttled call per real
        # submission -- and the throttle makes it far fewer than that.
        self.checkpoint(job_id)

    def add_container(self, job_id: str, container_info: dict):
        """Register a container discovered from log parsing.

        Deduped for the same reason add_slurm_job is: the workflow's log emits
        every line twice (two handlers on the workflow_tools logger, one
        timestamp-prefixed and one not), so each bapptainer __CONTAINER__ line
        was parsed twice and the Containers box listed every image twice over.
        add_slurm_job had this guard and this did not, which is why the SLURM
        table looked right while the container list did not.

        A container is identified by name + version + path: the same image used
        by ten rules is still one image, and this list answers "what did this
        run use", not "how many times was it used". Deduping here rather than
        at the parser covers every route into the list, and stays correct once
        the double logging is itself fixed.
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
