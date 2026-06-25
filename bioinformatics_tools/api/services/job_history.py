"""
Persistent job history, backed by the user's own main_database SQLite file.

job_store.py is deliberately in-memory only (fast, simple, good enough for
actively-streaming logs) -- but that means every dane-api restart wipes all
job records, and there's no way to list "my past jobs" at all. This module
is the persistence layer underneath it: every job gets a row here at
creation, kept in sync as job_store.update() reports status/work_dir
changes, so job history survives restarts and can be listed/resumed.

Deliberately a SEPARATE table from workflow_tools/output_cache.py's
run_log -- run_log is written by the workflow process itself (even for
direct CLI runs with no API involved at all) and keyed on its own run_id;
this table is keyed on the API's job_id specifically, the identifier the
front-end actually navigates with (/jobs/{job_id}).

main_database is a path on the user's cluster account, not necessarily
the machine dane-api itself runs on (each user has their own
cluster_host/cluster_username) -- so the API never opens this file
directly. It invokes this module's functions over SSH instead, as
`python -m bioinformatics_tools.api.services.job_history <action>` with a
JSON payload on stdin and a JSON result (if any) on stdout; see
api/services/job_history_client.py for the calling side.
"""
import json
import logging
import os
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

LOGGER = logging.getLogger(__name__)

CREATE_API_JOBS_SQL = """
CREATE TABLE IF NOT EXISTS api_jobs (
    job_id TEXT PRIMARY KEY,
    workflow TEXT NOT NULL,
    genome_path TEXT,
    output_dir TEXT,
    work_dir TEXT,
    status TEXT NOT NULL,
    phase TEXT,
    selected_tools TEXT,
    relaunched_from TEXT,
    logs TEXT,
    slurm_jobs TEXT,
    containers TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
"""

# Columns added after the table's initial release -- CREATE TABLE IF NOT
# EXISTS above only helps brand-new databases; existing deployed ones need
# an explicit ALTER (SQLite has no ADD COLUMN IF NOT EXISTS).
_ADDED_COLUMNS = ("selected_tools", "relaunched_from", "logs", "slurm_jobs", "containers")

# update() may be called with any subset of job_store's fields. status/
# phase/work_dir are persisted on every change (job_store.update()'s own
# change-detection decides when). logs/slurm_jobs/containers are NEVER
# persisted incrementally (that would mean an SSH round-trip per log
# line) -- they only ever arrive here via job_store.finalize()'s one-time
# snapshot at job completion/failure. sub_jobs/report/steps_done/total/
# progress remain pure live-session detail, never persisted at all.
_PERSISTED_UPDATE_FIELDS = ("status", "phase", "work_dir", "logs", "slurm_jobs", "containers")

# Fields whose value is a Python list/dict (slurm_jobs, containers) need
# JSON-encoding before they can be bound as a SQLite TEXT column -- plain
# strings (status/phase/work_dir/logs) pass through unchanged.
_JSON_ENCODED_FIELDS = ("slurm_jobs", "containers")


def _get_connection(db_path: str, timeout: float = 30.0) -> sqlite3.Connection:
    """SQLite connection configured for network filesystems (same pattern as
    workflow_tools/output_cache.py's _get_connection). This module always
    runs on the remote cluster account (invoked over SSH -- see the module
    docstring), so os.path.expanduser resolves a leading ~ against the
    right user every time, regardless of who dane-api itself runs as.
    """
    conn = sqlite3.connect(os.path.expanduser(db_path), timeout=timeout)
    conn.execute(f"PRAGMA busy_timeout={int(timeout * 1000)}")
    return conn


def ensure_table(db_path: str) -> None:
    if not db_path:
        return
    try:
        conn = _get_connection(db_path)
        try:
            conn.execute(CREATE_API_JOBS_SQL)
            existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(api_jobs)")}
            for col in _ADDED_COLUMNS:
                if col not in existing_cols:
                    conn.execute(f"ALTER TABLE api_jobs ADD COLUMN {col} TEXT")
            conn.commit()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not ensure api_jobs table at %s: %s", db_path, exc)


def record_job_created(db_path: str, job_id: str, workflow: str,
                       genome_path: str | None, output_dir: str | None,
                       selected_tools: str | None = None,
                       relaunched_from: str | None = None) -> None:
    if not db_path:
        return
    ensure_table(db_path)
    now = datetime.now(timezone.utc).isoformat()
    try:
        conn = _get_connection(db_path)
        try:
            conn.execute(
                "INSERT OR REPLACE INTO api_jobs "
                "(job_id, workflow, genome_path, output_dir, work_dir, status, phase, "
                " selected_tools, relaunched_from, created_at, updated_at) "
                "VALUES (?, ?, ?, ?, NULL, 'pending', 'Initializing', ?, ?, ?, ?)",
                (job_id, workflow, genome_path, output_dir, selected_tools, relaunched_from, now, now),
            )
            conn.commit()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not record job %s in history: %s", job_id, exc)


def record_job_updated(db_path: str, job_id: str, **fields) -> None:
    """Persist only the subset of fields that matter for history (see
    _PERSISTED_UPDATE_FIELDS); a no-op if nothing relevant changed."""
    if not db_path:
        return
    relevant = {k: v for k, v in fields.items() if k in _PERSISTED_UPDATE_FIELDS}
    if not relevant:
        return
    relevant = {
        k: (json.dumps(v) if k in _JSON_ENCODED_FIELDS else v)
        for k, v in relevant.items()
    }
    set_clause = ", ".join(f"{k} = ?" for k in relevant)
    values = list(relevant.values()) + [datetime.now(timezone.utc).isoformat(), job_id]
    try:
        conn = _get_connection(db_path)
        try:
            conn.execute(
                f"UPDATE api_jobs SET {set_clause}, updated_at = ? WHERE job_id = ?",
                values,
            )
            conn.commit()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not update job %s in history: %s", job_id, exc)


def _decode_row(row: dict) -> dict:
    """Decode JSON-encoded list columns (see _JSON_ENCODED_FIELDS) back into
    real Python lists before this dict crosses the SSH/JSON transport
    boundary -- decoding here (not on the caller's side) means the outer
    json.dumps() in _main() serializes them as proper nested arrays
    instead of double-encoded strings."""
    for field in _JSON_ENCODED_FIELDS:
        raw = row.get(field)
        if not raw:
            row[field] = []
            continue
        try:
            row[field] = json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            row[field] = []
    return row


def get_job(db_path: str, job_id: str) -> dict | None:
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return None
    try:
        conn = _get_connection(db_path)
        try:
            conn.row_factory = sqlite3.Row
            row = conn.execute("SELECT * FROM api_jobs WHERE job_id = ?", (job_id,)).fetchone()
            return _decode_row(dict(row)) if row else None
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not read job %s from history: %s", job_id, exc)
        return None


def list_jobs(db_path: str, workflow: str | None = None, limit: int = 100, offset: int = 0) -> list[dict]:
    """Most recent jobs first, optionally filtered to one workflow, paginated
    via limit/offset. Empty list (not an error) if the table doesn't exist
    yet -- a brand new user with no history at all is the normal case, not
    a failure. Pair with count_jobs() for the total matching the same
    workflow filter, to compute page count."""
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return []
    try:
        conn = _get_connection(db_path)
        try:
            conn.row_factory = sqlite3.Row
            if workflow:
                rows = conn.execute(
                    "SELECT * FROM api_jobs WHERE workflow = ? ORDER BY created_at DESC LIMIT ? OFFSET ?",
                    (workflow, limit, offset),
                ).fetchall()
            else:
                rows = conn.execute(
                    "SELECT * FROM api_jobs ORDER BY created_at DESC LIMIT ? OFFSET ?",
                    (limit, offset),
                ).fetchall()
            return [_decode_row(dict(r)) for r in rows]
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not list jobs from history at %s: %s", db_path, exc)
        return []


def count_jobs(db_path: str, workflow: str | None = None) -> int:
    """Total number of history rows matching workflow (or all rows if
    None) -- used alongside list_jobs's limit/offset to compute total
    page count on the API side."""
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return 0
    try:
        conn = _get_connection(db_path)
        try:
            if workflow:
                row = conn.execute("SELECT COUNT(*) FROM api_jobs WHERE workflow = ?", (workflow,)).fetchone()
            else:
                row = conn.execute("SELECT COUNT(*) FROM api_jobs").fetchone()
            return row[0] if row else 0
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not count jobs from history at %s: %s", db_path, exc)
        return 0


# ---------------------------------------------------------------------------
# CLI entry point: `python -m bioinformatics_tools.api.services.job_history
# <action>`, payload as JSON on stdin, result (if any) as JSON on stdout.
# This is how the API process (which may run on a different machine than
# the user's cluster account) reaches this module -- by running it
# remotely over SSH rather than importing it in-process. See the module
# docstring and api/services/job_history_client.py.
# ---------------------------------------------------------------------------

def _main() -> int:
    if len(sys.argv) != 2:
        print("usage: python -m bioinformatics_tools.api.services.job_history "
              "<create|update|get|list|count>", file=sys.stderr)
        return 2
    action = sys.argv[1]
    payload = json.loads(sys.stdin.read() or "{}")
    db_path = payload.get("db_path")

    if action == "create":
        record_job_created(
            db_path, payload["job_id"], payload["workflow"],
            payload.get("genome_path"), payload.get("output_dir"),
            payload.get("selected_tools"), payload.get("relaunched_from"),
        )
    elif action == "update":
        record_job_updated(db_path, payload["job_id"], **payload.get("fields", {}))
    elif action == "get":
        result = get_job(db_path, payload["job_id"])
        print(json.dumps(result))
    elif action == "list":
        result = list_jobs(
            db_path, payload.get("workflow"), payload.get("limit", 100), payload.get("offset", 0),
        )
        print(json.dumps(result))
    elif action == "count":
        result = count_jobs(db_path, payload.get("workflow"))
        print(json.dumps(result))
    else:
        print(f"unknown action: {action}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(_main())
