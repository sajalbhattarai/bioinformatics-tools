"""
Persistent job history in the user's main_database SQLite file (table api_jobs),
keyed on the API's job_id and kept in sync with the in-memory job_store.
Runs on the cluster account, called over SSH as `python -m ... job_history <action>`
with JSON on stdin/stdout (see job_history_client.py).
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
    owner_username TEXT,
    owner_cluster_username TEXT,
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

# Columns added after the first release; older databases get them by ALTER TABLE.
_ADDED_COLUMNS = (
    "owner_username",
    "owner_cluster_username",
    "selected_tools",
    "relaunched_from",
    "logs",
    "slurm_jobs",
    "containers",
)

# job_store fields kept in history; logs/slurm_jobs/containers arrive only in the final snapshot.
_PERSISTED_UPDATE_FIELDS = ("status", "phase", "work_dir", "logs", "slurm_jobs", "containers")

# List/dict fields stored as JSON text.
_JSON_ENCODED_FIELDS = ("slurm_jobs", "containers")


def _get_connection(db_path: str, timeout: float = 30.0) -> sqlite3.Connection:
    """Opens a SQLite connection with a busy timeout suited to network filesystems.
    Runs as the cluster account, so ~ expands to that user's home."""
    conn = sqlite3.connect(os.path.expanduser(db_path), timeout=timeout)
    conn.execute(f"PRAGMA busy_timeout={int(timeout * 1000)}")
    return conn


def ensure_table(db_path: str) -> None:
    """Creates api_jobs if needed and adds any missing later columns."""
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
                       owner_username: str | None = None,
                       owner_cluster_username: str | None = None,
                       selected_tools: str | None = None,
                       relaunched_from: str | None = None) -> None:
    """Inserts a new pending job row."""
    if not db_path:
        return
    ensure_table(db_path)
    now = datetime.now(timezone.utc).isoformat()
    try:
        conn = _get_connection(db_path)
        try:
            conn.execute(
                "INSERT OR REPLACE INTO api_jobs "
                "(job_id, owner_username, owner_cluster_username, workflow, genome_path, output_dir, "
                " work_dir, status, phase, selected_tools, relaunched_from, created_at, updated_at) "
                "VALUES (?, ?, ?, ?, ?, ?, NULL, 'pending', 'Initializing', ?, ?, ?, ?)",
                (
                    job_id,
                    owner_username,
                    owner_cluster_username,
                    workflow,
                    genome_path,
                    output_dir,
                    selected_tools,
                    relaunched_from,
                    now,
                    now,
                ),
            )
            conn.commit()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not record job %s in history: %s", job_id, exc)


def record_job_updated(db_path: str, job_id: str, **fields) -> None:
    """Writes the history fields among `fields`; does nothing if none are present."""
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
    """Decodes the JSON text columns back into lists so the output is not double-encoded."""
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


def _ownership_where(owner_username: str | None,
                     owner_cluster_username: str | None) -> tuple[str, list[str]]:
    """Returns an SQL ownership condition and its params.
    Rows without owner_username match on owner_cluster_username or a '/<cluster_user>/' path segment."""
    if not owner_username:
        return "", []

    if owner_cluster_username:
        like = f"%/{owner_cluster_username}/%"
        return (
            "(owner_username = ? "
            "OR (owner_username IS NULL AND owner_cluster_username = ?) "
            "OR (owner_username IS NULL AND owner_cluster_username IS NULL "
            "    AND (work_dir LIKE ? OR genome_path LIKE ?)))",
            [owner_username, owner_cluster_username, like, like],
        )

    return "(owner_username = ?)", [owner_username]


def get_job(db_path: str, job_id: str,
            owner_username: str | None = None,
            owner_cluster_username: str | None = None) -> dict | None:
    """Returns one job row the owner may see, or None."""
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return None
    try:
        conn = _get_connection(db_path)
        try:
            conn.row_factory = sqlite3.Row
            where_ownership, owner_params = _ownership_where(owner_username, owner_cluster_username)
            sql = "SELECT * FROM api_jobs WHERE job_id = ?"
            params: list[str] = [job_id]
            if where_ownership:
                sql += f" AND {where_ownership}"
                params.extend(owner_params)
            row = conn.execute(sql, params).fetchone()
            return _decode_row(dict(row)) if row else None
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not read job %s from history: %s", job_id, exc)
        return None


# Characters of each log's end included in list results.
LIST_LOG_TAIL = 8000


def list_jobs(db_path: str, workflow: str | None = None, limit: int = 100,
              offset: int = 0, owner_username: str | None = None,
              owner_cluster_username: str | None = None) -> list[dict]:
    """Returns one page of jobs, newest first, optionally for one workflow.
    Returns an empty list when the database or table does not exist yet."""
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return []
    try:
        conn = _get_connection(db_path)
        try:
            conn.row_factory = sqlite3.Row
            where_ownership, owner_params = _ownership_where(owner_username, owner_cluster_username)
            conditions: list[str] = []
            params: list[str | int] = []
            if workflow:
                conditions.append("workflow = ?")
                params.append(workflow)
            if where_ownership:
                conditions.append(where_ownership)
                params.extend(owner_params)

            where_clause = f"WHERE {' AND '.join(conditions)} " if conditions else ""
            # Lists carry only each log's tail and no SLURM jobs or containers; get_job has the full row.
            columns = [r[1] for r in conn.execute("PRAGMA table_info(api_jobs)")]
            select = ", ".join(
                f"substr(logs, -{LIST_LOG_TAIL}) AS logs" if c == "logs"
                else f'NULL AS "{c}"' if c in ("slurm_jobs", "containers")
                else f'"{c}"' for c in columns)
            rows = conn.execute(
                f"SELECT {select} FROM api_jobs {where_clause}ORDER BY created_at DESC LIMIT ? OFFSET ?",
                [*params, limit, offset],
            ).fetchall()
            return [_decode_row(dict(r)) for r in rows]
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not list jobs from history at %s: %s", db_path, exc)
        return []


def list_jobs_and_count(db_path: str, workflow: str | None = None,
                        limit: int = 100, offset: int = 0,
                        owner_username: str | None = None,
                        owner_cluster_username: str | None = None) -> dict:
    """Returns {"jobs", "total"} from list_jobs and count_jobs in one remote call."""
    return {
        "jobs": list_jobs(
            db_path,
            workflow=workflow,
            limit=limit,
            offset=offset,
            owner_username=owner_username,
            owner_cluster_username=owner_cluster_username,
        ),
        "total": count_jobs(
            db_path,
            workflow=workflow,
            owner_username=owner_username,
            owner_cluster_username=owner_cluster_username,
        ),
    }


def count_jobs(db_path: str, workflow: str | None = None,
               owner_username: str | None = None,
               owner_cluster_username: str | None = None) -> int:
    """Returns the number of jobs matching the workflow filter and owner."""
    if not db_path or not Path(os.path.expanduser(db_path)).exists():
        return 0
    try:
        conn = _get_connection(db_path)
        try:
            where_ownership, owner_params = _ownership_where(owner_username, owner_cluster_username)
            conditions: list[str] = []
            params: list[str] = []
            if workflow:
                conditions.append("workflow = ?")
                params.append(workflow)
            if where_ownership:
                conditions.append(where_ownership)
                params.extend(owner_params)
            where_clause = f" WHERE {' AND '.join(conditions)}" if conditions else ""
            row = conn.execute(f"SELECT COUNT(*) FROM api_jobs{where_clause}", params).fetchone()
            return row[0] if row else 0
        finally:
            conn.close()
    except sqlite3.Error as exc:
        LOGGER.warning("Could not count jobs from history at %s: %s", db_path, exc)
        return 0


# ---- CLI entry point: `<action>` argument, JSON payload on stdin, JSON result on stdout ----

def dispatch(action: str, payload: dict):
    """Runs one action for the CLI or the in-process client; returns None for create/update."""
    db_path = payload.get("db_path")
    if action == "create":
        record_job_created(
            db_path, payload["job_id"], payload["workflow"],
            payload.get("genome_path"), payload.get("output_dir"),
            payload.get("owner_username"), payload.get("owner_cluster_username"),
            payload.get("selected_tools"), payload.get("relaunched_from"),
        )
        return None
    if action == "update":
        record_job_updated(db_path, payload["job_id"], **payload.get("fields", {}))
        return None
    if action == "get":
        return get_job(db_path, payload["job_id"],
                       payload.get("owner_username"), payload.get("owner_cluster_username"))
    if action == "list":
        return list_jobs(db_path, payload.get("workflow"), payload.get("limit", 100), payload.get("offset", 0),
                         payload.get("owner_username"), payload.get("owner_cluster_username"))
    if action == "list_and_count":
        return list_jobs_and_count(db_path, payload.get("workflow"), payload.get("limit", 100),
                                   payload.get("offset", 0), payload.get("owner_username"),
                                   payload.get("owner_cluster_username"))
    if action == "count":
        return count_jobs(db_path, payload.get("workflow"),
                          payload.get("owner_username"), payload.get("owner_cluster_username"))
    raise ValueError(f"unknown action: {action}")


def _main() -> int:
    """Reads the action and payload, runs it and prints the JSON result."""
    if len(sys.argv) != 2:
        print("usage: python -m bioinformatics_tools.api.services.job_history "
              "<create|update|get|list|list_and_count|count>", file=sys.stderr)
        return 2
    action = sys.argv[1]
    payload = json.loads(sys.stdin.read() or "{}")
    try:
        result = dispatch(action, payload)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if action not in ("create", "update"):
        print(json.dumps(result))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
