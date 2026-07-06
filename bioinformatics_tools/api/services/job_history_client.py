"""
API-side client for job_history.py.
Saves and retrieves job history by SSHing into the cluster and running
job_history.py there. The job history database lives on the cluster, so
the API server cannot open it directly.
Errors here are ignored — if saving history fails, the job itself keeps
running.
"""
import json
import logging

from bioinformatics_tools.utilities.ssh_connection import SSHConnection

LOGGER = logging.getLogger(__name__)

_REMOTE_MODULE = "bioinformatics_tools.api.services.job_history"
_REMOTE_PYTHON = "~/bioinformatics-tools/.venv/bin/python"


def _run(connection: SSHConnection, action: str, payload: dict):
    ssh = connection.connect()
    try:
        stdin, stdout, stderr = ssh.exec_command(f"{_REMOTE_PYTHON} -m {_REMOTE_MODULE} {action}")
        stdin.write(json.dumps(payload))
        stdin.channel.shutdown_write()
        out = stdout.read().decode()
        exit_code = stdout.channel.recv_exit_status()
        if exit_code != 0:
            err = stderr.read().decode()
            LOGGER.warning("job_history %s failed (exit %d): %s", action, exit_code, err)
            return None
        out = out.strip()
        return json.loads(out) if out else None
    except Exception as exc:
        LOGGER.warning("job_history %s failed: %s", action, exc)
        return None
    finally:
        ssh.close()


def record_job_created(connection: SSHConnection, db_path: str, job_id: str, workflow: str,
                        genome_path: str | None, output_dir: str | None,
                        selected_tools: str | None = None,
                        relaunched_from: str | None = None) -> None:
    _run(connection, "create", {
        "db_path": db_path, "job_id": job_id, "workflow": workflow,
        "genome_path": genome_path, "output_dir": output_dir,
        "selected_tools": selected_tools, "relaunched_from": relaunched_from,
    })


def record_job_updated(connection: SSHConnection, db_path: str, job_id: str, **fields) -> None:
    _run(connection, "update", {"db_path": db_path, "job_id": job_id, "fields": fields})


def get_job(connection: SSHConnection, db_path: str, job_id: str) -> dict | None:
    return _run(connection, "get", {"db_path": db_path, "job_id": job_id})


def list_jobs(connection: SSHConnection, db_path: str, workflow: str | None = None,
              limit: int = 100, offset: int = 0) -> list[dict]:
    result = _run(connection, "list", {
        "db_path": db_path, "workflow": workflow, "limit": limit, "offset": offset,
    })
    return result if result is not None else []


def list_jobs_and_count(connection: SSHConnection, db_path: str, workflow: str | None = None,
                        limit: int = 100, offset: int = 0) -> tuple[list[dict], int]:
    """Single SSH call returning (jobs, total) — avoids the separate count round-trip."""
    result = _run(connection, "list_and_count", {
        "db_path": db_path, "workflow": workflow, "limit": limit, "offset": offset,
    })
    if result is None:
        return [], 0
    return result.get("jobs", []), result.get("total", 0)


def count_jobs(connection: SSHConnection, db_path: str, workflow: str | None = None) -> int:
    result = _run(connection, "count", {"db_path": db_path, "workflow": workflow})
    return result if result is not None else 0
