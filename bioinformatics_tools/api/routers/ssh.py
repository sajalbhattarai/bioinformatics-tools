"""
SSH and job management endpoints.

Thin routing layer — delegates to job_store, job_runner, and ssh utilities.

All endpoints (except /health) require a valid Bearer token. The token is
validated by get_current_user(), which returns the user's cluster credentials.
_build_connection() decrypts the stored private key and builds a per-user
SSHConnection for each request.
"""
import io
import logging
import re
import uuid
from datetime import datetime

import pandas as pd
from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import StreamingResponse

from dataclasses import asdict

from bioinformatics_tools.api.auth import decrypt_private_key, get_current_user
from bioinformatics_tools.api.models import GenomeSend, SlurmSend
from bioinformatics_tools.api.services import job_history_client, job_runner
from bioinformatics_tools.api.services.job_store import job_store
from bioinformatics_tools.utilities import ssh_sftp, ssh_slurm
from bioinformatics_tools.utilities.ssh_connection import make_user_connection
from bioinformatics_tools.workflow_tools.workflow_helpers import GENOME_EXTENSIONS
from bioinformatics_tools.workflow_tools.workflow_registry import (
    MARGIE_SB_PHASED_TOOLS,
    WORKFLOWS,
    REQUIRED_SYSTEM_PARAMS,
    workflow_path_params,
)

LOGGER = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/ssh", tags=["ssh"])

# Workflows visible on the frontend but not yet implemented.
STUB_WORKFLOWS: set[str] = {"custom_microbiome"}

# (job_id, path) -> (mtime, size, total_lines). Wiped on API restart, same
# as job_store's in-memory state. Avoids re-running wc -l (a full file
# scan) on every page click for a file whose content hasn't changed --
# output files don't change once a job completes.
_line_count_cache: dict[tuple[str, str], tuple[float, int, int]] = {}


def _validate_relative_path(path: str, *, label: str = "file") -> None:
    """Raise HTTPException(400) if path attempts directory traversal.

    Shared by job_files (subdir), download_file (path), and view_file
    (path) -- all three take a user-supplied relative path under a job's
    work_dir.
    """
    if path and (path.startswith("/") or ".." in path.split("/")):
        raise HTTPException(status_code=400, detail=f"Invalid {label} path")


def _resolve_job_work_dir(job_id: str, current_user: dict, conn) -> str:
    """Resolve a job's work_dir, falling back to persistent history if the
    job isn't in the in-memory job_store (e.g. after a dane-api restart) --
    the same fallback get_job_status already uses, applied here so file
    browsing/download/view work for resumed/historical jobs too, not just
    job_status itself.

    Raises HTTPException(404) if the job can't be found anywhere (live or
    history), or HTTPException(400) if found but has no work_dir yet.
    """
    job = job_store.get(job_id)
    if job is not None:
        if job.get("user_id") != current_user["user_id"]:
            raise HTTPException(status_code=403, detail="Access denied")
        work_dir = job.get("work_dir")
    else:
        try:
            user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
        except Exception:
            raise HTTPException(status_code=404, detail="Job not found")

        main_db = user_config.get('main_database')
        row = job_history_client.get_job(conn, main_db, job_id) if main_db else None
        if row is None:
            raise HTTPException(status_code=404, detail="Job not found")
        work_dir = row.get("work_dir")

    if not work_dir:
        raise HTTPException(status_code=400, detail="No working directory available for this job")
    return work_dir


def _detect_delimiter(path: str, header: str) -> str:
    """Pick a column delimiter: by extension first, else sniff the header.

    Naive -- no RFC4180 quote-handling. Every sampled real output file is
    quote-free TSV; a CSV with delimiters embedded in quoted fields would
    misparse. Acceptable v1 limitation, not silently papered over.
    """
    lower = path.lower()
    if lower.endswith(".csv"):
        return ","
    if lower.endswith(".tsv"):
        return "\t"
    return "\t" if "\t" in header else ","


def _get_available_workflows() -> list[dict]:
    """
    Build the list of available workflows from WORKFLOWS registry.
    Returns detailed metadata for each workflow including tools, params, etc.
    Automatically merges REQUIRED_SYSTEM_PARAMS with workflow-specific params.
    """
    workflows = []

    # Add workflows from WORKFLOWS registry
    for wf_id, wf_key in WORKFLOWS.items():
        # Skip internal test workflows
        if wf_id in ['example', 'selftest']:
            continue

        # Convert dataclass to dict and add computed fields
        wf_dict = asdict(wf_key)
        wf_dict['id'] = wf_key.cmd_identifier
        wf_dict['containers'] = [{'name': sif[0], 'version': sif[1]} for sif in wf_key.sif_files]

        # Merge system-wide required params, this workflow's own root-path settings,
        # and the workflow's other params. System params come first since they're
        # infra-level. input_path/output_path apply to every workflow; sif_path and
        # db_root are only included when this workflow actually has a local-folder
        # sif lookup / a unified db_root fallback to point at (see
        # workflow_path_params()'s docstring for why that distinction matters).
        path_params = workflow_path_params(
            wf_id,
            include_sif=wf_key.local_sif_only,
            include_db_root=wf_key.supports_db_root,
            supports_batch_input=wf_key.supports_batch_input,
        )
        wf_dict['configurable_params'] = REQUIRED_SYSTEM_PARAMS + path_params + (wf_key.configurable_params or [])

        workflows.append(wf_dict)

    # Add stub workflows (not yet implemented but visible)
    # Even stub workflows get system params since they'll need them when implemented
    workflows.append({
        'id': 'custom_microbiome',
        'label': 'Custom Microbiome',
        'description': 'Custom microbiome annotation workflow (coming soon)',
        'full_description': 'A specialized workflow for microbiome annotation. This workflow is currently under development.',
        'tools': [],
        'configurable_params': REQUIRED_SYSTEM_PARAMS,  # Stub still needs system params
        'database_deps': [],
        'docs_url': None,
        'containers': [],
        'cmd_identifier': 'custom_microbiome',
        'snakemake_file': '',
        'other': [],
        'sif_files': [],
    })

    return workflows


def _build_connection(current_user: dict):
    """Decrypt the user's stored private key and return a ready SSHConnection."""
    private_key = decrypt_private_key(current_user['private_key_encrypted'])
    return make_user_connection(
        current_user['cluster_host'],
        current_user['cluster_username'],
        private_key,
    )


def _config_path(home_dir: str) -> str:
    """Remote path to the user's BSP config file."""
    return f'{home_dir}/.config/bioinformatics-tools/config.yaml'


@router.get("/workflows")
def list_workflows(current_user: dict = Depends(get_current_user)):
    """Return the list of user-facing workflows with detailed metadata."""
    return _get_available_workflows()


@router.get("/health")
def health_check():
    """Test endpoint to verify API is working. No auth required."""
    return {"status": "success"}


@router.get("/status")
def ssh_status(current_user: dict = Depends(get_current_user)):
    """Check whether the BSP server can reach the user's cluster via SSH."""
    try:
        conn = _build_connection(current_user)
        ssh = conn.connect()
        ssh.close()
        return {"connected": True, "host": current_user["cluster_host"]}
    except Exception as exc:
        LOGGER.warning("SSH status check failed for user %s: %s", current_user["username"], exc)
        return {"connected": False, "host": current_user["cluster_host"]}


@router.get("/config")
def get_config(current_user: dict = Depends(get_current_user)):
    """Read the user's ~/.config/bioinformatics-tools/config.yaml from their cluster via SFTP."""
    conn = _build_connection(current_user)
    path = _config_path(current_user["home_dir"])
    try:
        data = ssh_sftp.read_remote_yaml(path, connection=conn)
        return data
    except Exception as exc:
        LOGGER.error("Failed to read remote config for %s: %s", current_user["username"], exc)
        raise HTTPException(status_code=500, detail=f"Failed to read remote config: {exc}")


@router.put("/config")
def save_config(config: dict, current_user: dict = Depends(get_current_user)):
    """Write a config dict back to the user's cluster as YAML via SFTP."""
    conn = _build_connection(current_user)
    path = _config_path(current_user["home_dir"])
    try:
        ssh_sftp.write_remote_yaml(path, config, connection=conn)
        return {"success": True}
    except Exception as exc:
        LOGGER.error("Failed to write remote config for %s: %s", current_user["username"], exc)
        raise HTTPException(status_code=500, detail=f"Failed to write remote config: {exc}")


@router.post("/config/create-default")
def create_default_config(current_user: dict = Depends(get_current_user)):
    """Create a default config file with all system defaults populated."""
    conn = _build_connection(current_user)
    path = _config_path(current_user["home_dir"])

    # Build default config from REQUIRED_SYSTEM_PARAMS
    default_config = {
        "main_database": "~/.local/share/bioinformatics-tools/my-db.db",
        "compute": {
            "cluster_default": {}
        }
    }

    # Populate compute.cluster_default with all defaults from REQUIRED_SYSTEM_PARAMS
    for param in REQUIRED_SYSTEM_PARAMS:
        if param['param'].startswith('compute.cluster_default.'):
            key = param['param'].split('.')[-1]  # Extract the last part (e.g., 'account', 'partition')
            default_value = param.get('default')
            # Use empty string for required fields with no default, otherwise use the default
            default_config['compute']['cluster_default'][key] = default_value if default_value is not None else ""

    try:
        ssh_sftp.write_remote_yaml(path, default_config, connection=conn)
        LOGGER.info("Created default config for user %s at %s", current_user["username"], path)
        return {"success": True, "config": default_config}
    except Exception as exc:
        LOGGER.error("Failed to create default config for %s: %s", current_user["username"], exc)
        raise HTTPException(status_code=500, detail=f"Failed to create default config: {exc}")


@router.post("/test-path-writable")
def test_path_writable(path_data: dict, current_user: dict = Depends(get_current_user)):
    """Test if a path on the cluster is writable by attempting to create parent directories and a test file."""
    conn = _build_connection(current_user)
    test_path = path_data.get("path", "").strip()

    if not test_path:
        raise HTTPException(status_code=400, detail="Path is required")

    try:
        ssh = conn.connect()

        # Expand ~ to actual home directory
        if test_path.startswith("~"):
            test_path = test_path.replace("~", current_user["home_dir"], 1)

        # Get the directory (remove filename if present)
        import posixpath
        test_dir = posixpath.dirname(test_path)

        # Try to create the directory structure
        _, stdout, stderr = ssh.exec_command(f'mkdir -p "{test_dir}" 2>&1 && echo "DIR_OK"')
        output = stdout.read().decode().strip()

        if "DIR_OK" not in output:
            ssh.close()
            return {
                "writable": False,
                "error": f"Cannot create directory: {test_dir}",
                "details": output
            }

        # Try to write a test file
        test_file = f"{test_path}.write_test"
        _, stdout, stderr = ssh.exec_command(f'touch "{test_file}" 2>&1 && rm -f "{test_file}" 2>&1 && echo "WRITE_OK"')
        output = stdout.read().decode().strip()

        ssh.close()

        if "WRITE_OK" in output:
            return {"writable": True}
        else:
            return {
                "writable": False,
                "error": f"Path is not writable: {test_path}",
                "details": output
            }

    except Exception as exc:
        LOGGER.error("Failed to test path writability for %s: %s", current_user["username"], exc)
        return {
            "writable": False,
            "error": f"Failed to test path: {str(exc)}"
        }


@router.post("/run_slurm")
def run_slurm(content: SlurmSend, current_user: dict = Depends(get_current_user)):
    """Submit a SLURM job and return the job ID immediately."""
    conn = _build_connection(current_user)
    job_id = ssh_slurm.submit_slurm_job(script_content=content.script, connection=conn)
    return {"success": True, "job_id": job_id, "message": "Job submitted successfully"}


@router.post("/run_ssh")
def run_ssh(content: SlurmSend, current_user: dict = Depends(get_current_user)):
    """Execute an SSH command and return output."""
    LOGGER.info('Running run_ssh for user %s', current_user["username"])
    conn = _build_connection(current_user)
    std_txt = ssh_slurm.submit_ssh_job(cmd=content.script, connection=conn)
    return {"success": True, "std_txt": std_txt, "message": "Job submitted successfully"}


def _check_genome_path_exists(genome_path: str, workflow: str, conn) -> None:
    """Raise HTTPException(400) if genome_path doesn't exist on the cluster,
    or (for batch-input workflows) the folder has no recognized genome file.

    Shared by run_workflow (fresh submissions) and resume_job/restart_job
    (re-validating a previously-known-good path, since the remote file
    could have been deleted or moved since the original run)."""
    wf_key = WORKFLOWS.get(workflow)
    supports_batch_input = bool(wf_key and wf_key.supports_batch_input)
    try:
        if supports_batch_input:
            attr = ssh_sftp.check_remote_path_kind(genome_path, conn)
            if attr == 'directory':
                entries = ssh_sftp.list_remote_dir(genome_path, conn)
                has_genome_file = any(
                    e['type'] == 'file' and e['name'].lower().endswith(GENOME_EXTENSIONS)
                    for e in entries
                )
                if not has_genome_file:
                    raise HTTPException(
                        status_code=400,
                        detail=f"No recognized genome files (e.g. .fasta, .fa, .fna) found in folder: '{genome_path}'",
                    )
        else:
            ssh_sftp.check_remote_file(genome_path, conn)
    except FileNotFoundError:
        raise HTTPException(
            status_code=400,
            detail=f"Path not found on the cluster: '{genome_path}'. "
                   "Make sure the path is a Negishi path, not a path on your local machine.",
        )
    except IsADirectoryError:
        raise HTTPException(
            status_code=400,
            detail=f"Path points to a directory, not a file: '{genome_path}'",
        )
    except HTTPException:
        raise
    except Exception as exc:
        LOGGER.warning("File pre-check failed for %s: %s", genome_path, exc)
        raise HTTPException(
            status_code=400,
            detail=f"Could not verify path on cluster: {exc}",
        )


def _launch_job(
    *, genome_path: str, workflow: str, base_output_dir: str,
    selected_tools: list[str] | None, current_user: dict, conn, main_db: str | None,
    relaunched_from: str | None = None,
    copy_from_work_dir: str | None = None,
) -> dict:
    """Shared job-launch sequence: generate job_id/timestamp/output_dir,
    optionally copy a previous run's output_dir into the new one first
    (Resume), persist to job_store, build the dane_wf command (adding
    margie_sb.resume: true when copy_from_work_dir is set, so Snakemake's
    mtime-only rerun triggers recognize the copied-forward outputs as
    already done -- see workflow_tools/workflow.py's build_executable), and
    submit.

    Does NOT do workflow-id/stub/config/genome-path pre-flight validation --
    callers (run_workflow, resume_job, restart_job) each do whatever subset
    of that is appropriate for their entry point before calling this.
    """
    selected_tools_csv = ",".join(selected_tools) if selected_tools is not None else None
    selected_tools_arg = f" {workflow}.selected_tools: {selected_tools_csv}" if selected_tools_csv else ""

    job_id = str(uuid.uuid4())
    timestamp = datetime.now().strftime('%Y-%m-%d-%H%M')
    output_dir = f"{base_output_dir.rstrip('/')}/{timestamp}"

    if copy_from_work_dir:
        ssh_sftp.copy_remote_directory(copy_from_work_dir, output_dir, connection=conn)
        try:
            ssh_sftp.rewrite_path_references(output_dir, copy_from_work_dir, output_dir, connection=conn)
        except Exception as exc:
            # Cosmetic cleanup only (confirmed: nothing downstream reads the
            # stale provenance columns this fixes up) -- never let a failure
            # here block the resumed job from launching.
            LOGGER.warning("Could not rewrite stale path references for resumed job: %s", exc)

    job_store.create(
        job_id, genome_path, user_id=current_user["user_id"],
        workflow=workflow, output_dir=output_dir,
        selected_tools=selected_tools_csv, relaunched_from=relaunched_from,
        persist_db_path=main_db, persist_connection=conn,
    )
    job_store.update(job_id, work_dir=output_dir)

    # The CLI dispatcher (caragols) matches do_<a>_<b> against the SEPARATE
    # tokens "<a> <b>", not the underscore-joined string -- e.g. do_margie_sb
    # is invoked as "margie sb", not "margie_sb". workflow itself (and every
    # config key built from it, e.g. selected_tools_arg/resume_arg) stays
    # underscore-joined, matching the registry's cmd_identifier.
    dispatch_tokens = workflow.replace('_', ' ')
    resume_arg = f" {workflow}.resume: true" if copy_from_work_dir else ""
    # Invoke dane_wf directly from ~/bioinformatics-tools/.venv (an editable
    # install -- code changes there are picked up instantly, no reinstall
    # needed) rather than through `uvx --from`. uvx re-resolves/caches the
    # local package on its own schedule, independent of whether
    # ~/bioinformatics-tools (see _ensure_remote_deployment_symlink in
    # api/main.py) actually points at fresh code -- confirmed adding 3+
    # minutes per run and still serving a stale build without --refresh.
    # The venv binary has neither problem: ~0.4s overhead, always current.
    command = (
        f"~/bioinformatics-tools/.venv/bin/dane_wf {dispatch_tokens}"
        f" input: {genome_path} output_dir: {output_dir}{selected_tools_arg}{resume_arg}"
    )
    job_runner.submit_job(job_id, command, connection=conn)

    return {"success": True, "job_id": job_id, "output_dir": output_dir, "message": "Job submitted successfully"}


@router.post("/run_workflow")
def run_workflow(genome_data: GenomeSend, current_user: dict = Depends(get_current_user)):
    """Submit a genome analysis workflow by name."""
    available_workflows = _get_available_workflows()
    allowed_ids = {wf["id"] for wf in available_workflows}

    if genome_data.workflow not in allowed_ids:
        raise HTTPException(status_code=400, detail=f"Unknown workflow '{genome_data.workflow}'. Available: {sorted(allowed_ids)}")

    if genome_data.workflow in STUB_WORKFLOWS:
        raise HTTPException(status_code=501, detail=f"Workflow '{genome_data.workflow}' is not yet implemented. Check back soon!")

    conn = _build_connection(current_user)

    # Pre-flight: validate required config values are set
    config_path = _config_path(current_user["home_dir"])
    try:
        user_config = ssh_sftp.read_remote_yaml(config_path, connection=conn)
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="Configuration file not found. Please create a configuration in your Profile settings first."
        )

    # Validate required fields
    missing_fields = []

    # Check main_database
    main_db = user_config.get('main_database')
    if not main_db or str(main_db).strip() == '':
        missing_fields.append('main_database')

    # Check compute.cluster_default.account
    account = user_config.get('compute', {}).get('cluster_default', {}).get('account')
    if not account or str(account).strip() == '':
        missing_fields.append('compute.cluster_default.account (SLURM account)')

    if missing_fields:
        raise HTTPException(
            status_code=400,
            detail=f"Required configuration missing: {', '.join(missing_fields)}. "
                   "Please configure these in your Profile settings before running workflows."
        )

    # Resolve genome path / output dir, falling back to the user's global config
    # defaults (input_path / output_path) when the request didn't specify one.
    genome_path = genome_data.genome_path or user_config.get(genome_data.workflow, {}).get('input_path')
    if not genome_path or str(genome_path).strip() == '':
        raise HTTPException(
            status_code=400,
            detail="No genome file or folder specified, and no input_path default is configured. "
                   "Set one in your Profile settings or pass genome_path explicitly.",
        )

    _check_genome_path_exists(genome_path, genome_data.workflow, conn)

    # Validate phase/tool selection, if given -- catch typos here rather than
    # have them silently no-op as an unrecognized run_<tool> config key.
    # `is not None` (not a truthiness check) matters here: an explicit empty
    # list means "run nothing", which must NOT be treated the same as
    # omitting the field entirely ("run everything").
    if genome_data.selected_tools is not None:
        if not genome_data.selected_tools:
            raise HTTPException(
                status_code=400,
                detail="selected_tools was empty -- select at least one tool/phase to run, "
                       "or omit the field entirely to run everything.",
            )
        valid_tool_keys = {tool['key'] for tool in MARGIE_SB_PHASED_TOOLS}
        unknown = set(genome_data.selected_tools) - valid_tool_keys
        if unknown:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown tool key(s) in selected_tools: {sorted(unknown)}. "
                       f"Available: {sorted(valid_tool_keys)}",
            )

    base_dir = (genome_data.output_dir or user_config.get(genome_data.workflow, {}).get('output_path') or current_user['home_dir']).rstrip('/')

    return _launch_job(
        genome_path=genome_path, workflow=genome_data.workflow,
        base_output_dir=base_dir, selected_tools=genome_data.selected_tools,
        current_user=current_user, conn=conn, main_db=main_db,
    )


def _job_from_history_row(row: dict) -> dict:
    """Shape a persisted api_jobs row like an in-memory job_store entry, so
    the front-end's job page can render it the same way whether the job is
    still live or was resumed after a dane-api restart.

    logs/slurm_jobs/containers are only ever a final snapshot, taken once
    by job_store.finalize() at job completion/failure -- never persisted
    incrementally (that would mean an SSH round-trip per log line), so a
    job that's still mid-run when dane-api restarts has none of this yet,
    and a job that died without ever reaching finalize() (e.g. dane-api
    itself crashed) never gets one at all. sub_jobs/report/progress/
    steps_done/total remain pure live-session detail, never persisted.
    job_history_client already JSON-decodes slurm_jobs/containers back
    into real lists before this row ever reaches here.
    """
    return {
        "job_id": row["job_id"],
        "status": row["status"],
        "phase": row.get("phase"),
        "genome_path": row.get("genome_path"),
        "workflow": row.get("workflow"),
        "work_dir": row.get("work_dir"),
        "selected_tools": row.get("selected_tools"),
        "relaunched_from": row.get("relaunched_from"),
        "start_time": row.get("created_at"),
        "sub_jobs": [],
        "slurm_jobs": row.get("slurm_jobs") or [],
        "containers": row.get("containers") or [],
        "logs": row.get("logs") or "",
        "resumed_from_history": True,
    }


def _load_job_for_action(job_id: str, current_user: dict, conn) -> dict:
    """Resolve a job_id to enough info to relaunch it (genome_path, workflow,
    work_dir, selected_tools, status), whether it's still live in job_store
    or only in persisted history. Raises 404/403 the same way get_job_status
    does.

    Deliberately separate from get_job_status itself: that endpoint also
    does SLURM-reconciliation (still_active/status_note) this lookup
    doesn't need, and from _resolve_job_work_dir (file-serving endpoints
    only need work_dir; resume/restart need the full launch-relevant
    shape). Minor duplication of the job_store-then-history-fallback
    pattern across these three is accepted -- collapsing them would couple
    endpoints with different response-shape needs."""
    job = job_store.get(job_id)
    if job is not None:
        if job.get("user_id") != current_user["user_id"]:
            raise HTTPException(status_code=403, detail="Access denied")
        return job

    try:
        user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
    except Exception:
        raise HTTPException(status_code=404, detail="Job not found")
    main_db = user_config.get('main_database')
    row = job_history_client.get_job(conn, main_db, job_id) if main_db else None
    if row is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return _job_from_history_row(row)


_TERMINAL_STATUSES = {"completed", "failed", "cancelled"}
_POTENTIALLY_STALE_STATUSES = {"pending", "running", "snakemake"}


@router.get("/job_status/{job_id}")
def get_job_status(job_id: str, current_user: dict = Depends(get_current_user)):
    """Get status of a running job. Falls back to persistent history (e.g.
    after a dane-api restart wiped the in-memory job_store) before giving
    up. Returns 403 if the in-memory job belongs to a different user --
    history rows can't make that check since they live in each user's own
    main_database already, with no cross-user data to separate.

    For a non-terminal history row (nothing live is tracking it anymore,
    so its persisted status could be stale), also checks squeue to see
    whether it's genuinely still active on the cluster -- adds
    still_active/status_note to the response, purely additive, so this
    never changes the existing status field or breaks older clients."""
    job = job_store.get(job_id)
    if job is not None:
        if job.get("user_id") != current_user["user_id"]:
            raise HTTPException(status_code=403, detail="Access denied")
        return {**job, "cluster_host": current_user["cluster_host"]}

    conn = _build_connection(current_user)
    try:
        user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
    except Exception:
        raise HTTPException(status_code=404, detail="Job not found")

    main_db = user_config.get('main_database')
    row = job_history_client.get_job(conn, main_db, job_id) if main_db else None
    if row is None:
        raise HTTPException(status_code=404, detail="Job not found")

    result = {**_job_from_history_row(row), "cluster_host": current_user["cluster_host"]}

    if row["status"] not in _TERMINAL_STATUSES and row.get("work_dir"):
        try:
            matches = ssh_slurm.find_active_jobs_in_workdir(
                row["work_dir"], current_user["cluster_username"], connection=conn,
            )
        except Exception as exc:
            LOGGER.warning("SLURM reconciliation check failed for job %s: %s", job_id, exc)
            matches = []
        still_active = any(m["state"] in ("RUNNING", "PENDING") for m in matches)
        result["still_active"] = still_active
        if still_active:
            m = matches[0]
            result["status_note"] = (
                f"Still {m['state']} on the cluster (SLURM job {m['job_id']}) -- "
                "history shows a stale status from before a server restart."
            )
        else:
            result["status_note"] = (
                f"No longer active on the cluster -- last recorded phase: {row.get('phase') or 'unknown'}."
            )

    return result


@router.get("/jobs")
def list_jobs(
    workflow: str | None = None,
    page: int = Query(1, ge=1),
    page_size: int = Query(20, ge=1, le=100),
    current_user: dict = Depends(get_current_user),
):
    """List this user's persistent job history, optionally filtered to one
    workflow, paginated most-recent-first. A brand new user with no history
    at all gets an empty list, not an error."""
    empty_response = {"jobs": [], "page": page, "page_size": page_size, "total_jobs": 0, "total_pages": 1}

    conn = _build_connection(current_user)
    try:
        user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
    except Exception:
        return empty_response

    main_db = user_config.get('main_database')
    if not main_db:
        return empty_response

    offset = (page - 1) * page_size
    rows = job_history_client.list_jobs(conn, main_db, workflow=workflow, limit=page_size, offset=offset)
    total_jobs = job_history_client.count_jobs(conn, main_db, workflow=workflow)
    total_pages = max((total_jobs + page_size - 1) // page_size, 1)

    return {
        "jobs": [_job_from_history_row(row) for row in rows],
        "page": page,
        "page_size": page_size,
        "total_jobs": total_jobs,
        "total_pages": total_pages,
    }


@router.post("/cancel_job/{job_id}")
def cancel_job(job_id: str, current_user: dict = Depends(get_current_user)):
    """Emergency stop - cancel all SLURM jobs, kill remote process, and mark job as cancelled.

    Falls back to persistent history when the job isn't in the in-memory
    job_store -- e.g. after a dane-api restart -- the same fallback
    _resolve_job_work_dir already uses, so Emergency Stop still works for
    jobs the status page can still show via that fallback.
    """
    conn = _build_connection(current_user)
    job = job_store.get(job_id)
    main_db = None

    if job is not None:
        if job.get("user_id") != current_user["user_id"]:
            raise HTTPException(status_code=403, detail="Access denied")
        slurm_ids = [sj["job_id"] for sj in job_store.get_slurm_jobs(job_id)]
    else:
        try:
            user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
        except Exception:
            raise HTTPException(status_code=404, detail="Job not found")

        main_db = user_config.get('main_database')
        row = job_history_client.get_job(conn, main_db, job_id) if main_db else None
        if row is None:
            raise HTTPException(status_code=404, detail="Job not found")
        slurm_ids = [sj["job_id"] for sj in (row.get("slurm_jobs") or [])]

    # Cancel all SLURM subjobs
    if slurm_ids:
        ssh_slurm.cancel_slurm_jobs(slurm_ids, connection=conn)
        LOGGER.info("Cancelled %d SLURM jobs for job %s", len(slurm_ids), job_id)

    # Kill the remote dane_wf process on the login node
    # This ensures the SSH task stops immediately instead of waiting for Snakemake to notice
    ssh_slurm.kill_remote_process("dane_wf", connection=conn)
    LOGGER.info("Killed remote dane_wf process for job %s", job_id)

    # Mark job as cancelled (this will also stop the status checker daemon, if any)
    if job is not None:
        job_store.cancel(job_id)
    else:
        job_history_client.record_job_updated(
            conn, main_db, job_id, status="cancelled", phase="Cancelled by user",
        )

    return {
        "success": True,
        "message": f"Cancelled job {job_id}",
        "slurm_jobs_cancelled": len(slurm_ids)
    }


def _main_db_for(current_user: dict, conn) -> str:
    """Fetch main_database from the user's config, raising HTTPException(400)
    if the config or that field is missing. Shared pre-flight for resume_job/
    restart_job (run_workflow does the same check inline as part of its
    richer missing_fields validation, which these two intentionally don't
    repeat in full -- they're relaunching an already-known-good prior
    submission, not validating a fresh one)."""
    try:
        user_config = ssh_sftp.read_remote_yaml(_config_path(current_user["home_dir"]), connection=conn)
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="Configuration file not found. Please create a configuration in your Profile settings first.",
        )
    main_db = user_config.get('main_database')
    if not main_db or str(main_db).strip() == '':
        raise HTTPException(
            status_code=400,
            detail="main_database is not configured. Please configure it in your Profile settings.",
        )
    return main_db


def _base_output_dir_from(path: str) -> str:
    """Strip a trailing /YYYY-MM-DD-HHMM timestamp segment off a previous
    job's work_dir/output_dir, recovering the base directory run_workflow
    originally built it from (see its own timestamp = ...strftime('%Y-%m-%d-%H%M'))."""
    return re.sub(r'/\d{4}-\d{2}-\d{2}-\d{4}$', '', path)


@router.post("/resume_job/{job_id}")
def resume_job(job_id: str, current_user: dict = Depends(get_current_user)):
    """Resume a failed (or stale-but-no-longer-active) job: copy its
    work_dir into a new timestamped folder, then relaunch the same
    genome_path/workflow/selected_tools there with margie_sb.resume: true
    so Snakemake's mtime-based rebuild skips whatever already completed."""
    conn = _build_connection(current_user)
    original = _load_job_for_action(job_id, current_user, conn)

    status = original.get("status")
    work_dir = original.get("work_dir")
    if status not in ("failed", "cancelled"):
        if status not in _POTENTIALLY_STALE_STATUSES:
            raise HTTPException(
                status_code=400,
                detail=f"Cannot resume a job with status '{status}' -- only failed/cancelled jobs, "
                       "or non-terminal jobs no longer actually active on the cluster, can be resumed.",
            )
        # status looks non-terminal (pending/running/snakemake) -- this
        # could be a stale label left over from before a dane-api restart
        # (see get_job_status's still_active reconciliation), so actually
        # check the cluster before trusting it.
        still_active = True
        if work_dir:
            try:
                matches = ssh_slurm.find_active_jobs_in_workdir(
                    work_dir, current_user["cluster_username"], connection=conn,
                )
                still_active = any(m["state"] in ("RUNNING", "PENDING") for m in matches)
            except Exception as exc:
                LOGGER.warning("SLURM active-check failed while resuming job %s: %s", job_id, exc)
                still_active = True  # fail safe: don't resume something we couldn't confirm is dead
        if still_active:
            raise HTTPException(
                status_code=400,
                detail=f"Cannot resume a job with status '{status}' while it appears to still be "
                       "active on the cluster.",
            )

    if not work_dir:
        raise HTTPException(status_code=400, detail="Original job has no recorded working directory to resume from")

    genome_path = original.get("genome_path")
    workflow = original.get("workflow")
    if not genome_path or not workflow:
        raise HTTPException(status_code=400, detail="Original job is missing genome_path/workflow, cannot resume")

    main_db = _main_db_for(current_user, conn)
    _check_genome_path_exists(genome_path, workflow, conn)

    selected_tools_csv = original.get("selected_tools")
    selected_tools = selected_tools_csv.split(",") if selected_tools_csv else None

    return _launch_job(
        genome_path=genome_path, workflow=workflow,
        base_output_dir=_base_output_dir_from(work_dir), selected_tools=selected_tools,
        current_user=current_user, conn=conn, main_db=main_db,
        relaunched_from=job_id, copy_from_work_dir=work_dir,
    )


@router.post("/restart_job/{job_id}")
def restart_job(job_id: str, current_user: dict = Depends(get_current_user)):
    """Restart a job from scratch: same genome_path/workflow/selected_tools
    as a brand-new job/timestamp, with no copying -- functionally "Start
    New Analysis" auto-filled from a prior job. Works from any status."""
    conn = _build_connection(current_user)
    original = _load_job_for_action(job_id, current_user, conn)

    base_for_dir = original.get("work_dir") or original.get("output_dir")
    if not base_for_dir:
        raise HTTPException(status_code=400, detail="Original job has no recorded output directory to restart from")

    genome_path = original.get("genome_path")
    workflow = original.get("workflow")
    if not genome_path or not workflow:
        raise HTTPException(status_code=400, detail="Original job is missing genome_path/workflow, cannot restart")

    main_db = _main_db_for(current_user, conn)
    _check_genome_path_exists(genome_path, workflow, conn)

    selected_tools_csv = original.get("selected_tools")
    selected_tools = selected_tools_csv.split(",") if selected_tools_csv else None

    return _launch_job(
        genome_path=genome_path, workflow=workflow,
        base_output_dir=_base_output_dir_from(base_for_dir), selected_tools=selected_tools,
        current_user=current_user, conn=conn, main_db=main_db,
        relaunched_from=job_id,
    )


@router.get("/job_files/{job_id}")
def get_job_files(
    job_id: str,
    subdir: str = "",
    current_user: dict = Depends(get_current_user),
):
    """List output files for a job via SFTP."""
    _validate_relative_path(subdir, label="subdirectory")

    conn = _build_connection(current_user)
    work_dir = _resolve_job_work_dir(job_id, current_user, conn)

    target_dir = f"{work_dir}/{subdir}".rstrip("/") if subdir else work_dir

    try:
        entries = ssh_sftp.list_remote_dir(target_dir, connection=conn)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Directory not found on remote")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to list remote directory: {str(e)}")

    return {"work_dir": work_dir, "subdir": subdir, "entries": entries}


@router.get("/download_file/{job_id}")
def download_file(
    job_id: str,
    path: str,
    format: str = Query("raw", pattern="^(raw|excel)$"),
    current_user: dict = Depends(get_current_user),
):
    """Download a file from a job's working directory via SFTP.

    format=raw (default) streams the file unmodified. format=excel reads
    the whole delimited text file into memory and converts it to .xlsx
    before sending -- only sensible for the TSV/CSV outputs the viewer
    already supports, not arbitrary binary files.
    """
    _validate_relative_path(path)

    conn = _build_connection(current_user)
    work_dir = _resolve_job_work_dir(job_id, current_user, conn)

    remote_path = f"{work_dir}/{path}"
    filename = path.split("/")[-1]

    if format == "excel":
        try:
            content = b"".join(ssh_sftp.stream_remote_file(remote_path, connection=conn))
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail="File not found on remote")
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Failed to download file: {str(e)}")

        text = content.decode("utf-8", errors="replace")
        delimiter = _detect_delimiter(path, text.splitlines()[0] if text else "")
        try:
            df = pd.read_csv(io.StringIO(text), sep=delimiter, dtype=str, keep_default_na=False)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Could not parse file as a table: {str(e)}")

        buffer = io.BytesIO()
        df.to_excel(buffer, index=False, engine="openpyxl")
        buffer.seek(0)

        xlsx_filename = re.sub(r"\.(tsv|csv)$", "", filename, flags=re.IGNORECASE) + ".xlsx"
        return StreamingResponse(
            buffer,
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={"Content-Disposition": f'attachment; filename="{xlsx_filename}"'},
        )

    try:
        return StreamingResponse(
            ssh_sftp.stream_remote_file(remote_path, connection=conn),
            media_type="application/octet-stream",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="File not found on remote")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to download file: {str(e)}")


@router.get("/view_file/{job_id}")
def view_file(
    job_id: str,
    path: str,
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=500),
    current_user: dict = Depends(get_current_user),
):
    """Read a paginated slice of a delimited text file from a job's
    working directory, for in-browser viewing without downloading the
    whole file. Cost scales with page position, not file size -- see
    ssh_sftp.read_remote_file_page."""
    _validate_relative_path(path)

    conn = _build_connection(current_user)
    work_dir = _resolve_job_work_dir(job_id, current_user, conn)

    remote_path = f"{work_dir}/{path}"

    cache_key = (job_id, path)
    known_total_lines = None
    try:
        mtime, size = ssh_sftp.stat_remote_file(remote_path, connection=conn)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="File not found on remote")

    cached = _line_count_cache.get(cache_key)
    if cached and cached[0] == mtime and cached[1] == size:
        known_total_lines = cached[2]

    start_row = 2 + (page - 1) * page_size  # row 1 is always the header
    end_row = start_row + page_size - 1

    try:
        result = ssh_sftp.read_remote_file_page(
            remote_path, start_row, end_row, connection=conn,
            known_total_lines=known_total_lines,
        )
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="File not found on remote")
    except IsADirectoryError:
        raise HTTPException(status_code=400, detail="Path is a directory, not a file")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to read file: {str(e)}")

    total_lines = result["total_lines"]
    if known_total_lines is None:
        _line_count_cache[cache_key] = (mtime, size, total_lines)
    total_rows = max(total_lines - 1, 0)
    total_pages = max((total_rows + page_size - 1) // page_size, 1)

    delimiter = _detect_delimiter(path, result["header"])
    columns = result["header"].split(delimiter) if result["header"] else []
    rows = [line.split(delimiter) for line in result["lines"]]

    return {
        "columns": columns,
        "rows": rows,
        "page": page,
        "page_size": page_size,
        "total_rows": total_rows,
        "total_pages": total_pages,
    }


@router.get("/job_status/{job_id}/stream")
def stream_job_status(job_id: str, current_user: dict = Depends(get_current_user)):
    """SSE endpoint that streams real-time job status updates."""
    job = job_store.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.get("user_id") != current_user["user_id"]:
        raise HTTPException(status_code=403, detail="Access denied")

    return StreamingResponse(
        job_runner.job_status_generator(job_id),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
        }
    )


@router.get("/all_genomes")
def all_genomes(path: str, current_user: dict = Depends(get_current_user)):
    """List genome files at a remote path on the user's cluster."""
    conn = _build_connection(current_user)
    genomes = ssh_slurm.get_genomes(path, connection=conn)
    return {"success": True, "Genomes": genomes}
