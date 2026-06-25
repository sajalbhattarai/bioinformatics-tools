"""
Pytest-based API tests using FastAPI TestClient.

No live server required — TestClient handles requests in-process.

Tiers:
  1. Health/root endpoints (no mocks)
  2. JobStore unit tests (in-memory, no mocks)
  3. SSH endpoints (mocked external calls)
  4. Path traversal security tests
  5. Auth endpoints (register / login / me)
"""
import io
import uuid
from unittest.mock import MagicMock, patch

import paramiko
import pytest
from fastapi.testclient import TestClient

from bioinformatics_tools.api.auth import get_current_user
from bioinformatics_tools.api.main import app
from bioinformatics_tools.api.services.job_store import job_store


# ---------------------------------------------------------------------------
# Shared test helpers
# ---------------------------------------------------------------------------

#: Fake user dict returned by the get_current_user dependency override.
#: user_id=1 must match what we pass to job_store.create(user_id=...) in
#: tests that exercise ownership checks.
FAKE_USER = {
    "user_id": 1,
    "username": "testuser",
    "cluster_host": "test.cluster.edu",
    "cluster_username": "testuser",
    "home_dir": "/home/testuser",
    "private_key_encrypted": "fake_encrypted_key",
}


@pytest.fixture(scope="session")
def test_rsa_key():
    """Generate a 1024-bit RSA key once per session for auth tests."""
    key = paramiko.RSAKey.generate(1024)
    buf = io.StringIO()
    key.write_private_key(buf)
    return buf.getvalue()


@pytest.fixture()
def client(tmp_path, monkeypatch):
    """
    Fresh TestClient with an isolated SQLite DB and cleared job store.

    monkeypatch sets BSP_DB_PATH to a per-test temp dir so that the startup
    event creates a clean users table each time.
    """
    monkeypatch.setenv("BSP_DB_PATH", str(tmp_path / "test.db"))
    job_store._jobs.clear()
    with TestClient(app) as c:
        yield c
    job_store._jobs.clear()


@pytest.fixture()
def authed_client(client):
    """
    client with get_current_user dependency overridden to return FAKE_USER.

    Use this for any test that hits an authenticated endpoint but isn't
    specifically testing the auth flow itself.
    """
    app.dependency_overrides[get_current_user] = lambda: FAKE_USER
    yield client
    app.dependency_overrides.pop(get_current_user, None)


# ---------------------------------------------------------------------------
# Tier 1 — Health / root endpoints (no mocks needed)
# ---------------------------------------------------------------------------

class TestHealthEndpoints:
    """Smoke tests for every health and info endpoint."""

    def test_root(self, client):
        resp = client.get("/")
        assert resp.status_code == 200
        body = resp.json()
        assert body["status"] == "success"
        assert "endpoints" in body

    def test_health(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "success"

    def test_fasta_health(self, client):
        resp = client.get("/v1/fasta/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "success"

    def test_ssh_health(self, client):
        resp = client.get("/v1/ssh/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "success"


# ---------------------------------------------------------------------------
# Tier 2 — JobStore unit tests (purely in-memory)
# ---------------------------------------------------------------------------

class TestJobStore:
    """Direct tests of the JobStore singleton (no HTTP involved)."""

    def setup_method(self):
        job_store._jobs.clear()

    def teardown_method(self):
        job_store._jobs.clear()

    def test_create(self):
        job = job_store.create("j1", "/genomes/ecoli.fasta")
        assert job["job_id"] == "j1"
        assert job["status"] == "pending"
        assert job["genome_path"] == "/genomes/ecoli.fasta"
        assert job["work_dir"] is None
        assert "start_time" in job

    def test_create_with_user_id(self):
        job = job_store.create("j1b", "/genomes/ecoli.fasta", user_id=42)
        assert job["user_id"] == 42

    def test_get_missing(self):
        assert job_store.get("nonexistent") is None

    def test_update(self):
        job_store.create("j2", "/g")
        job_store.update("j2", status="running", phase="Aligning")
        job = job_store.get("j2")
        assert job["status"] == "running"
        assert job["phase"] == "Aligning"

    def test_append_log(self):
        job_store.create("j3", "/g")
        job_store.append_log("j3", "line one")
        job_store.append_log("j3", "line two")
        logs = job_store.get("j3")["logs"]
        assert "line one" in logs
        assert "line two" in logs

    def test_add_slurm_job(self):
        job_store.create("j4", "/g")
        job_store.add_slurm_job("j4", slurm_id="12345", rule="fastp")
        slurm_jobs = job_store.get_slurm_jobs("j4")
        assert len(slurm_jobs) == 1
        assert slurm_jobs[0]["job_id"] == "12345"
        assert slurm_jobs[0]["rule"] == "fastp"
        assert slurm_jobs[0]["status"] == "SUBMITTED"

    def test_finalize_updates_in_memory_status_and_phase(self):
        job_store.create("j5", "/g")
        job_store.finalize("j5", status="completed", phase="Done")
        job = job_store.get("j5")
        assert job["status"] == "completed"
        assert job["phase"] == "Done"

    def test_finalize_on_missing_job_is_a_noop(self):
        job_store.finalize("nonexistent-job", status="completed", phase="Done")  # must not raise

    @patch("bioinformatics_tools.api.services.job_store.job_history_client")
    def test_finalize_persists_full_snapshot(self, mock_history):
        mock_conn = MagicMock()
        job_store.create(
            "j6", "/g", persist_db_path="~/my-db.db", persist_connection=mock_conn,
        )
        job_store.append_log("j6", "line one")
        job_store.add_slurm_job("j6", slurm_id="111", rule="quast")
        job_store.add_container("j6", {"name": "quast", "version": "5.0"})

        job_store.finalize("j6", status="completed", phase="Done")

        mock_history.record_job_updated.assert_called_once()
        _, kwargs = mock_history.record_job_updated.call_args
        assert kwargs["status"] == "completed"
        assert kwargs["phase"] == "Done"
        assert "line one" in kwargs["logs"]
        assert kwargs["slurm_jobs"][0]["job_id"] == "111"
        assert kwargs["containers"][0]["name"] == "quast"

    @patch("bioinformatics_tools.api.services.job_store.job_history_client")
    def test_finalize_without_persistence_configured_skips_history_call(self, mock_history):
        job_store.create("j7", "/g")  # no persist_db_path/persist_connection
        job_store.finalize("j7", status="failed", phase="Error")
        mock_history.record_job_updated.assert_not_called()


class TestJobStatusEndpoint:
    """HTTP-level tests for job_status (requires auth, enforces ownership)."""

    def test_job_status_endpoint(self, authed_client):
        # user_id must match FAKE_USER["user_id"] (1) to pass ownership check
        job_store.create("test-job-1", "/genomes/test.fasta", user_id=1)
        resp = authed_client.get("/v1/ssh/job_status/test-job-1")
        assert resp.status_code == 200
        body = resp.json()
        assert body["job_id"] == "test-job-1"
        assert body["status"] == "pending"

    def test_job_status_404(self, authed_client):
        resp = authed_client.get("/v1/ssh/job_status/nonexistent")
        assert resp.status_code == 404

    def test_job_status_ownership_denied(self, authed_client):
        # Job belongs to a different user — should get 403
        job_store.create("other-users-job", "/genomes/test.fasta", user_id=999)
        resp = authed_client.get("/v1/ssh/job_status/other-users-job")
        assert resp.status_code == 403

    def test_job_status_requires_auth(self, client):
        resp = client.get("/v1/ssh/job_status/any-job")
        assert resp.status_code == 401


class TestJobStatusHistoryReconciliation:
    """get_job_status's history-fallback branch (job not in job_store, e.g.
    after a dane-api restart) -- for non-terminal statuses, checks squeue
    to see whether the job is genuinely still active on the cluster."""

    @staticmethod
    def _history_row(status="running", work_dir="/scratch/x/2026-06-21-1118"):
        return {
            "job_id": "resumed-job", "status": status, "phase": "Running Snakemake",
            "work_dir": work_dir, "workflow": "margie_sb",
            "genome_path": "/g/e.fasta", "created_at": "2026-06-21T11:18:00Z",
        }

    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_still_active_true_when_squeue_matches(
        self, mock_build_conn, mock_sftp, mock_history, mock_slurm, authed_client,
    ):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_slurm.find_active_jobs_in_workdir.return_value = [{"job_id": "39600517", "state": "RUNNING"}]

        resp = authed_client.get("/v1/ssh/job_status/resumed-job")
        assert resp.status_code == 200
        body = resp.json()
        assert body["status"] == "running"  # unchanged -- additive only
        assert body["still_active"] is True
        assert "39600517" in body["status_note"]

    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_still_active_false_when_no_squeue_match(
        self, mock_build_conn, mock_sftp, mock_history, mock_slurm, authed_client,
    ):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_slurm.find_active_jobs_in_workdir.return_value = []

        resp = authed_client.get("/v1/ssh/job_status/resumed-job")
        body = resp.json()
        assert body["status"] == "running"  # status field itself untouched
        assert body["still_active"] is False
        assert "Running Snakemake" in body["status_note"]

    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_terminal_status_skips_squeue_call_entirely(
        self, mock_build_conn, mock_sftp, mock_history, mock_slurm, authed_client,
    ):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row(status="completed")

        resp = authed_client.get("/v1/ssh/job_status/resumed-job")
        body = resp.json()
        assert body["status"] == "completed"
        assert "still_active" not in body
        mock_slurm.find_active_jobs_in_workdir.assert_not_called()

    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_squeue_failure_does_not_break_endpoint(
        self, mock_build_conn, mock_sftp, mock_history, mock_slurm, authed_client,
    ):
        """A transient SSH/squeue failure during reconciliation must not
        turn an otherwise-200 history-fallback response into a 500."""
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_slurm.find_active_jobs_in_workdir.side_effect = Exception("ssh timeout")

        resp = authed_client.get("/v1/ssh/job_status/resumed-job")
        assert resp.status_code == 200
        assert resp.json()["still_active"] is False

    def test_live_job_in_job_store_unaffected(self, authed_client):
        """Confirms the in-memory branch (job_store hit) never reaches the
        new SLURM-reconciliation code at all -- no still_active/status_note
        key should appear on a normal in-memory job response."""
        job_store.create("live-job-1", "/genomes/test.fasta", user_id=1)
        job_store.update("live-job-1", status="running")
        resp = authed_client.get("/v1/ssh/job_status/live-job-1")
        body = resp.json()
        assert body["status"] == "running"
        assert "still_active" not in body
        assert "status_note" not in body


class TestFileEndpointsHistoryFallback:
    """job_files/download_file/view_file must work for a job resumed from
    history (not in job_store), not just live in-memory jobs -- previously
    all three 404'd unconditionally for any job dane-api wasn't actively
    tracking, even though job_status already resolved a valid work_dir for
    the same job via history fallback."""

    @staticmethod
    def _history_row(work_dir="/scratch/x/2026-06-21-1118"):
        return {
            "job_id": "resumed-job", "status": "completed", "phase": "Done",
            "work_dir": work_dir, "workflow": "margie_sb",
            "genome_path": "/g/e.fasta", "created_at": "2026-06-21T11:18:00Z",
        }

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_job_files_resolves_via_history(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_sftp.list_remote_dir.return_value = [{"name": "results.tsv", "type": "file", "size": 123}]

        resp = authed_client.get("/v1/ssh/job_files/resumed-job")
        assert resp.status_code == 200
        body = resp.json()
        assert body["work_dir"] == "/scratch/x/2026-06-21-1118"
        assert body["entries"][0]["name"] == "results.tsv"

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_download_file_resolves_via_history(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_sftp.stream_remote_file.return_value = iter([b"chunk"])

        resp = authed_client.get("/v1/ssh/download_file/resumed-job", params={"path": "results.tsv"})
        assert resp.status_code == 200

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_view_file_resolves_via_history(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = self._history_row()
        mock_sftp.stat_remote_file.return_value = (1700000000.0, 1234)
        mock_sftp.read_remote_file_page.return_value = {
            "total_lines": 2, "header": "a\tb", "lines": ["1\t2"],
        }

        resp = authed_client.get("/v1/ssh/view_file/resumed-job", params={"path": "results.tsv"})
        assert resp.status_code == 200
        assert resp.json()["columns"] == ["a", "b"]

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_job_files_404_when_not_in_store_or_history(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = None

        resp = authed_client.get("/v1/ssh/job_files/nonexistent-job")
        assert resp.status_code == 404

    def test_path_traversal_still_checked_before_any_lookup(self, authed_client):
        """Reordering validation ahead of the job/connection resolution
        (needed so history fallback can build a connection) must not weaken
        the traversal guard -- it should still reject before touching SSH
        at all, with no mocks needed since nothing real should be called."""
        resp = authed_client.get(
            "/v1/ssh/job_files/any-job", params={"subdir": "../../etc"},
        )
        assert resp.status_code == 400


# ---------------------------------------------------------------------------
# Tier 3 — SSH endpoints with mocks
# ---------------------------------------------------------------------------

class TestSSHEndpointsMocked:
    """Endpoints that call external SSH/SLURM services — all mocked."""

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.job_runner")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_run_margie(self, mock_sftp, mock_runner, mock_build_conn, authed_client):
        # Mock the config file read and genome file check
        mock_sftp.read_remote_yaml.return_value = {
            "main_database": "~/my-db.db",
            "compute": {"cluster_default": {"account": "test-account"}}
        }
        mock_sftp.check_remote_file.return_value = None  # No exception = file exists

        resp = authed_client.post(
            "/v1/ssh/run_workflow",
            json={"genome_path": "/depot/genomes/ecoli.fasta"},
        )
        assert resp.status_code == 200
        body = resp.json()
        assert body["success"] is True
        assert "job_id" in body

        # Verify job_runner.submit_job was called with the created job_id
        mock_runner.submit_job.assert_called_once()
        call_args = mock_runner.submit_job.call_args
        assert call_args[0][0] == body["job_id"]

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    def test_run_slurm(self, mock_slurm, mock_build_conn, authed_client):
        mock_slurm.submit_slurm_job.return_value = "99999"
        resp = authed_client.post(
            "/v1/ssh/run_slurm",
            json={"script": "#!/bin/bash\necho hello"},
        )
        assert resp.status_code == 200
        body = resp.json()
        assert body["success"] is True
        assert body["job_id"] == "99999"
        mock_slurm.submit_slurm_job.assert_called_once()

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    def test_all_genomes(self, mock_slurm, mock_build_conn, authed_client):
        mock_slurm.get_genomes.return_value = ["genome1.fasta", "genome2.fasta"]
        resp = authed_client.get("/v1/ssh/all_genomes", params={"path": "/depot/genomes"})
        assert resp.status_code == 200
        body = resp.json()
        assert body["success"] is True
        assert len(body["Genomes"]) == 2
        mock_slurm.get_genomes.assert_called_once()

    def test_all_genomes_requires_auth(self, client):
        resp = client.get("/v1/ssh/all_genomes", params={"path": "/depot/genomes"})
        assert resp.status_code == 401


class TestJobsHistoryPagination:
    """GET /v1/ssh/jobs -- paginated job history listing."""

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_response_includes_pagination_fields(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.list_jobs.return_value = [
            {"job_id": "j1", "status": "completed", "workflow": "margie_sb", "created_at": "t1"},
        ]
        mock_history.count_jobs.return_value = 45

        resp = authed_client.get("/v1/ssh/jobs", params={"page": 2, "page_size": 20})
        assert resp.status_code == 200
        body = resp.json()
        assert body["page"] == 2
        assert body["page_size"] == 20
        assert body["total_jobs"] == 45
        assert body["total_pages"] == 3  # ceil(45 / 20)
        assert len(body["jobs"]) == 1

    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_offset_computed_from_page(self, mock_build_conn, mock_sftp, mock_history, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.list_jobs.return_value = []
        mock_history.count_jobs.return_value = 0

        authed_client.get("/v1/ssh/jobs", params={"page": 3, "page_size": 20})

        _, kwargs = mock_history.list_jobs.call_args
        assert kwargs["offset"] == 40  # (3 - 1) * 20
        assert kwargs["limit"] == 20

    def test_page_size_bounds(self, authed_client):
        resp = authed_client.get("/v1/ssh/jobs", params={"page_size": 9999})
        assert resp.status_code == 422

    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    def test_empty_when_no_config(self, mock_build_conn, mock_sftp, authed_client):
        mock_sftp.read_remote_yaml.side_effect = Exception("no config")
        resp = authed_client.get("/v1/ssh/jobs")
        assert resp.status_code == 200
        body = resp.json()
        assert body["jobs"] == []
        assert body["total_pages"] == 1


# ---------------------------------------------------------------------------
# Tier 4 — Path traversal security tests
# ---------------------------------------------------------------------------

class TestResumeAndRestartJob:
    """/v1/ssh/resume_job and /v1/ssh/restart_job -- relaunch a prior job,
    with or without copying its work_dir forward."""

    def _create_job(self, status="failed", work_dir="/remote/work/2026-06-21-1118", selected_tools=None):
        jid = str(uuid.uuid4())
        job_store.create(
            jid, "/genomes/test.fasta", user_id=1, workflow="margie_sb",
            selected_tools=selected_tools,
        )
        job_store.update(jid, status=status, work_dir=work_dir)
        return jid

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.job_runner")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_resume_job_copies_and_relaunches(self, mock_sftp, mock_runner, mock_build_conn, authed_client):
        jid = self._create_job(status="failed", selected_tools="quast,gtdbtk")
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}

        resp = authed_client.post(f"/v1/ssh/resume_job/{jid}")
        assert resp.status_code == 200
        body = resp.json()
        assert body["success"] is True
        assert body["job_id"] != jid

        mock_sftp.copy_remote_directory.assert_called_once()
        copy_args = mock_sftp.copy_remote_directory.call_args[0]
        assert copy_args[0] == "/remote/work/2026-06-21-1118"

        mock_runner.submit_job.assert_called_once()
        command = mock_runner.submit_job.call_args[0][1]
        assert "margie_sb.resume: true" in command
        assert "margie_sb.selected_tools: quast,gtdbtk" in command

        new_job = job_store.get(body["job_id"])
        assert new_job["relaunched_from"] == jid
        assert new_job["selected_tools"] == "quast,gtdbtk"

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_resume_job_rejects_completed_status(self, mock_sftp, mock_build_conn, authed_client):
        jid = self._create_job(status="completed")
        resp = authed_client.post(f"/v1/ssh/resume_job/{jid}")
        assert resp.status_code == 400
        mock_sftp.copy_remote_directory.assert_not_called()

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.job_runner")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_resume_job_allows_stale_non_terminal_when_not_active(
        self, mock_sftp, mock_runner, mock_slurm, mock_build_conn, authed_client,
    ):
        """A job stuck showing 'running' from before a dane-api restart
        should still be resumable once confirmed not actually active."""
        jid = self._create_job(status="running")
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_slurm.find_active_jobs_in_workdir.return_value = []

        resp = authed_client.post(f"/v1/ssh/resume_job/{jid}")
        assert resp.status_code == 200

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_resume_job_rejects_when_still_genuinely_active(
        self, mock_sftp, mock_slurm, mock_build_conn, authed_client,
    ):
        jid = self._create_job(status="running")
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_slurm.find_active_jobs_in_workdir.return_value = [{"job_id": "1", "state": "RUNNING"}]

        resp = authed_client.post(f"/v1/ssh/resume_job/{jid}")
        assert resp.status_code == 400
        mock_sftp.copy_remote_directory.assert_not_called()

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.job_runner")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_restart_job_works_from_completed_status_with_no_copy(
        self, mock_sftp, mock_runner, mock_build_conn, authed_client,
    ):
        jid = self._create_job(status="completed", selected_tools="quast")
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}

        resp = authed_client.post(f"/v1/ssh/restart_job/{jid}")
        assert resp.status_code == 200
        body = resp.json()
        assert body["job_id"] != jid

        mock_sftp.copy_remote_directory.assert_not_called()
        command = mock_runner.submit_job.call_args[0][1]
        assert "margie_sb.resume: true" not in command
        assert "margie_sb.selected_tools: quast" in command

        new_job = job_store.get(body["job_id"])
        assert new_job["relaunched_from"] == jid

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.job_history_client")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_resume_job_404_for_unknown_job(self, mock_sftp, mock_history, mock_build_conn, authed_client):
        mock_sftp.read_remote_yaml.return_value = {"main_database": "~/my-db.db"}
        mock_history.get_job.return_value = None
        resp = authed_client.post("/v1/ssh/resume_job/nonexistent-job-id")
        assert resp.status_code == 404


class TestPathTraversalSecurity:
    """Ensure path-based endpoints reject directory traversal attempts."""

    def _create_job_with_workdir(self):
        jid = str(uuid.uuid4())
        # user_id=1 matches FAKE_USER so ownership check passes
        job_store.create(jid, "/genomes/test.fasta", user_id=1)
        job_store.update(jid, work_dir="/remote/work/dir")
        return jid

    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_job_files_path_traversal(self, mock_sftp, authed_client):
        jid = self._create_job_with_workdir()
        resp = authed_client.get(
            f"/v1/ssh/job_files/{jid}",
            params={"subdir": "../../etc"},
        )
        assert resp.status_code == 400
        assert "Invalid" in resp.json()["detail"]
        mock_sftp.list_remote_dir.assert_not_called()

    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_download_file_path_traversal(self, mock_sftp, authed_client):
        jid = self._create_job_with_workdir()
        resp = authed_client.get(
            f"/v1/ssh/download_file/{jid}",
            params={"path": "../../etc/passwd"},
        )
        assert resp.status_code == 400
        assert "Invalid" in resp.json()["detail"]
        mock_sftp.stream_remote_file.assert_not_called()

    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_view_file_path_traversal(self, mock_sftp, authed_client):
        jid = self._create_job_with_workdir()
        resp = authed_client.get(
            f"/v1/ssh/view_file/{jid}",
            params={"path": "../../etc/passwd"},
        )
        assert resp.status_code == 400
        assert "Invalid" in resp.json()["detail"]
        mock_sftp.read_remote_file_page.assert_not_called()
        mock_sftp.stat_remote_file.assert_not_called()


# ---------------------------------------------------------------------------
# Tier 4b — Paginated file viewer (view_file)
# ---------------------------------------------------------------------------

class TestViewFile:
    """/v1/ssh/view_file -- paginated read of a job's output file."""

    def _create_job_with_workdir(self):
        jid = str(uuid.uuid4())
        job_store.create(jid, "/genomes/test.fasta", user_id=1)
        job_store.update(jid, work_dir="/remote/work/dir")
        return jid

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_view_file_happy_path(self, mock_sftp, mock_build_conn, authed_client):
        jid = self._create_job_with_workdir()
        mock_sftp.stat_remote_file.return_value = (1700000000.0, 1234)
        mock_sftp.read_remote_file_page.return_value = {
            "total_lines": 5,  # header + 4 data rows
            "header": "gene_id\tphase",
            "lines": ["g1\t3", "g2\t4"],
        }

        resp = authed_client.get(
            f"/v1/ssh/view_file/{jid}",
            params={"path": "results.tsv", "page": 1, "page_size": 2},
        )
        assert resp.status_code == 200
        body = resp.json()
        assert body["columns"] == ["gene_id", "phase"]
        assert body["rows"] == [["g1", "3"], ["g2", "4"]]
        assert body["total_rows"] == 4
        assert body["total_pages"] == 2
        assert body["page"] == 1
        assert body["page_size"] == 2

        # Cache was empty, so the first call must not claim a known total.
        _, kwargs = mock_sftp.read_remote_file_page.call_args
        assert kwargs["known_total_lines"] is None

    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_sftp")
    def test_view_file_reuses_cached_total_lines(self, mock_sftp, mock_build_conn, authed_client):
        """Second page request for the same unchanged file should skip the
        wc -l pass entirely -- this is what keeps paging through a huge
        file from re-paying a full scan on every click."""
        jid = self._create_job_with_workdir()
        mock_sftp.stat_remote_file.return_value = (1700000000.0, 1234)
        mock_sftp.read_remote_file_page.return_value = {
            "total_lines": 5,
            "header": "gene_id\tphase",
            "lines": ["g1\t3"],
        }

        resp1 = authed_client.get(
            f"/v1/ssh/view_file/{jid}",
            params={"path": "results.tsv", "page": 1, "page_size": 1},
        )
        assert resp1.status_code == 200
        assert mock_sftp.read_remote_file_page.call_args.kwargs["known_total_lines"] is None

        resp2 = authed_client.get(
            f"/v1/ssh/view_file/{jid}",
            params={"path": "results.tsv", "page": 2, "page_size": 1},
        )
        assert resp2.status_code == 200
        assert mock_sftp.read_remote_file_page.call_args.kwargs["known_total_lines"] == 5

    def test_view_file_page_size_bounds(self, authed_client):
        # FastAPI's Query(..., le=500) must reject this before the handler
        # body (and therefore any SSH call) ever runs -- this is what
        # actually prevents a client from forcing a full-file sed range.
        resp = authed_client.get(
            "/v1/ssh/view_file/any-job-id",
            params={"path": "results.tsv", "page_size": 99999},
        )
        assert resp.status_code == 422


# ---------------------------------------------------------------------------
# Tier 5 — Auth endpoints (register / login / me)
# ---------------------------------------------------------------------------

class TestAuth:
    """
    End-to-end tests for the auth flow.

    Registration requires a valid SSH private key and a reachable cluster.
    We mock make_user_connection at the router level so no real SSH happens.
    """

    #: Base registration payload — private_key added per-test from the
    #: session-scoped test_rsa_key fixture.
    BASE_REG = {
        "username": "authtest",
        "password": "S3cur3P@ss!",
        "cluster_host": "test.cluster.edu",
        "cluster_username": "authtest",
    }

    def _mock_ssh_conn(self, home_dir: str = "/home/authtest"):
        """Return a MagicMock SSHConnection whose connect() returns a stub SSH client."""
        mock_stdout = MagicMock()
        mock_stdout.read.return_value = home_dir.encode() + b"\n"
        mock_ssh = MagicMock()
        mock_ssh.exec_command.return_value = (None, mock_stdout, None)
        mock_conn = MagicMock()
        mock_conn.connect.return_value = mock_ssh
        return mock_conn

    # --- register ---

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_register_success(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value = self._mock_ssh_conn()
        resp = client.post(
            "/v1/auth/register",
            json={**self.BASE_REG, "private_key": test_rsa_key},
        )
        assert resp.status_code == 201
        body = resp.json()
        assert body["username"] == "authtest"
        assert "user_id" in body

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_register_duplicate_username(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value = self._mock_ssh_conn()
        data = {**self.BASE_REG, "private_key": test_rsa_key}
        client.post("/v1/auth/register", json=data)  # first succeeds
        resp = client.post("/v1/auth/register", json=data)  # duplicate fails
        assert resp.status_code == 400
        assert "taken" in resp.json()["detail"].lower()

    def test_register_bad_private_key(self, client):
        resp = client.post(
            "/v1/auth/register",
            json={**self.BASE_REG, "private_key": "this-is-not-a-valid-ssh-key"},
        )
        assert resp.status_code == 400
        assert "parse" in resp.json()["detail"].lower()

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_register_ssh_connection_fails(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value.connect.side_effect = Exception("Connection refused")
        resp = client.post(
            "/v1/auth/register",
            json={**self.BASE_REG, "private_key": test_rsa_key},
        )
        assert resp.status_code == 400
        # Error message mentions the cluster host
        assert "test.cluster.edu" in resp.json()["detail"]

    # --- login ---

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_login_success(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value = self._mock_ssh_conn()
        client.post("/v1/auth/register", json={**self.BASE_REG, "private_key": test_rsa_key})

        resp = client.post(
            "/v1/auth/login",
            json={"username": "authtest", "password": "S3cur3P@ss!"},
        )
        assert resp.status_code == 200
        body = resp.json()
        assert "access_token" in body
        assert body["token_type"] == "bearer"

    def test_login_wrong_password(self, client):
        resp = client.post(
            "/v1/auth/login",
            json={"username": "nonexistent", "password": "wrong"},
        )
        assert resp.status_code == 401
        assert resp.json()["detail"] == "Invalid credentials"

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_login_wrong_password_for_real_user(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value = self._mock_ssh_conn()
        client.post("/v1/auth/register", json={**self.BASE_REG, "private_key": test_rsa_key})

        resp = client.post(
            "/v1/auth/login",
            json={"username": "authtest", "password": "wrongpassword"},
        )
        assert resp.status_code == 401
        assert resp.json()["detail"] == "Invalid credentials"

    # --- /me ---

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    def test_me_with_valid_token(self, mock_make_conn, client, test_rsa_key):
        mock_make_conn.return_value = self._mock_ssh_conn()
        client.post("/v1/auth/register", json={**self.BASE_REG, "private_key": test_rsa_key})
        login_resp = client.post(
            "/v1/auth/login",
            json={"username": "authtest", "password": "S3cur3P@ss!"},
        )
        token = login_resp.json()["access_token"]

        resp = client.get("/v1/auth/me", headers={"Authorization": f"Bearer {token}"})
        assert resp.status_code == 200
        body = resp.json()
        assert body["username"] == "authtest"
        assert body["cluster_host"] == "test.cluster.edu"
        assert body["home_dir"] == "/home/authtest"
        # Sensitive fields must never be exposed
        assert "password_hash" not in body
        assert "private_key_encrypted" not in body
        assert "private_key" not in body

    def test_me_without_token(self, client):
        resp = client.get("/v1/auth/me")
        assert resp.status_code == 401

    def test_me_with_invalid_token(self, client):
        resp = client.get("/v1/auth/me", headers={"Authorization": "Bearer not-a-real-token"})
        assert resp.status_code == 401

    # --- protected endpoints require auth ---

    def test_protected_endpoint_rejects_no_token(self, client):
        resp = client.get("/v1/ssh/all_genomes", params={"path": "/depot/genomes"})
        assert resp.status_code == 401

    @patch("bioinformatics_tools.api.routers.auth.make_user_connection")
    @patch("bioinformatics_tools.api.routers.ssh._build_connection")
    @patch("bioinformatics_tools.api.routers.ssh.ssh_slurm")
    def test_protected_endpoint_with_real_token(
        self, mock_slurm, mock_build_conn, mock_make_conn, client, test_rsa_key
    ):
        """Full round-trip: register → login → hit a protected endpoint."""
        mock_make_conn.return_value = self._mock_ssh_conn()
        mock_slurm.get_genomes.return_value = ["genome1.fasta"]

        client.post("/v1/auth/register", json={**self.BASE_REG, "private_key": test_rsa_key})
        login_resp = client.post(
            "/v1/auth/login",
            json={"username": "authtest", "password": "S3cur3P@ss!"},
        )
        token = login_resp.json()["access_token"]

        resp = client.get(
            "/v1/ssh/all_genomes",
            params={"path": "/depot/genomes"},
            headers={"Authorization": f"Bearer {token}"},
        )
        assert resp.status_code == 200
        assert resp.json()["Genomes"] == ["genome1.fasta"]
