"""
Unit tests for bioinformatics_tools.utilities.ssh_slurm.

These test the SSH-exec helpers directly (no FastAPI TestClient, no real
network) by mocking the SSHConnection's paramiko client.
"""
from unittest.mock import MagicMock

from bioinformatics_tools.utilities import ssh_slurm


def _mock_connection_for_exec(exec_stdout: bytes):
    mock_ssh = MagicMock()
    mock_stdout = MagicMock()
    mock_stdout.read.return_value = exec_stdout
    mock_ssh.exec_command.return_value = (None, mock_stdout, MagicMock())

    mock_connection = MagicMock()
    mock_connection.connect.return_value = mock_ssh
    return mock_connection, mock_ssh


class TestFindActiveJobsInWorkdir:
    def test_matches_running_job(self):
        mock_connection, _ = _mock_connection_for_exec(
            b"39600517|RUNNING|/scratch/negishi/u/margie/output/2026-06-21-1118\n"
        )
        result = ssh_slurm.find_active_jobs_in_workdir(
            "/scratch/negishi/u/margie/output/2026-06-21-1118", "u", connection=mock_connection,
        )
        assert result == [{"job_id": "39600517", "state": "RUNNING"}]

    def test_trailing_slash_normalized_both_sides(self):
        mock_connection, _ = _mock_connection_for_exec(b"1|RUNNING|/scratch/x/job/\n")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == [{"job_id": "1", "state": "RUNNING"}]

    def test_no_match_returns_empty(self):
        mock_connection, _ = _mock_connection_for_exec(b"1|RUNNING|/scratch/x/other_job\n")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == []

    def test_empty_queue_returns_empty(self):
        mock_connection, _ = _mock_connection_for_exec(b"")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == []

    def test_multiple_jobs_only_matching_workdir_returned(self):
        mock_connection, _ = _mock_connection_for_exec(
            b"1|RUNNING|/scratch/x/job\n2|PENDING|/scratch/x/other\n3|RUNNING|/scratch/x/job\n"
        )
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == [
            {"job_id": "1", "state": "RUNNING"},
            {"job_id": "3", "state": "RUNNING"},
        ]
