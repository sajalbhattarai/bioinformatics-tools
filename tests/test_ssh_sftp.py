"""
Unit tests for bioinformatics_tools.utilities.ssh_sftp, with the SSHConnection's
paramiko client mocked (no network).
"""
import shlex
import stat as stat_module
import subprocess
import sys
from unittest.mock import MagicMock

import pytest

from bioinformatics_tools.utilities import ssh_sftp


def _mock_connection_for_exec(exec_stdout: bytes, exec_stderr: bytes = b""):
    """Builds an SSHConnection mock whose SFTP stat() reports a regular file and
    whose exec_command() returns the given output."""
    mock_ssh = MagicMock()

    mock_sftp_client = MagicMock()
    mock_stat_result = MagicMock()
    mock_stat_result.st_mode = stat_module.S_IFREG | 0o644  # regular file
    mock_sftp_client.stat.return_value = mock_stat_result
    mock_ssh.open_sftp.return_value = mock_sftp_client

    mock_stdout = MagicMock()
    mock_stdout.read.return_value = exec_stdout
    mock_stderr = MagicMock()
    mock_stderr.read.return_value = exec_stderr
    mock_ssh.exec_command.return_value = (None, mock_stdout, mock_stderr)

    mock_connection = MagicMock()
    mock_connection.connect.return_value = mock_ssh
    return mock_connection, mock_ssh


class TestReadRemoteFilePage:
    def test_escapes_path_with_shell_metacharacters(self):
        """Checks that a path with spaces and parentheses is shell-quoted in the remote command."""
        sentinel = ssh_sftp._PAGE_SENTINEL
        stdout = (
            f"3\n{sentinel}\ncol_a\tcol_b\n{sentinel}\nval1\tval2\n"
        ).encode()
        dangerous_path = "/work/dir/interpro results (final).tsv"
        mock_connection, mock_ssh = _mock_connection_for_exec(stdout)

        result = ssh_sftp.read_remote_file_page(
            dangerous_path, start_row=2, end_row=2, connection=mock_connection,
        )

        command = mock_ssh.exec_command.call_args[0][0]
        quoted = shlex.quote(dangerous_path)
        # One quoted use each for wc, the header sed and the range sed.
        assert command.count(quoted) == 3

        assert result == {
            "total_lines": 3,
            "header": "col_a\tcol_b",
            "lines": ["val1\tval2"],
        }

    def test_skips_wc_when_total_lines_known(self):
        """Checks that known_total_lines skips the wc -l pass."""
        sentinel = ssh_sftp._PAGE_SENTINEL
        stdout = f"col_a\tcol_b\n{sentinel}\nval1\tval2\n".encode()
        mock_connection, mock_ssh = _mock_connection_for_exec(stdout)

        result = ssh_sftp.read_remote_file_page(
            "/work/dir/results.tsv", start_row=2, end_row=2,
            connection=mock_connection, known_total_lines=42,
        )

        command = mock_ssh.exec_command.call_args[0][0]
        assert "wc -l" not in command
        assert result["total_lines"] == 42
        assert result["header"] == "col_a\tcol_b"
        assert result["lines"] == ["val1\tval2"]

    def test_print_before_quit_keeps_last_row_of_page(self):
        """Checks the sed range is `START,ENDp;ENDq`, so the last row prints before sed quits."""
        sentinel = ssh_sftp._PAGE_SENTINEL
        stdout = f"42\n{sentinel}\nh\n{sentinel}\nrow2\trow2b\nrow3\trow3b\n".encode()
        mock_connection, mock_ssh = _mock_connection_for_exec(stdout)

        ssh_sftp.read_remote_file_page(
            "/work/dir/results.tsv", start_row=2, end_row=3, connection=mock_connection,
        )

        command = mock_ssh.exec_command.call_args[0][0]
        assert "2,3p;3q" in command

    def test_empty_page_past_end_of_file(self):
        """Checks that a page beyond EOF returns an empty page."""
        sentinel = ssh_sftp._PAGE_SENTINEL
        stdout = f"3\n{sentinel}\ncol_a\n{sentinel}\n".encode()
        mock_connection, _ = _mock_connection_for_exec(stdout)

        result = ssh_sftp.read_remote_file_page(
            "/work/dir/results.tsv", start_row=100, end_row=199, connection=mock_connection,
        )

        assert result["total_lines"] == 3
        assert result["lines"] == []


def _mock_connection_for_copy(exit_code: int, stderr: bytes = b""):
    mock_ssh = MagicMock()
    mock_stdout = MagicMock()
    mock_stdout.channel.recv_exit_status.return_value = exit_code
    mock_stderr = MagicMock()
    mock_stderr.read.return_value = stderr
    mock_ssh.exec_command.return_value = (None, mock_stdout, mock_stderr)

    mock_connection = MagicMock()
    mock_connection.connect.return_value = mock_ssh
    return mock_connection, mock_ssh


class TestCopyRemoteDirectory:
    def test_builds_rsync_with_excludes_and_succeeds(self):
        mock_connection, mock_ssh = _mock_connection_for_copy(exit_code=0)

        ssh_sftp.copy_remote_directory("/old/run", "/new/run", connection=mock_connection)

        command = mock_ssh.exec_command.call_args[0][0]
        assert "mkdir -p" in command
        assert "rsync -a" in command
        assert "--exclude='.snakemake'" in command
        assert "--exclude='original_container_outputs/*/stage'" in command

    def test_raises_on_nonzero_exit(self):
        mock_connection, _ = _mock_connection_for_copy(exit_code=1, stderr=b"rsync error")

        with pytest.raises(RuntimeError, match="rsync failed"):
            ssh_sftp.copy_remote_directory("/old/run", "/new/run", connection=mock_connection)

    def test_escapes_paths_with_spaces(self):
        mock_connection, mock_ssh = _mock_connection_for_copy(exit_code=0)

        ssh_sftp.copy_remote_directory("/old/run with spaces", "/new run", connection=mock_connection)

        command = mock_ssh.exec_command.call_args[0][0]
        assert shlex.quote("/old/run with spaces/") in command
        assert shlex.quote("/new run") in command


class TestRewritePathReferences:
    """Tests the script that rewrites old-run paths in copied provenance columns after a Resume copy."""

    def test_script_rewrites_old_path_in_text_files(self, tmp_path):
        """Runs the generated script as a subprocess against fixture files."""
        old_dir = str(tmp_path / "2026-06-21-1118")
        new_dir = str(tmp_path / "2026-06-21-1300")

        nested = tmp_path / "2026-06-21-1118" / "genome" / "kegg"
        nested.mkdir(parents=True)
        results_file = nested / "kegg_results.tsv"
        results_file.write_text(
            f"input_path\toutput_path\n{old_dir}/genome/rasttk/rast.faa\t{old_dir}/genome/kegg\n"
        )
        untouched_file = nested / "no_path_here.tsv"
        untouched_file.write_text("col_a\tcol_b\n1\t2\n")

        script = ssh_sftp._build_path_rewrite_script(old_dir, old_dir, new_dir)
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)

        assert result.returncode == 0
        assert result.stdout.strip() == "1"  # only kegg_results.tsv contained old_dir
        assert new_dir in results_file.read_text()
        assert old_dir not in results_file.read_text()
        assert untouched_file.read_text() == "col_a\tcol_b\n1\t2\n"  # unmodified, no match

    def test_script_skips_binary_files_without_crashing(self, tmp_path):
        target_dir = tmp_path / "run"
        target_dir.mkdir()
        binary_file = target_dir / "blob.bin"
        binary_file.write_bytes(b"\xff\xfe\x00\x01not valid utf-8 \xfa")

        script = ssh_sftp._build_path_rewrite_script(str(target_dir), str(tmp_path / "old"), str(tmp_path / "new"))
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)

        assert result.returncode == 0
        assert result.stdout.strip() == "0"

    def test_script_handles_paths_with_dots_and_dashes_literally(self, tmp_path):
        """Checks that '.' and '-' in paths are matched literally."""
        old_dir = str(tmp_path / "2026-06-21-1118")
        new_dir = str(tmp_path / "2026-06-21-9999")
        target = tmp_path / "data"
        target.mkdir()
        f = target / "result.tsv"
        f.write_text(f"{old_dir}/Genome_GCF_000027325.1/rasttk/rast.faa\n")

        script = ssh_sftp._build_path_rewrite_script(str(target), old_dir, new_dir)
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)

        assert result.returncode == 0
        assert f.read_text() == f"{new_dir}/Genome_GCF_000027325.1/rasttk/rast.faa\n"

    def _mock_connection_for_rewrite(self, stdout_text: str, exit_code: int = 0, stderr: bytes = b""):
        mock_ssh = MagicMock()
        mock_stdout = MagicMock()
        mock_stdout.read.return_value = stdout_text.encode()
        mock_stdout.channel.recv_exit_status.return_value = exit_code
        mock_stderr = MagicMock()
        mock_stderr.read.return_value = stderr
        mock_ssh.exec_command.return_value = (None, mock_stdout, mock_stderr)
        mock_connection = MagicMock()
        mock_connection.connect.return_value = mock_ssh
        return mock_connection, mock_ssh

    def test_returns_modified_count_from_stdout(self):
        mock_connection, _ = self._mock_connection_for_rewrite("3\n")

        result = ssh_sftp.rewrite_path_references("/new/dir", "/old/dir", "/new/dir", connection=mock_connection)

        assert result == 3

    def test_returns_zero_on_nonzero_exit_without_raising(self):
        """Checks that a failed rewrite returns 0 instead of raising."""
        mock_connection, _ = self._mock_connection_for_rewrite("", exit_code=1, stderr=b"python error")

        result = ssh_sftp.rewrite_path_references("/new/dir", "/old/dir", "/new/dir", connection=mock_connection)

        assert result == 0

    def test_invokes_remote_venv_python_with_quoted_script(self):
        mock_connection, mock_ssh = self._mock_connection_for_rewrite("0\n")

        ssh_sftp.rewrite_path_references("/new/dir", "/old/dir", "/new/dir", connection=mock_connection)

        command = mock_ssh.exec_command.call_args[0][0]
        assert command.startswith("~/bioinformatics-tools/.venv/bin/python -c ")
