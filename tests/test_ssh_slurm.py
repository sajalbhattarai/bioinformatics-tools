"""
Unit tests for bioinformatics_tools.utilities.ssh_slurm.

These test the SSH-exec helpers directly (no FastAPI TestClient, no real
network) by mocking the SSHConnection's paramiko client.
"""
import pytest
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
    """squeue is asked for "%i|%T|%Z|%M" -- id, state, workdir, elapsed. These
    fixtures carry all four; a three-field line is silently skipped by the
    parser, so omitting elapsed made every one of these assert on []."""

    def test_matches_running_job(self):
        mock_connection, _ = _mock_connection_for_exec(
            b"39600517|RUNNING|/scratch/negishi/u/margie/output/2026-06-21-1118|1:02:03\n"
        )
        result = ssh_slurm.find_active_jobs_in_workdir(
            "/scratch/negishi/u/margie/output/2026-06-21-1118", "u", connection=mock_connection,
        )
        assert result == [{"job_id": "39600517", "state": "RUNNING", "time": "1:02:03"}]

    def test_trailing_slash_normalized_both_sides(self):
        mock_connection, _ = _mock_connection_for_exec(b"1|RUNNING|/scratch/x/job/|0:05\n")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == [{"job_id": "1", "state": "RUNNING", "time": "0:05"}]

    def test_no_match_returns_empty(self):
        mock_connection, _ = _mock_connection_for_exec(b"1|RUNNING|/scratch/x/other_job|0:05\n")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == []

    def test_empty_queue_returns_empty(self):
        mock_connection, _ = _mock_connection_for_exec(b"")
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == []

    def test_multiple_jobs_only_matching_workdir_returned(self):
        mock_connection, _ = _mock_connection_for_exec(
            b"1|RUNNING|/scratch/x/job|0:01\n"
            b"2|PENDING|/scratch/x/other|0:00\n"
            b"3|RUNNING|/scratch/x/job|0:03\n"
        )
        result = ssh_slurm.find_active_jobs_in_workdir("/scratch/x/job", "u", connection=mock_connection)
        assert result == [
            {"job_id": "1", "state": "RUNNING", "time": "0:01"},
            {"job_id": "3", "state": "RUNNING", "time": "0:03"},
        ]


class TestProbeRun:
    """probe_run decides whether a detached run's log is worth tailing after a
    dane-api restart. Getting this wrong in the "still alive" direction parks
    one of job_runner's four workers on a `tail -F` that never returns."""

    def test_live_run(self):
        """Log being written, no exit sentinel."""
        mock_connection, _ = _mock_connection_for_exec(b"6806541|-|178\n")
        probe = ssh_slurm.probe_run("job-1", connection=mock_connection)
        assert probe == {"has_log": True, "exit_code": None, "log_idle": 178.0,
                         "driver_state": None}
        assert ssh_slurm.is_replayable(probe)

    def test_finished_run_is_replayable_however_old(self):
        """The sentinel bounds the replay, so age does not matter: tailing
        ends at the exit code instead of following forever. This is how a run
        that finished while the API was down stops being stuck at 'running'."""
        mock_connection, _ = _mock_connection_for_exec(b"628807|0|40432\n")
        probe = ssh_slurm.probe_run("job-2", connection=mock_connection)
        assert probe == {"has_log": True, "exit_code": "0", "log_idle": 40432.0,
                         "driver_state": None}
        assert ssh_slurm.is_replayable(probe)

    def test_interrupted_run_is_not_replayable(self):
        """No sentinel and a log untouched for 11 hours: the driver died
        without recording an exit code. Nothing will ever arrive."""
        mock_connection, _ = _mock_connection_for_exec(b"627147|-|41231\n")
        probe = ssh_slurm.probe_run("job-3", connection=mock_connection)
        assert probe["exit_code"] is None
        assert not ssh_slurm.is_replayable(probe)

    def test_missing_log_is_not_replayable(self):
        mock_connection, _ = _mock_connection_for_exec(b"0|-|0\n")
        assert not ssh_slurm.is_replayable(ssh_slurm.probe_run("job-4", connection=mock_connection))

    def test_unparseable_output_fails_safe(self):
        mock_connection, _ = _mock_connection_for_exec(b"bash: stat: command not found\n")
        probe = ssh_slurm.probe_run("job-5", connection=mock_connection)
        assert probe == {"has_log": False, "exit_code": None, "log_idle": float("inf"),
                         "driver_state": None}
        assert not ssh_slurm.is_replayable(probe)

    def test_probe_does_not_consult_the_process_table(self):
        """The job_id is necessarily in this probe's own command line (it
        stats that job's log), so any pgrep for it matches the probe itself
        and every run reads as alive forever -- the bracket trick included,
        since the unbracketed copy sits right there in the log path."""
        mock_connection, mock_ssh = _mock_connection_for_exec(b"1|-|1\n")
        ssh_slurm.probe_run("job-6", connection=mock_connection)
        cmd = mock_ssh.exec_command.call_args[0][0]
        assert "pgrep" not in cmd and "ps " not in cmd
        assert "stat -c %Y" in cmd          # liveness read off the log instead

    def test_stale_threshold_brackets_a_snakemake_status_cycle(self):
        """Snakemake logs a status-check cycle every ~30s while it has jobs in
        flight, so the cutoff must sit far above that and far below the age of
        an abandoned run."""
        assert 300 < ssh_slurm.RUN_STALE_AFTER < 3600


class TestLaunchDoesNotWaitForTheRun:
    """The launch must not read the exec channel to EOF.

    `A && B && nohup setsid ... &` backgrounds the whole and-list, so a shell
    sits waiting for the workflow with the exec channel still open on its
    stdout. EOF therefore means "the run has finished" -- hours away. Reading
    to EOF hung the generator before its first yield: no tail, no parsed
    lines, and a job page stuck on "Submitting via SSH" with an empty Slurm
    Jobs table for the whole run.
    """

    @staticmethod
    def _connection_whose_stdout_never_ends(pid_line=b"1330075\n"):
        """stdout yields the pid line, then blocks forever -- exactly what the
        real channel does while the launcher is still alive."""
        import threading

        never = threading.Event()          # never set: stands in for "no EOF"

        mock_stdout = MagicMock()
        mock_stdout.readline.return_value = pid_line.decode()
        mock_stdout.read.side_effect = lambda *a, **k: (never.wait(), b"")[1]

        mock_ssh = MagicMock()
        mock_ssh.exec_command.return_value = (None, mock_stdout, MagicMock())
        mock_connection = MagicMock()
        mock_connection.connect.return_value = mock_ssh
        return mock_connection, mock_stdout

    def test_launch_yields_without_waiting_for_eof(self):
        conn, stdout = self._connection_whose_stdout_never_ends()
        gen = ssh_slurm.submit_ssh_job(cmd="dane_wf margie sb", connection=conn,
                                       job_id="job-1")

        # Before the fix this call never returned.
        assert next(gen) == "__LAUNCHED__"
        stdout.read.assert_not_called()
        gen.close()

    def test_launch_survives_a_pid_that_never_arrives(self):
        """A silent launcher must not fail a run that has already started."""
        import socket

        conn, stdout = self._connection_whose_stdout_never_ends()
        stdout.readline.side_effect = socket.timeout("timed out")

        gen = ssh_slurm.submit_ssh_job(cmd="dane_wf margie sb", connection=conn,
                                       job_id="job-2", in_slurm=False)
        assert next(gen) == "__LAUNCHED__"
        gen.close()

    def test_a_refused_sbatch_is_a_real_failure(self):
        """The opposite of the login-node case, and deliberately so: there,
        the run is already going and a missing pid costs only bookkeeping.
        Here nothing has started, so yielding __LAUNCHED__ would leave the
        runner tailing a log that will never exist."""
        conn, stdout = self._connection_whose_stdout_never_ends()
        stdout.readline.return_value = "sbatch: error: Invalid account\n"

        gen = ssh_slurm.submit_ssh_job(cmd="dane_wf margie sb", connection=conn,
                                       job_id="job-2b", in_slurm=True)
        with pytest.raises(RuntimeError, match="Could not submit"):
            next(gen)
        gen.close()

    def test_launch_never_blocks_on_the_exit_status(self):
        """recv_exit_status() waits for the same EOF and is just as fatal."""
        conn, stdout = self._connection_whose_stdout_never_ends()
        gen = ssh_slurm.submit_ssh_job(cmd="dane_wf margie sb", connection=conn,
                                       job_id="job-3")
        next(gen)
        stdout.channel.recv_exit_status.assert_not_called()
        gen.close()


class TestDriverInSlurm:
    """The workflow driver runs as a SLURM job, not on the login node.

    setsid + nohup already survive the SSH session dying, so a closed laptop
    was safe. What they cannot survive is the login node -- a reboot, a
    maintenance window, or an administrator reaping long-running processes,
    which is what Snakemake running for hours looks like. The run then stalls
    half-done: queued SLURM jobs finish, nothing further is ever submitted.
    """

    def _launch(self, cmd="dane_wf margie sb"):
        base = "$HOME/.local/share/bsp/jobs/j1"
        return ssh_slurm.build_driver_launch(
            cmd=cmd, base=base, safe="j1", log=f"{base}.log", rcf=f"{base}.rc",
            jobidf=f"{base}.jobid", driversh=f"{base}.driver.sh")

    def test_it_is_submitted_with_sbatch_and_the_id_recorded(self):
        s = self._launch()
        assert "sbatch --parsable" in s
        assert "$HOME/.local/share/bsp/jobs/j1.jobid" in s

    def test_log_and_sentinel_paths_are_unchanged(self):
        """Everything downstream -- tail -F, probe_run, is_replayable, the
        reattach -- keys off these two files and nothing else."""
        s = self._launch()
        assert "> $HOME/.local/share/bsp/jobs/j1.log 2>&1" in s
        assert "echo $? > $HOME/.local/share/bsp/jobs/j1.rc" in s

    def test_the_workflow_is_not_backgrounded(self):
        """On the login node it had to be, so the SSH call could return. Inside
        a batch job the opposite holds: if the script exits, SLURM tears the
        allocation down and takes the run with it."""
        s = self._launch()
        run_line = [l for l in s.splitlines() if l.startswith("nohup bash ")][0]
        assert not run_line.rstrip().endswith("&")

    def test_nohup_is_used(self):
        assert "nohup bash " in self._launch()

    def test_a_command_containing_quotes_survives(self):
        """The real command carries MARGIE_LICENSE_ACCEPTED='...' -- wrapping
        it in another layer of single quotes would break it, which is why it
        goes into a file of its own via a quoted heredoc."""
        cmd = "MARGIE_LICENSE_ACCEPTED='2026-07-31' MARGIE_USAGE_TYPE='academic' dane_wf margie sb"
        s = self._launch(cmd)
        assert cmd in s
        assert "<<'MARGIE_WORKFLOW_EOF'" in s      # quoted: nothing expands early

    def test_a_queued_driver_is_replayable_even_with_no_log(self):
        """A driver waiting for a node has no log and no sentinel, which under
        the old two-case rule was indistinguishable from a dead run."""
        assert ssh_slurm.is_replayable(
            {"has_log": False, "exit_code": None, "log_idle": float("inf"),
             "driver_state": "PENDING"}) is True

    def test_a_dead_run_with_no_driver_is_still_not_replayable(self):
        assert ssh_slurm.is_replayable(
            {"has_log": False, "exit_code": None, "log_idle": float("inf"),
             "driver_state": None}) is False

    def test_an_interrupted_run_whose_driver_is_gone_is_not_replayable(self):
        assert ssh_slurm.is_replayable(
            {"has_log": True, "exit_code": None, "log_idle": 99999.0,
             "driver_state": None}) is False
