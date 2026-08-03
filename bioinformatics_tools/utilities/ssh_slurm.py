'''
SSH-based SLURM operations for interfacing with HPC clusters.

Uses paramiko to run login node and SLURM commands. Since we use Snakemake,
which controls SLURM batching and queueing, we can mainly run on login node.

All functions are API-layer only. Pass a per-user SSHConnection built
with make_user_connection() for every call.
'''
import logging
import socket
import time
import re
import shlex

from bioinformatics_tools.utilities.ssh_connection import SSHConnection

LOGGER = logging.getLogger(__name__)


def get_genomes(
    location,
    connection: SSHConnection,
):
    """List genome files at a remote path via SSH ls."""
    ssh = connection.connect()
    LOGGER.info('ls -lah %s', location)
    stdin, stdout, stderr = ssh.exec_command(f'ls -lah {location}')
    output = stdout.read().decode()
    error = stderr.read().decode()
    pass  # pooled client: closing it would defeat SSHConnection's pool

    if error:
        LOGGER.warning('Error listing genomes: %s', error)

    files = [line.strip() for line in output.split('\n') if line.strip()]
    return files


RUN_FILE_DIR = '.local/share/bsp/jobs'


def run_file_stem(job_id: str | None) -> str:
    """The sanitized stem naming a run's .log/.rc/.pid on the cluster.

    Shared so a reattach looks for exactly the files the launch wrote --
    see submit_ssh_job and the API's job_status reattach path.
    """
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(job_id or 'run'))


# A live run's log is never quiet for this long: snakemake logs a status-check
# cycle roughly every 30s for as long as it has jobs in flight.
RUN_STALE_AFTER = 900.0

# How long to wait for the launching shell to echo the detached run's pid.
# Generous, because it is only a guard against hanging: the pid arrives in
# milliseconds when it arrives at all, and the run is already going either way.
_LAUNCH_ACK_TIMEOUT = 30.0

# ---------------------------------------------------------------------------
# Hosting the workflow driver in SLURM rather than on the login node.
#
# setsid + nohup already stop the driver dying with the SSH session, so closing
# a laptop is safe. What they cannot survive is the login node itself: a reboot,
# a maintenance window, or an administrator reaping long-running login-node
# processes. The driver is not a small process -- it is Snakemake, running for
# the whole annotation -- and login nodes are explicitly not where that belongs.
# When it dies, SLURM jobs already queued still finish, but nothing further is
# ever submitted and the run stalls half-done.
#
# Submitting the driver as a job of its own puts it on a compute node under the
# scheduler's protection. Verified on Negishi: a compute node has /usr/bin/sbatch
# and a job submitted from one is accepted, so Snakemake's own SLURM executor
# keeps working from in there.
#
# Deliberately unchanged: the log path and the .rc sentinel. Everything
# downstream -- tail -F, probe_run, is_replayable, the reattach -- keys off
# those two files and nothing else, so none of it needs to know where the
# driver is hosted.
DRIVER_PARTITION = 'cpu'
DRIVER_ACCOUNT = None          # None -> let SLURM pick the default account
# No cap to respect: the cpu partition is MaxTime=UNLIMITED and the normal QOS
# sets no MaxWall (only standby does, at 4h). SLURM charges what a job actually
# uses, not what it asked for, so a generous limit costs nothing but a slightly
# harder backfill -- and with one CPU that is easy to place. Runs measured so
# far take one to three hours; this is a wide margin, not a prediction.
DRIVER_TIME = '3-00:00:00'
DRIVER_CPUS = 2                # Snakemake plus the shell that waits on it
DRIVER_MEM_MB = 8000


def probe_run(job_id: str, connection: SSHConnection) -> dict:
    """One SSH round-trip answering: is this detached run's log worth tailing?

    Returns {"has_log": bool, "exit_code": str|None, "log_idle": float}, where
    log_idle is seconds since the log was last written.

    Liveness is read off the LOG, not the process table. The obvious probe --
    pgrep for the job_id, which is in the driver's command line -- cannot work:
    the job_id is also in this probe's own command line (the log path it stats),
    so the probe matches itself and every run, however long dead, reads as
    alive. The `ps | grep [f]oo` bracket trick does not save it either, because
    the unbracketed copy in the log path is right there in the same command.
    Log mtime answers the question being asked anyway -- "will `tail -F` ever
    produce anything?" -- rather than a proxy for it.

    Anything unparseable reads as "nothing there", the safe answer: callers use
    this to decide whether to start tailing, and `tail -F` on a file that will
    never grow never returns.
    """
    base = f'$HOME/{RUN_FILE_DIR}/{run_file_stem(job_id)}'
    probe = (
        f'now=$(date +%s); '
        f'mt=$(stat -c %Y {base}.log 2>/dev/null || echo 0); '
        f'sz=$(stat -c %s {base}.log 2>/dev/null || echo 0); '
        f'rc=$(cat {base}.rc 2>/dev/null); '
        # A driver hosted in SLURM has a fourth state the log cannot show: it
        # can be sitting in the queue, not yet started, with no log at all. That
        # is emphatically not a dead run, and without this it read as one -- the
        # reattach would refuse a job that was simply waiting for a node.
        f'dj=$(cat {base}.jobid 2>/dev/null); '
        f'ds=$(squeue -h -j "${{dj:-0}}" -o %T 2>/dev/null | head -1); '
        f'echo "$sz|${{rc:--}}|$((now - mt))|${{ds:--}}"'
    )
    ssh = connection.connect()
    stdin, stdout, stderr = ssh.exec_command(probe)
    lines = [ln for ln in stdout.read().decode().strip().splitlines() if ln.strip()]
    parts = lines[-1].split('|') if lines else []
    # Tolerate the 3-field answer too: a run launched before driver jobs
    # existed, or by an older dane-api, still has to be probeable.
    if len(parts) not in (3, 4):
        LOGGER.warning('Unreadable run probe for %s: %r', job_id, lines)
        return {"has_log": False, "exit_code": None, "log_idle": float('inf'),
                "driver_state": None}
    size, rc, idle = (p.strip() for p in parts[:3])
    driver_state = parts[3].strip() if len(parts) == 4 else '-'
    try:
        idle_s = float(idle)
    except ValueError:
        idle_s = float('inf')
    return {
        "has_log": size.isdigit() and int(size) > 0,
        "exit_code": None if rc in ('-', '') else rc,
        "log_idle": idle_s,
        "driver_state": None if driver_state in ('-', '') else driver_state,
    }


# squeue states meaning the driver job will still produce output. COMPLETING is
# excluded on purpose: its log is already written and the log checks below judge
# it correctly, whereas treating it as live would start a tail on a file that is
# about to stop growing.
_DRIVER_PENDING_STATES = ('PENDING', 'CONFIGURING', 'RUNNING', 'RESIZING',
                          'REQUEUED', 'SUSPENDED')


def is_replayable(probe: dict) -> bool:
    """Whether a probe_run() result means `tail -F` will terminate or deliver.

    True in three cases: the driver job is still queued or running (so a log is
    coming, even if there is none yet -- `tail -F` waits for the file, which is
    exactly what it is for), the run finished (sentinel present, so the replay
    ends at it and recovers the real exit code), or its log is still being
    written (so the replay catches up and then follows it live).

    The first case only exists because the driver moved into SLURM. A queued
    driver has no log and no sentinel, which under the old two-case rule read
    identically to a run that had died -- so reopening the page during the queue
    wait refused to reattach, and the job sat there looking abandoned until it
    happened to be polled again after the node was allocated.
    """
    if probe.get("driver_state") in _DRIVER_PENDING_STATES:
        return True
    if not probe.get("has_log"):
        return False
    return probe.get("exit_code") is not None or probe.get("log_idle", float('inf')) < RUN_STALE_AFTER


def driver_job_id(job_id: str, connection: SSHConnection) -> str | None:
    """The SLURM id of this run's driver job, or None if it has no driver job.

    None is the correct answer for a run launched before the driver moved into
    SLURM, and for one launched with in_slurm=False. Callers must treat it as
    "nothing extra to cancel", not as an error.
    """
    base = f'$HOME/{RUN_FILE_DIR}/{run_file_stem(job_id)}'
    ssh = connection.connect()
    _in, out, _err = ssh.exec_command(f'cat {base}.jobid 2>/dev/null')
    lines = [ln.strip() for ln in out.read().decode().splitlines() if ln.strip()]
    value = lines[-1] if lines else ''
    return value if value.isdigit() else None


def build_driver_launch(cmd: str, base: str, safe: str, log: str, rcf: str,
                        jobidf: str, driversh: str,
                        partition: str = DRIVER_PARTITION,
                        account: str | None = DRIVER_ACCOUNT,
                        time_limit: str = DRIVER_TIME,
                        cpus: int = DRIVER_CPUS,
                        mem_mb: int = DRIVER_MEM_MB) -> str:
    """Shell that writes the driver's two scripts and submits the batch job.

    The workflow command goes into a file of its OWN rather than inside a
    `bash -c '...'`, because it legitimately contains single quotes
    (MARGIE_LICENSE_ACCEPTED='2026-07-31' ...) and would not survive being
    wrapped in more of them. Both heredocs are quoted, so nothing in the
    command is expanded by the shell that writes it -- it is stored verbatim
    and interpreted only when the batch job runs it.

    nohup is kept even though a batch job has no controlling terminal to be
    hung up on. It costs nothing and it means the driver script is equally safe
    if it is ever run outside SLURM again.

    The workflow is NOT backgrounded here. On the login node it had to be, so
    the SSH call could return; inside a batch job the opposite is true -- if the
    script exits, SLURM tears the allocation down and takes the run with it. So
    the script waits, and the job lives exactly as long as the run.
    """
    sbatch_directives = [
        f'#SBATCH --job-name=margie-{safe}',
        f'#SBATCH --partition={partition}',
        f'#SBATCH --time={time_limit}',
        '#SBATCH --nodes=1',
        '#SBATCH --ntasks=1',
        f'#SBATCH --cpus-per-task={cpus}',
        f'#SBATCH --mem={mem_mb}',
        # The driver's own stdout is noise; the workflow's real output is
        # redirected to $log below, which is what the GUI tails.
        '#SBATCH --output=/dev/null',
        '#SBATCH --error=/dev/null',
    ]
    if account:
        sbatch_directives.insert(1, f'#SBATCH --account={account}')

    return (
        f'mkdir -p $HOME/.local/share/bsp/jobs && '
        f'rm -f {rcf} {jobidf} && '
        f"cat > {driversh} <<'MARGIE_WORKFLOW_EOF'\n"
        f'#!/bin/bash\n'
        f'export PATH=$HOME/.local/bin:$PATH\n'
        f'{cmd}\n'
        f'MARGIE_WORKFLOW_EOF\n'
        f'chmod +x {driversh}\n'
        f"cat > {base}.sbatch <<'MARGIE_SBATCH_EOF'\n"
        f'#!/bin/bash\n'
        + '\n'.join(sbatch_directives) + '\n'
        f'nohup bash {driversh} > {log} 2>&1\n'
        f'echo $? > {rcf}\n'
        f'MARGIE_SBATCH_EOF\n'
        f'sbatch --parsable {base}.sbatch | tee {jobidf}\n'
    )


def submit_ssh_job(
    cmd,
    connection: SSHConnection,
    job_id: str | None = None,
    poll: float = 1.0,
    reattach: bool = False,
    in_slurm: bool = True,
):
    '''Run a workflow command on the login node, DETACHED, and stream its log.

    reattach=True skips the launch entirely and only streams: the run is
    already going, started by an earlier (now dead) dane-api, and its log and
    exit sentinel are still on disk under the same job_id. cmd is ignored in
    that mode. The tail has always started at line 1 precisely so this would
    replay everything the lost session saw -- until now nothing called it, so a
    dane-api restart meant a live run's job page went permanently blank: no
    logs, no SLURM jobs, phase frozen wherever it was when the API died.

    Yields each output line as it arrives, then a final __EXIT_CODE__: line.
    The contract is unchanged; how the process is hosted is not.

    Previously this did exec_command(..., get_pty=True) and read the channel
    directly, which tied the run's lifetime to the SSH session. Because
    margie.sh traps EXIT and kills the remote dane-api, quitting the GUI closed
    that session, tore down the PTY, and SIGHUP'd dane_wf -- the Snakemake
    driver. Already-submitted SLURM jobs kept running, but nothing further was
    ever submitted, so a run silently stalled half-finished whenever the user
    closed their laptop. A genome annotation takes hours; that is not a
    reasonable thing to require.

    Now the command is started with setsid + nohup, writing to a log file, and
    its exit status to a sentinel. Streaming is a SEPARATE concern: we tail the
    log. Losing the tail (closed GUI, dropped VPN) loses only the live output --
    the run continues, and reattaching later replays the log from the top.
    '''
    # pooled=False: this client is held for the ENTIRE run -- hours -- while the
    # pooled one is shared with every status poll and file listing the GUI makes
    # in the meantime, and is evicted on a 600s TTL. A run must not depend on
    # either. See SSHConnection.connect's docstring.
    ssh = connection.connect(pooled=False)

    safe = run_file_stem(job_id)
    base = f'$HOME/.local/share/bsp/jobs/{safe}'
    log, rcf, pidf = f'{base}.log', f'{base}.rc', f'{base}.pid'
    jobidf, driversh = f'{base}.jobid', f'{base}.driver.sh'

    try:
        if reattach:
            LOGGER.info('Reattaching to run %s, replaying log %s', safe, log)
        elif in_slurm:
            # The driver goes to a compute node under the scheduler, not onto
            # the login node -- see the DRIVER_* constants for why.
            launch = build_driver_launch(cmd, base, safe, log, rcf, jobidf, driversh)
            _in, _out, _err = ssh.exec_command(launch)
            driver_job = ''
            _out.channel.settimeout(_LAUNCH_ACK_TIMEOUT)
            try:
                driver_job = (_out.readline() or '').strip()
            except Exception as exc:
                LOGGER.warning('No SLURM id from the driver submission for %s (%s)',
                               safe, exc)
            if not driver_job.isdigit():
                # sbatch refused. Unlike a login-node launch there is nothing
                # running yet, so this IS a failed start and must be reported
                # as one rather than yielding __LAUNCHED__ and waiting for a
                # log that will never appear.
                err = ''
                try:
                    _err.channel.settimeout(_LAUNCH_ACK_TIMEOUT)
                    err = (_err.read().decode() or '').strip()
                except Exception:
                    pass
                raise RuntimeError(
                    f'Could not submit the workflow driver to SLURM: '
                    f'{err or driver_job or "no job id returned"}')
            LOGGER.info('Driver submitted as SLURM job %s, log %s', driver_job, log)
        else:
            # setsid detaches from the session so no SIGHUP reaches it; nohup
            # covers the gap before setsid takes effect. The exit code is
            # written by the same shell that runs the command, so it is
            # recorded even though nobody is attached.
            #
            # The pid is recorded by the launching shell itself. It used to be a
            # SECOND exec_command issued from here, which bought nothing and
            # added a failure point squarely between "the run has started" and
            # "we are watching it": when that call was the one that failed, the
            # workflow was already running on the cluster but the caller saw
            # only an exception.
            launch = (
                f'mkdir -p $HOME/.local/share/bsp/jobs && '
                f'rm -f {rcf} && '
                f'export PATH=$HOME/.local/bin:$PATH && '
                f"nohup setsid bash -c '{{ {cmd} ; }} > {log} 2>&1; echo $? > {rcf}' "
                f'>/dev/null 2>&1 & echo $! | tee {pidf}'
            )
            _in, _out, _err = ssh.exec_command(launch)
            # ONE line, with a deadline -- never .read(), and never
            # recv_exit_status().
            #
            # Both of those wait for the channel to reach EOF, and EOF here does
            # not mean "the pid has been printed". `A && B && nohup setsid ... &`
            # backgrounds the whole and-list, so a shell sits there waiting for
            # the workflow with this channel still open on its fd 1 -- confirmed
            # on a live run: the launcher was 13 minutes old, /proc/<pid>/fd/1
            # still pointing at the channel pipe. EOF therefore arrives when the
            # RUN ends, hours later.
            #
            # So the generator blocked here forever: it never yielded, never
            # started the tail, and never parsed a line. The job page showed the
            # phase frozen at "Submitting via SSH", an empty Slurm Jobs table and
            # an empty log for the entire run, while the run itself went on
            # perfectly well. Every run since detached launches were introduced
            # was affected; the last run with provenance predates them.
            #
            # A missed pid is survivable -- the launch shell tees it to $pidf
            # anyway -- so a timeout here logs and carries on rather than
            # failing a run that has already started.
            pid = ''
            _out.channel.settimeout(_LAUNCH_ACK_TIMEOUT)
            try:
                pid = (_out.readline() or '').strip()
            except Exception as exc:
                LOGGER.warning('No pid from the launch of %s within %ss (%s); '
                               'the run is started regardless',
                               safe, _LAUNCH_ACK_TIMEOUT, exc)
            LOGGER.info('Detached run started (pid %s), log %s', pid, log)

        # From here on the workflow IS running on the cluster. Everything after
        # this point only decides how well we can watch it, so the caller is
        # told now -- see job_runner.run_ssh_task, which uses this to tell "the
        # run never started" (a real failure) apart from "we lost sight of a run
        # that is still going" (not a failure at all).
        yield '__LAUNCHED__'

        # Tail from the beginning so a reattach replays everything already
        # written. -F rather than -f: the log may not exist for a moment after
        # launch.
        tail_cmd = f'tail -n +1 -F {log} 2>/dev/null'
        t_in, t_out, t_err = ssh.exec_command(tail_cmd)
        chan = t_out.channel
        chan.settimeout(poll)
    except BaseException:
        try:
            ssh.close()
        except Exception:
            pass
        raise

    # Completion is detected by stat-ing the sentinel over ONE long-lived SFTP
    # session. The previous version called exec_command() once per poll to `cat`
    # it -- a new SSH session every second, never closed. sshd's MaxSessions is
    # 10 by default, so the transport refused new sessions within ~10s, the loop
    # died, and the generator returned __DETACHED__. job_runner read that as
    # "stopped watching" and broke out, so streaming ended seconds into every
    # run: no logs stored, phase frozen at its first value, and no SLURM jobs
    # ever recorded. Runs before this change stored ~2MB of log; after, zero.
    sftp = None
    try:
        sftp = ssh.open_sftp()
    except Exception as exc:
        LOGGER.warning('Could not open SFTP for completion checks: %s', exc)

    def _finished():
        """Exit code if the run has finished, else None."""
        if sftp is None:
            return None
        try:
            with sftp.open(rcf, 'r') as fh:
                txt = fh.read().decode('utf-8', 'replace').strip()
            return txt or None
        except IOError:
            return None                       # not there yet
        except Exception:
            return None

    exit_code = None
    buf = ''
    last_check = 0.0
    CHECK_EVERY = 5.0                         # the sentinel is not urgent
    try:
        while True:
            try:
                data = chan.recv(65536)
                if not data:
                    raise EOFError
                buf += data.decode('utf-8', 'replace')
                while '\n' in buf:
                    line, buf = buf.split('\n', 1)
                    LOGGER.debug('[remote] %s', line.rstrip())
                    yield line.rstrip()
                continue                       # more may be waiting; drain first
            except socket.timeout:
                pass
            except EOFError:
                break

            now = time.time()
            if now - last_check < CHECK_EVERY:
                continue
            last_check = now
            got = _finished()
            if got:
                # Drain whatever the tail has not delivered yet.
                try:
                    chan.settimeout(1.0)
                    while True:
                        data = chan.recv(65536)
                        if not data:
                            break
                        buf += data.decode('utf-8', 'replace')
                except Exception:
                    pass
                for line in buf.split('\n'):
                    if line.strip():
                        yield line.rstrip()
                try:
                    exit_code = int(got.split()[0])
                except (ValueError, IndexError):
                    exit_code = 1
                break
    finally:
        # ssh is closed here too, unlike everywhere else in this module: this
        # client is not the pool's (connect(pooled=False) above), so nothing
        # else can be using it and leaving it open would leak one connection
        # per run.
        for closer in (chan, sftp, ssh):
            try:
                if closer is not None:
                    closer.close()
            except Exception:
                pass

    if exit_code is None:
        # The tail ended without a sentinel: the run is still going, we just
        # stopped watching. Do NOT report an exit code -- that would mark a
        # live run as finished.
        LOGGER.info('Detached from run %s; it continues in the background', safe)
        yield '__DETACHED__'
    else:
        LOGGER.info('Remote execution completed with exit code: %d', exit_code)
        yield f'__EXIT_CODE__:{exit_code}'


def submit_slurm_job(
    script_content,
    connection: SSHConnection,
    nodes=1,
    cpus=4,
    mem='4G',
    time='00:30:00',
):
    """Write a SLURM batch script and submit it via sbatch."""
    ssh = connection.connect()

    stdin, stdout, stderr = ssh.exec_command('touch im-here.flag')

    # Create SLURM script
    slurm_script = f"""#!/bin/bash
#SBATCH -A lindems
#SBATCH --partition=cpu
#SBATCH --nodes={nodes}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --job-name=remote_job

# TODO: Here, we need script_content

source /etc/profile

{script_content}
    """
    # Write script and submit
    stdin, stdout, stderr = ssh.exec_command(
        f'cat > ~/job.sh << "EOF"\n{slurm_script}\nEOF\n'
        f'sbatch ~/job.sh'
    )

    job_id = stdout.read().decode().strip()
    try:
        stderr_content = stderr.read().decode().strip()
    except OSError:
        stderr_content = 'None'

    LOGGER.info('submit_slurm_job stdout: %s, stderr: %s', job_id, stderr_content)
    pass  # pooled client: closing it would defeat SSHConnection's pool

    # Extract just the job number (sbatch returns "Submitted batch job 12345")
    if "Submitted batch job" in job_id:
        job_id = job_id.split()[-1]
    return job_id


def check_slurm_job_status(
    job_id,
    connection: SSHConnection,
):
    """Check the status of a single SLURM job via squeue then sacct.

    Returns: dict with status info (state, elapsed_time, etc.)
    """
    ssh = connection.connect()

    # Use squeue to check if job is running/pending
    stdin, stdout, stderr = ssh.exec_command(f'squeue -j {job_id} --format="%T %M %j %a %l" --noheader')
    squeue_output = stdout.read().decode().strip()

    if squeue_output:
        parts = squeue_output.split()
        state = parts[0] if len(parts) > 0 else "UNKNOWN"
        elapsed = parts[1] if len(parts) > 1 else "0:00"
        job_name = parts[2] if len(parts) > 2 else "0:00"
        account = parts[3] if len(parts) > 3 else "0:00"
        limit = parts[4] if len(parts) > 4 else "0:00"
        pass  # pooled client: closing it would defeat SSHConnection's pool
        return {"state": state, "elapsed_time": elapsed, "job_name": job_name, "account": account, "time limit": limit, "exists": True}

    # Job not in queue, check sacct for completed/failed jobs
    stdin, stdout, stderr = ssh.exec_command(f'sacct -j {job_id} --format=JobName,State,Elapsed --noheader | head -1')
    sacct_output = stdout.read().decode().strip()

    pass  # pooled client: closing it would defeat SSHConnection's pool

    if sacct_output:
        parts = sacct_output.split()
        job_name = parts[0] if len(parts) > 0 else "UNKNOWN"
        state = parts[1] if len(parts) > 1 else "UNKNOWN"
        elapsed = parts[2] if len(parts) > 2 else "0:00"
        return {"job_name": job_name, "state": state, "elapsed_time": elapsed, "exists": True}

    return {"state": "NOT_FOUND", "elapsed_time": "0:00", "exists": False}


def check_multiple_slurm_jobs(
    job_ids: list[str],
    connection: SSHConnection,
) -> dict[str, dict]:
    """Check status of multiple SLURM jobs in a single SSH call.

    Returns a dict mapping each job_id to {"state": ..., "time": ...}.
    """
    if not job_ids:
        return {}

    ssh = connection.connect()

    results = {}
    ids_str = ",".join(job_ids)

    # Try squeue first for active jobs
    stdin, stdout, stderr = ssh.exec_command(
        f'squeue -j {ids_str} --format="%i %T %M" --noheader 2>/dev/null'
    )
    squeue_output = stdout.read().decode().strip()

    found_ids = set()
    if squeue_output:
        for line in squeue_output.splitlines():
            parts = line.split()
            if len(parts) >= 3:
                jid, state, elapsed = parts[0], parts[1], parts[2]
                results[jid] = {"state": state, "time": elapsed}
                found_ids.add(jid)

    # For any IDs not found in squeue, check sacct
    missing = [jid for jid in job_ids if jid not in found_ids]
    if missing:
        missing_str = ",".join(missing)
        stdin, stdout, stderr = ssh.exec_command(
            f'sacct -j {missing_str} --format=JobID,State,Elapsed --noheader --parsable2 2>/dev/null'
        )
        sacct_output = stdout.read().decode().strip()
        if sacct_output:
            for line in sacct_output.splitlines():
                parts = line.split("|")
                if len(parts) >= 3:
                    jid = parts[0].split(".")[0]  # strip .batch/.extern suffix
                    if jid in missing and jid not in results:
                        results[jid] = {"state": parts[1], "time": parts[2]}

    pass  # pooled client: closing it would defeat SSHConnection's pool
    return results


def get_job_genome(log_path: str, connection: SSHConnection) -> str:
    """Reads a SLURM job's own log file for its "wildcards: genome=<name>"
    line. Fallback only: job_runner.run_ssh_task now reads this same line
    directly from the live orchestrator stream (Snakemake's --verbose
    output, see WILDCARDS_GENOME_RE), since this remote file gets cleaned
    up shortly after the job finishes and isn't always still around by the
    time _slurm_status_checker's polling loop gets to it. Returns "" if the
    file doesn't exist (already cleaned up, or job hasn't started) or has
    no genome wildcard (e.g. quast_batch/gtdbtk_batch, which process every
    genome at once).
    """
    ssh = connection.connect()
    stdin, stdout, stderr = ssh.exec_command(
        f"grep -m1 'wildcards:' {shlex.quote(log_path)} 2>/dev/null"
    )
    line = stdout.read().decode().strip()
    pass  # pooled client: closing it would defeat SSHConnection's pool
    match = re.search(r'\bgenome=([^\s,]+)', line)
    return match.group(1) if match else ""


def find_active_jobs_in_workdir(
    work_dir: str,
    username: str,
    connection: SSHConnection,
) -> list[dict]:
    """Lists this user's SLURM jobs (RUNNING or PENDING) whose working
    directory matches work_dir exactly (trailing-slash-normalized on both
    sides).

    Used to check whether a job that looks "running"/"pending" in
    persisted history (because dane-api restarted and lost live track of
    it) is actually still active on the cluster -- work_dir is set once
    at job creation and never changes, and squeue's WorkDir column
    reflects the directory a job's driver process was launched from.

    Returns a list of {"job_id": ..., "state": ..., "time": ...} dicts --
    usually 0 or 1 entries; empty means nothing currently active matches
    this work_dir.
    """
    ssh = connection.connect()
    stdin, stdout, stderr = ssh.exec_command(
        f'squeue -u {username} --format="%i|%T|%Z|%M" --noheader'
    )
    output = stdout.read().decode().strip()
    pass  # pooled client: closing it would defeat SSHConnection's pool

    target = work_dir.rstrip('/')
    matches = []
    for line in output.splitlines():
        parts = line.split('|')
        if len(parts) != 4:
            continue
        slurm_job_id, state, workdir, elapsed = parts
        if workdir.rstrip('/') == target:
            matches.append({"job_id": slurm_job_id, "state": state, "time": elapsed})
    return matches


_RULE_FROM_LOG_RE = re.compile(r'/slurm_logs/(?:rule_(\w+)|group_([^_]+)_)')


def enrich_slurm_jobs_from_logs(
    work_dir: str,
    matches: list[dict],
    connection: SSHConnection,
) -> list[dict]:
    """Adds 'rule' and 'genome' to each match dict from find_active_jobs_in_workdir
    by scanning snakemake's slurm_logs directory. Single SSH call regardless of
    how many jobs. No-op if matches is empty."""
    if not matches:
        return matches

    logs_dir = shlex.quote(f"{work_dir}/.snakemake/slurm_logs")
    cmd = (
        f"find {logs_dir} -name '*.log' 2>/dev/null | "
        "while read f; do "
        "  id=$(basename \"$f\" .log); "
        "  genome=$(grep -m1 'wildcards: genome=' \"$f\" 2>/dev/null "
        "           | sed 's/.*genome=//;s/[[:space:]].*//'); "
        "  echo \"$id|$f|$genome\"; "
        "done"
    )
    ssh = connection.connect()
    try:
        _, stdout, _ = ssh.exec_command(cmd)
        output = stdout.read().decode().strip()
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool

    log_info: dict[str, tuple[str, str]] = {}
    for line in output.splitlines():
        parts = line.split("|", 2)
        if len(parts) == 3:
            jid, path, genome = parts
            log_info[jid] = (path, genome.strip())

    enriched = []
    for m in matches:
        m = dict(m)
        jid = m["job_id"]
        if jid in log_info:
            path, genome = log_info[jid]
            rule_match = _RULE_FROM_LOG_RE.search(path)
            if rule_match:
                m["rule"] = rule_match.group(1) or rule_match.group(2) or None
            if genome:
                m["genome"] = genome
        enriched.append(m)
    return enriched


def read_latest_snakemake_log(
    work_dir: str,
    connection: SSHConnection,
    tail_lines: int = 300,
) -> str:
    """Returns the tail of the most recent Snakemake master log from
    {work_dir}/.snakemake/log/. Used as a fallback when the API restarted
    and the in-memory log buffer was lost. Returns "" if no log is found."""
    log_glob = shlex.quote(f"{work_dir}/.snakemake/log")
    cmd = (
        f"latest=$(ls -t {log_glob}/*.log 2>/dev/null | head -1); "
        f"[ -n \"$latest\" ] && tail -{tail_lines} \"$latest\" 2>/dev/null || true"
    )
    ssh = connection.connect()
    try:
        _, stdout, _ = ssh.exec_command(cmd)
        return stdout.read().decode()
    except Exception:
        return ""
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool


def cancel_slurm_jobs(
    job_ids: list[str],
    connection: SSHConnection,
) -> None:
    """Cancel multiple SLURM jobs via scancel.

    Args:
        job_ids: List of SLURM job IDs to cancel
        connection: SSH connection to the cluster
    """
    if not job_ids:
        LOGGER.info('No SLURM jobs to cancel')
        return

    ssh = connection.connect()
    ids_str = ",".join(job_ids)

    LOGGER.info('Cancelling SLURM jobs: %s', ids_str)
    stdin, stdout, stderr = ssh.exec_command(f'scancel {ids_str}')

    # Read output to ensure command completes
    stdout.read()
    error = stderr.read().decode().strip()

    if error:
        LOGGER.warning('scancel stderr: %s', error)
    else:
        LOGGER.info('Successfully cancelled %d SLURM job(s)', len(job_ids))

    pass  # pooled client: closing it would defeat SSHConnection's pool


def kill_remote_process(
    process_pattern: str,
    connection: SSHConnection,
) -> None:
    """Kill remote processes matching a pattern.

    Args:
        process_pattern: Pattern to match in process command line (for pkill -f)
        connection: SSH connection to the cluster
    """
    ssh = connection.connect()

    # Use pkill -f to kill processes matching the pattern
    # The -f flag matches against the full command line
    LOGGER.info('Killing remote processes matching: %s', process_pattern)
    stdin, stdout, stderr = ssh.exec_command(f'pkill -f "{process_pattern}"')

    # Read output to ensure command completes
    stdout.read()
    error = stderr.read().decode().strip()

    # pkill returns 0 if at least one process was killed, 1 if none matched
    # So we don't treat non-zero exit as an error
    if error:
        LOGGER.debug('pkill stderr: %s', error)

    LOGGER.info('Sent kill signal to processes matching: %s', process_pattern)
    pass  # pooled client: closing it would defeat SSHConnection's pool
