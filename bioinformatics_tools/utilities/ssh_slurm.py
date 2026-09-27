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
    """Returns the sanitised stem naming a run's .log/.rc/.pid on the cluster.

    Shared so a reattach finds the files the launch wrote.
    """
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(job_id or 'run'))


# A live run's log is never quiet this long: snakemake logs a status check about every 30s.
RUN_STALE_AFTER = 900.0

# How long to wait for the launching shell to echo the detached run's pid (a hang guard).
_LAUNCH_ACK_TIMEOUT = 30.0

# ---- Workflow driver hosted in SLURM ----
# The Snakemake driver runs as its own batch job on a compute node, so login-node
# reboots or reaping cannot stall a run; compute nodes can call sbatch. The log
# path and .rc sentinel are unchanged, so tailing and reattach work as before.
DRIVER_PARTITION = 'cpu'
DRIVER_ACCOUNT = None          # None -> let SLURM pick the default account
# The cpu partition and normal QOS set no wall-time cap and SLURM charges actual
# use, so this is a wide margin over the one to three hours runs take.
DRIVER_TIME = '7-00:00:00'
DRIVER_CPUS = 2                # Snakemake plus the shell that waits on it
DRIVER_MEM_MB = 8000

# The six --time formats sbatch accepts; a user-edited value is validated here
# so a typo falls back to the default instead of failing the start.
_SLURM_TIME_RE = re.compile(
    r'^(?:\d+'                       # minutes
    r'|\d+:\d{1,2}'                  # minutes:seconds
    r'|\d+:\d{1,2}:\d{1,2}'          # hours:minutes:seconds
    r'|\d+-\d{1,2}'                  # days-hours
    r'|\d+-\d{1,2}:\d{1,2}'          # days-hours:minutes
    r'|\d+-\d{1,2}:\d{1,2}:\d{1,2}'  # days-hours:minutes:seconds
    r')$'
)


def _valid_driver_time(value: str | None) -> str | None:
    """Returns *value* if sbatch accepts it as --time, else None (use the default)."""
    if value is None:
        return None
    candidate = str(value).strip()
    if not candidate:
        return None
    if not _SLURM_TIME_RE.match(candidate):
        LOGGER.warning('Ignoring invalid driver walltime %r (expected D-HH:MM:SS or HH:MM:SS); '
                       'falling back to %s', value, DRIVER_TIME)
        return None
    return candidate


def probe_run(job_id: str, connection: SSHConnection) -> dict:
    """Checks in one SSH round-trip whether a detached run's log is worth tailing.

    Returns {"has_log": bool, "exit_code": str|None, "log_idle": float}, where
    log_idle is seconds since the log was last written. Liveness comes from log
    mtime, since a pgrep for the job_id would match the probe itself. Anything
    unparseable reads as "nothing there", because `tail -F` on a file that never
    grows never returns.
    """
    base = f'$HOME/{RUN_FILE_DIR}/{run_file_stem(job_id)}'
    probe = (
        f'now=$(date +%s); '
        f'mt=$(stat -c %Y {base}.log 2>/dev/null || echo 0); '
        f'sz=$(stat -c %s {base}.log 2>/dev/null || echo 0); '
        f'rc=$(cat {base}.rc 2>/dev/null); '
        # A driver job still queued has no log yet but is not dead.
        f'dj=$(cat {base}.jobid 2>/dev/null); '
        f'ds=$(squeue -h -j "${{dj:-0}}" -o %T 2>/dev/null | head -1); '
        f'echo "$sz|${{rc:--}}|$((now - mt))|${{ds:--}}"'
    )
    ssh = connection.connect()
    stdin, stdout, stderr = ssh.exec_command(probe)
    lines = [ln for ln in stdout.read().decode().strip().splitlines() if ln.strip()]
    parts = lines[-1].split('|') if lines else []
    # Also accepts the 3-field answer from runs launched before driver jobs existed.
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


# squeue states in which the driver will still produce output. COMPLETING is
# excluded: its log is complete and the log checks judge it correctly.
_DRIVER_PENDING_STATES = ('PENDING', 'CONFIGURING', 'RUNNING', 'RESIZING',
                          'REQUEUED', 'SUSPENDED')


def is_replayable(probe: dict) -> bool:
    """Returns whether a probe_run() result means `tail -F` will terminate or deliver.

    True when the driver job is queued or running, the run has finished
    (sentinel present), or its log is still being written.
    """
    if probe.get("driver_state") in _DRIVER_PENDING_STATES:
        return True
    if not probe.get("has_log"):
        return False
    return probe.get("exit_code") is not None or probe.get("log_idle", float('inf')) < RUN_STALE_AFTER


def driver_job_id(job_id: str, connection: SSHConnection) -> str | None:
    """Returns the SLURM id of this run's driver job, or None if it has none
    (older runs and in_slurm=False); None means nothing extra to cancel.
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
    """Returns the shell that writes the driver's two scripts and submits the batch job.

    The workflow command goes into its own file via a quoted heredoc, since it
    contains single quotes. The batch script waits on the workflow rather than
    backgrounding it, so the job lives exactly as long as the run.
    """
    sbatch_directives = [
        f'#SBATCH --job-name=margie-{safe}',
        f'#SBATCH --partition={partition}',
        f'#SBATCH --time={time_limit}',
        '#SBATCH --nodes=1',
        '#SBATCH --ntasks=1',
        f'#SBATCH --cpus-per-task={cpus}',
        f'#SBATCH --mem={mem_mb}',
        # The driver's stdout is noise; the workflow's output goes to $log, which the GUI tails.
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
    driver_account: str | None = None,
    driver_partition: str | None = None,
    driver_time: str | None = None,
):
    '''Runs a workflow command detached and streams its log.

    The command is started with setsid + nohup (or as a SLURM driver job),
    writing to a log file and its exit status to a sentinel, so closing the GUI
    never stops the run. reattach=True skips the launch (cmd is ignored) and
    replays the existing log from line 1.

    Yields each output line as it arrives, then a final __EXIT_CODE__: line.
    '''
    # pooled=False: this client is held for the whole run, independent of the
    # shared pool (see SSHConnection.connect).
    ssh = connection.connect(pooled=False)

    safe = run_file_stem(job_id)
    base = f'$HOME/.local/share/bsp/jobs/{safe}'
    log, rcf, pidf = f'{base}.log', f'{base}.rc', f'{base}.pid'
    jobidf, driversh = f'{base}.jobid', f'{base}.driver.sh'

    try:
        if reattach:
            LOGGER.info('Reattaching to run %s, replaying log %s', safe, log)
        elif in_slurm:
            # The driver goes to a compute node (see the DRIVER_* constants).
            launch = build_driver_launch(
                cmd, base, safe, log, rcf, jobidf, driversh,
                partition=(driver_partition or DRIVER_PARTITION),
                account=(driver_account if driver_account is not None else DRIVER_ACCOUNT),
                time_limit=(_valid_driver_time(driver_time) or DRIVER_TIME),
            )
            _in, _out, _err = ssh.exec_command(launch)
            driver_job = ''
            _out.channel.settimeout(_LAUNCH_ACK_TIMEOUT)
            try:
                driver_job = (_out.readline() or '').strip()
            except Exception as exc:
                LOGGER.warning('No SLURM id from the driver submission for %s (%s)',
                               safe, exc)
            if not driver_job.isdigit():
                # sbatch refused: nothing is running, so this is a failed start.
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
            # covers the gap before setsid takes effect. The same shell writes
            # the exit code and records the pid.
            launch = (
                f'mkdir -p $HOME/.local/share/bsp/jobs && '
                f'rm -f {rcf} && '
                f'export PATH=$HOME/.local/bin:$PATH && '
                f"nohup setsid bash -c '{{ {cmd} ; }} > {log} 2>&1; echo $? > {rcf}' "
                f'>/dev/null 2>&1 & echo $! | tee {pidf}'
            )
            _in, _out, _err = ssh.exec_command(launch)
            # Reads one line with a deadline, never .read() or recv_exit_status():
            # the backgrounded shell keeps the channel open until the run ends.
            # A missed pid is survivable (it is also written to $pidf).
            pid = ''
            _out.channel.settimeout(_LAUNCH_ACK_TIMEOUT)
            try:
                pid = (_out.readline() or '').strip()
            except Exception as exc:
                LOGGER.warning('No pid from the launch of %s within %ss (%s); '
                               'the run is started regardless',
                               safe, _LAUNCH_ACK_TIMEOUT, exc)
            LOGGER.info('Detached run started (pid %s), log %s', pid, log)

        # The workflow is now running; job_runner.run_ssh_task uses this marker to
        # tell a failed start from a lost view of a running job.
        yield '__LAUNCHED__'

        # Tails from the beginning so a reattach replays everything; -F because
        # the log may not exist yet.
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

    # Completion is detected by stat-ing the sentinel over one long-lived SFTP
    # session, which stays within sshd's MaxSessions limit.
    sftp = None
    try:
        sftp = ssh.open_sftp()
    except Exception as exc:
        LOGGER.warning('Could not open SFTP for completion checks: %s', exc)

    # SFTP does not expand `$HOME`, so the sentinel path is given relative to the home directory.
    rcf_sftp = rcf[len('$HOME/'):] if rcf.startswith('$HOME/') else rcf

    def _finished():
        """Returns the exit code if the run has finished, else None."""
        if sftp is None:
            return None
        try:
            with sftp.open(rcf_sftp, 'r') as fh:
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
                # Drains whatever the tail has not delivered yet.
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
        # This client is not pooled (connect(pooled=False)), so it is closed here.
        for closer in (chan, sftp, ssh):
            try:
                if closer is not None:
                    closer.close()
            except Exception:
                pass

    if exit_code is None:
        # The tail ended without a sentinel: the run is still going, so no exit
        # code is reported.
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
    """Reads a SLURM job's log for its "wildcards: genome=<name>" line.

    Fallback for job_runner's live parsing; the file may already be cleaned up.
    Returns "" if the file is missing or has no genome wildcard (batch rules).
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
    """Lists this user's RUNNING or PENDING SLURM jobs whose working directory
    matches work_dir (trailing slashes normalised).

    Used to check whether a job that looks active in history is still on the
    cluster. Returns a list of {"job_id": ..., "state": ..., "time": ...} dicts.
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
    """Adds 'rule' and 'genome' to each match from find_active_jobs_in_workdir by
    scanning snakemake's slurm_logs directory in one SSH call. No-op if matches is empty."""
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
    """Returns the tail of the newest Snakemake log in {work_dir}/.snakemake/log/,
    a fallback when the in-memory log was lost. Returns "" if none is found."""
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
