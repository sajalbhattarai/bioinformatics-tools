'''
SSH-based SLURM operations for interfacing with HPC clusters.

Uses paramiko to run login node and SLURM commands. Since we use Snakemake,
which controls SLURM batching and queueing, we can mainly run on login node.

All functions are API-layer only. Pass a per-user SSHConnection built
with make_user_connection() for every call.
'''
import logging
import socket
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


def submit_ssh_job(
    cmd,
    connection: SSHConnection,
    job_id: str | None = None,
    poll: float = 1.0,
):
    '''Run a workflow command on the login node, DETACHED, and stream its log.

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
    ssh = connection.connect()

    tag = job_id or 'run'
    safe = re.sub(r'[^A-Za-z0-9_.-]', '_', str(tag))
    base = f'$HOME/.local/share/bsp/jobs/{safe}'
    log, rcf, pidf = f'{base}.log', f'{base}.rc', f'{base}.pid'

    # setsid detaches from the session so no SIGHUP reaches it; nohup covers the
    # gap before setsid takes effect. The exit code is written by the same shell
    # that runs the command, so it is recorded even though nobody is attached.
    launch = (
        f'mkdir -p $HOME/.local/share/bsp/jobs && '
        f'rm -f {rcf} && '
        f'export PATH=$HOME/.local/bin:$PATH && '
        f"nohup setsid bash -c '{{ {cmd} ; }} > {log} 2>&1; echo $? > {rcf}' "
        f'>/dev/null 2>&1 & echo $!'
    )
    _in, _out, _err = ssh.exec_command(launch)
    pid = (_out.read().decode() or '').strip().splitlines()
    pid = pid[-1] if pid else ''
    _out.channel.recv_exit_status()
    ssh.exec_command(f'echo {pid} > {pidf}')
    LOGGER.info('Detached run started (pid %s), log %s', pid, log)

    # Tail from the beginning so a reattach replays everything already written.
    # -F rather than -f: the log may not exist for a moment after launch.
    tail_cmd = f'tail -n +1 -F {log} 2>/dev/null'
    t_in, t_out, t_err = ssh.exec_command(tail_cmd)
    chan = t_out.channel
    chan.settimeout(poll)

    exit_code = None
    buf = ''
    try:
        while True:
            try:
                data = chan.recv(65536)
                if not data:
                    raise EOFError
                buf += data.decode('utf-8', 'replace')
                while '\n' in buf:
                    line, buf = buf.split('\n', 1)
                    LOGGER.info('[remote] %s', line.rstrip())
                    yield line.rstrip()
            except socket.timeout:
                pass
            except EOFError:
                break

            # Finished? The sentinel is the only reliable signal -- tail -F
            # never ends on its own.
            c_in, c_out, c_err = ssh.exec_command(f'cat {rcf} 2>/dev/null')
            got = (c_out.read().decode() or '').strip()
            if got:
                # Drain whatever the tail has not delivered yet.
                try:
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
        try:
            chan.close()
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

    pass  # pooled client: closing it would defeat SSHConnection's pool


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
