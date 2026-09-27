"""
Single-machine stand-in for the sbatch, squeue, sacct and scancel commands the API calls.

A job is a detached bash process; its record lives in ~/.margie-dev/slurm/<id>/.
Only the flags and format fields used by utilities/ssh_slurm.py are supported.
"""
from __future__ import annotations

import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get('MARGIE_DEV_SLURM', Path.home() / '.margie-dev' / 'slurm'))
USER = os.environ.get('USER') or os.environ.get('LOGNAME') or 'user'


# ---- job records ----

def _dir(job: str) -> Path:
    return ROOT / str(job)


def _read(job: str, key: str, default: str = '') -> str:
    try:
        return (_dir(job) / key).read_text().strip()
    except OSError:
        return default


def _write(job: str, key: str, value: str) -> None:
    (_dir(job) / key).write_text(f'{value}\n')


def _alive(pid: str) -> bool:
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def _state(job: str) -> str:
    """Returns RUNNING while the process lives, otherwise its final state."""
    end = _read(job, 'state')
    if end in ('COMPLETED', 'FAILED', 'CANCELLED'):
        return end
    if _alive(_read(job, 'pid')):
        return 'RUNNING'
    # The process ended without recording a state, e.g. it was killed.
    rc = _read(job, 'exit')
    state = 'COMPLETED' if rc == '0' else 'FAILED'
    _write(job, 'state', state)
    return state


def _elapsed(job: str) -> str:
    start = float(_read(job, 'start', '0') or 0)
    end = float(_read(job, 'end', '0') or 0) or time.time()
    s = max(0, int(end - start))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f'{h:02d}:{m:02d}:{sec:02d}' if h else f'{m:02d}:{sec:02d}'


def _jobs() -> list[str]:
    if not ROOT.is_dir():
        return []
    return sorted((p.name for p in ROOT.iterdir() if p.name.isdigit()), key=int)


# ---- sbatch ----

def sbatch(argv: list[str]) -> int:
    parsable = '--parsable' in argv
    args = [a for a in argv if a != '--parsable']
    if not args:
        print('sbatch: no script given', file=sys.stderr)
        return 1
    script = Path(os.path.expanduser(args[-1])).resolve()
    if not script.is_file():
        print(f'sbatch: {script}: no such file', file=sys.stderr)
        return 1

    # Reads #SBATCH lines, then command-line options, which take precedence.
    opts: dict[str, str] = {}
    for line in script.read_text(errors='replace').splitlines():
        m = re.match(r'#SBATCH\s+(--?[\w-]+)(?:[= ]\s*(\S.*))?', line)
        if m:
            opts[m.group(1)] = (m.group(2) or '').strip()
    it = iter(args[:-1])
    for a in it:
        if '=' in a:
            k, v = a.split('=', 1)
            opts[k] = v
        elif a.startswith('-'):
            opts[a] = next(it, '')

    ROOT.mkdir(parents=True, exist_ok=True)
    counter = ROOT / '.next'
    job = str(int(counter.read_text()) if counter.exists() else 1000)
    counter.write_text(str(int(job) + 1))
    _dir(job).mkdir()

    workdir = os.path.expanduser(opts.get('--chdir') or opts.get('-D') or str(script.parent))
    name = opts.get('--job-name') or opts.get('-J') or script.name
    out = (opts.get('--output') or opts.get('-o') or f'slurm-{job}.out').replace('%j', job).replace('%x', name)
    out = out if os.path.isabs(os.path.expanduser(out)) else os.path.join(workdir, out)
    err = (opts.get('--error') or opts.get('-e') or '').replace('%j', job).replace('%x', name)

    for key, value in (('name', name), ('workdir', workdir), ('script', str(script)), ('start', str(time.time()))):
        _write(job, key, value)

    # The detached job records its own exit code, end time and state.
    record = _dir(job)
    wrapper = (
        f'cd {sh(workdir)} || exit 1\n'
        f'SLURM_JOB_ID={job} SLURM_JOB_NAME={sh(name)} SLURM_SUBMIT_DIR={sh(workdir)} bash {sh(str(script))}\n'
        f'rc=$?\n'
        f'echo $rc > {sh(str(record / "exit"))}\n'
        f'date +%s > {sh(str(record / "end"))}\n'
        f'[ -f {sh(str(record / "state"))} ] || {{ [ $rc -eq 0 ] && echo COMPLETED || echo FAILED; }} > {sh(str(record / "state"))}\n'
    )
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    stdout = open(out, 'ab')
    stderr = open(os.path.join(workdir, err), 'ab') if err else subprocess.STDOUT
    proc = subprocess.Popen(['bash', '-c', wrapper], stdout=stdout, stderr=stderr, stdin=subprocess.DEVNULL, start_new_session=True)
    _write(job, 'pid', str(proc.pid))
    print(job if parsable else f'Submitted batch job {job}')
    return 0


def sh(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


# ---- squeue / sacct ----

FIELDS = {
    '%i': lambda j: j,
    '%T': _state,
    '%M': _elapsed,
    '%j': lambda j: _read(j, 'name'),
    '%a': lambda j: 'dev',
    '%l': lambda j: 'UNLIMITED',
    '%Z': lambda j: _read(j, 'workdir'),
    '%u': lambda j: USER,
}


def _pick(argv: list[str]) -> tuple[list[str], dict[str, str]]:
    ids: list[str] = []
    opts: dict[str, str] = {}
    it = iter(argv)
    for a in it:
        if a in ('-h', '--noheader', '--parsable2'):
            opts[a] = '1'
        elif a in ('-j', '--jobs', '-u', '--user', '-o', '--format'):
            opts[a] = next(it, '')
        elif a.startswith(('--format=', '--jobs=', '--user=')):
            k, v = a.split('=', 1)
            opts[k] = v
        else:
            print(f'dev-slurm: unsupported argument {a!r}', file=sys.stderr)
            sys.exit(2)
    for key in ('-j', '--jobs'):
        if opts.get(key):
            ids = [x for x in re.split(r'[,\s]+', opts[key]) if x]
    return ids, opts


def squeue(argv: list[str]) -> int:
    ids, opts = _pick(argv)
    fmt = (opts.get('-o') or opts.get('--format') or '%i %j %T %M').strip('"')
    jobs = [j for j in (ids or _jobs()) if _dir(j).is_dir() and _state(j) == 'RUNNING']
    if not ('-h' in opts or '--noheader' in opts):
        print(re.sub(r'%\w', lambda m: m.group(0)[1:].upper(), fmt))
    for j in jobs:
        print(re.sub(r'%\w', lambda m: str(FIELDS.get(m.group(0), lambda _: '')(j)), fmt))
    return 0


def sacct(argv: list[str]) -> int:
    ids, opts = _pick(argv)
    cols = (opts.get('--format') or 'JobID,JobName,State,Elapsed').split(',')
    getters = {'JobID': lambda j: j, 'JobName': lambda j: _read(j, 'name'), 'State': _state,
               'Elapsed': _elapsed, 'ExitCode': lambda j: f"{_read(j, 'exit', '0')}:0"}
    sep = '|' if '--parsable2' in opts else ' '
    if not opts.get('--noheader'):
        print(sep.join(cols))
    for j in ids or _jobs():
        if _dir(j).is_dir():
            print(sep.join(str(getters.get(c, lambda _: '')(j)) for c in cols))
    return 0


def scancel(argv: list[str]) -> int:
    for j in argv:
        pid = _read(j, 'pid')
        if pid and _alive(pid):
            try:
                os.killpg(int(pid), signal.SIGTERM)
            except OSError:
                pass
        if _dir(j).is_dir():
            _write(j, 'state', 'CANCELLED')
            _write(j, 'end', str(time.time()))
    return 0


if __name__ == '__main__':
    tool = Path(sys.argv[1]).name
    sys.exit({'sbatch': sbatch, 'squeue': squeue, 'sacct': sacct, 'scancel': scancel}[tool](sys.argv[2:]))
