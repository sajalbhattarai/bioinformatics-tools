"""
Manages each user's working copies of MARGIE's growing databases on scratch and their backups on depot.
Base copies on depot are only read; a user's first run copies the base (or newest backup) to
<scratch>/margie-2026/<store>/<user>-<name>-vN, and "Back up" copies vN to depot and renames the working copy to vN+1.
Copies run as a SLURM job (or detached on the login node) driven by a plan file; nothing is ever deleted.
"""

from __future__ import annotations

import json
import os
import logging
import posixpath
import re
import shlex
import time

from bioinformatics_tools.utilities import ssh_sftp

LOGGER = logging.getLogger(__name__)

# Shared storage for the group (MARGIE_SHARED_ROOT on the server); empty: none.
DEPOT = os.environ.get('MARGIE_SHARED_ROOT', '').rstrip('/')
# Base copies, one folder per database; each user's backups go beside them.
GENERATED = f'{DEPOT}/databases/margie-generated-databases' if DEPOT else ''
# Folder under the user's scratch that holds every store and run output folder.
ROOT_NAME = 'margie-2026'
# Bookkeeping folder inside it: the store manifest and the copy under way.
STATE_DIR = '.margie'

# ---- Stores ----
# kind 'file' is one file plus companions; kind 'dir' is a folder copied whole.
# config maps config.yaml keys to the store (for a 'dir', to a file inside it; '' is the folder).
STORES: list[dict] = [
    {
        'id': 'job_db',
        'label': 'Job database',
        'note': 'Every run’s results and the job history, in SQLite.',
        'kind': 'file',
        'base': f'{GENERATED}/sqlite/margie-thesis-22-prokaryotes-base.db' if GENERATED else '',
        'depot_sub': 'sqlite',
        'sub': 'sqlite',
        'name': 'margie-thesis-22-prokaryotes-base',
        'ext': '.db',
        'companions': [],
        'config': {'main_database': None},
    },
    {
        'id': 'occ',
        'label': 'Operon reference',
        'note': 'The cross-genome operon (OCC) reference scoring reads and adds each genome to.',
        'kind': 'file',
        'base': f'{GENERATED}/operon-database/occ_reference.pkl' if GENERATED else '',
        'depot_sub': 'operon-database',
        'sub': 'operon-database',
        'name': 'occ_reference',
        'ext': '.pkl',
        # Written beside the reference by update-occ-reference-depot.py.
        'companions': ['.members.tsv', '.genome_stats.tsv'],
        'config': {'margie_sb.operon_database.occ_reference_pkl': None},
    },
    {
        'id': 'fingerprint',
        'label': 'Fingerprint databases',
        'note': 'The gene fingerprint database and the four operon fingerprint databases beside it.',
        'kind': 'dir',
        'base': f'{GENERATED}/fingerprint-database' if GENERATED else '',
        # Only the database files, not the backups and locks beside them.
        'base_files': [
            'fingerprint-database.tsv',
            'fingerprint-database-metadata.json',
            'operon-fingerprint-database-evidence-ordered.tsv',
            'operon-fingerprint-database-evidence-composition.tsv',
            'operon-fingerprint-database-label-ordered.tsv',
            'operon-fingerprint-database-label-composition.tsv',
            'operon-fingerprint-database-evidence_ordered-metadata.json',
            'operon-fingerprint-database-evidence_composition-metadata.json',
            'operon-fingerprint-database-label_ordered-metadata.json',
            'operon-fingerprint-database-label_composition-metadata.json',
        ],
        'depot_sub': 'fingerprint-database',
        'sub': 'fingerprint-database',
        'name': 'fingerprint-database',
        'config': {
            'margie_sb.fingerprint_database.path': 'fingerprint-database.tsv',
            # The report figures read the label-ordered operon fingerprints.
            'margie_sb.report_figures.operon_db': 'operon-fingerprint-database-label-ordered.tsv',
        },
    },
    {
        'id': 'genome_pool',
        'label': 'Genome pool',
        'note': 'Genomes ANI and AAI compare against; every annotated genome joins it.',
        'kind': 'dir',
        'base': f'{GENERATED}/genome-pool' if GENERATED else '',
        # Only the genome folders, not the backups beside them.
        'base_files': ['fna', 'faa'],
        'depot_sub': 'genome-pool',
        'sub': 'genome-pool',
        'name': 'genome-pool',
        'config': {'margie_sb.genome_pool.path': ''},
    },
]

# Run output folders on scratch; never backed up.
OUTPUT_DIRS = {
    'margie_sb.scoring_results_historical.path': 'scoring-archive',
    'margie_sb.final_tables_depot.path': 'final-tables',
    'margie_sb.sqlite_pipeline_snapshot.path': 'sqlite/snapshots',
}

STORE_BY_ID = {s['id']: s for s in STORES}


def backup_root(cfg: dict | None) -> str:
    """Returns margie_sb.backup_root if set, else the base folder on depot."""
    chosen = _cfg_get(cfg or {}, 'margie_sb.backup_root')
    return chosen.strip().rstrip('/') if isinstance(chosen, str) and chosen.strip() else GENERATED


def backup_dir(store: dict, cfg: dict | None) -> str:
    """Returns the store's backup folder."""
    return f"{backup_root(cfg)}/{store['depot_sub']}"


class StoreError(Exception):
    """A request that cannot be done now; the message and HTTP status go to the page."""

    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


# ---- Remote helpers ----

def _run(conn, command: str, timeout: float = 60.0) -> tuple[int, str]:
    """Runs a shell command over SSH; returns (exit code, stdout + stderr)."""
    ssh = conn.connect()
    _, stdout, stderr = ssh.exec_command(command, timeout=timeout)
    code = stdout.channel.recv_exit_status()
    out = stdout.read().decode(errors='replace') + stderr.read().decode(errors='replace')
    return code, out


def _cfg_get(cfg: dict, dotted: str):
    """Returns the value at a dotted config key, or None."""
    cur = cfg
    for part in dotted.split('.'):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _cfg_set(cfg: dict, dotted: str, value) -> None:
    """Sets a dotted config key, creating intermediate dicts."""
    parts = dotted.split('.')
    cur = cfg
    for part in parts[:-1]:
        if not isinstance(cur.get(part), dict):
            cur[part] = {}
        cur = cur[part]
    cur[parts[-1]] = value


def scratch_root(conn, cfg: dict) -> str:
    """Returns margie_sb.stores_root if set, else /scratch/<cluster>/<user>/margie-2026.
    The cluster short name comes from the login node's domain (login01.mycluster.example.edu -> mycluster).
    """
    chosen = _cfg_get(cfg, 'margie_sb.stores_root')
    if isinstance(chosen, str) and chosen.strip():
        return chosen.strip().rstrip('/')
    code, out = _run(conn, 'printf "%s %s" "$(hostname -d 2>/dev/null | cut -d. -f1)" "$USER"')
    cluster, _, user = out.strip().partition(' ')
    if code != 0 or not cluster or not user:
        raise StoreError('Could not work out your scratch folder on the cluster. Set margie_sb.stores_root in Settings.')
    return f'/scratch/{cluster}/{user}/{ROOT_NAME}'


# ---- Names ----

def versioned(store: dict, user: str, version: int) -> str:
    """Returns <user>-<name>-v<N><ext>, the name used on scratch and depot."""
    return f"{user}-{store['name']}-v{version}{store.get('ext', '') if store['kind'] == 'file' else ''}"


def _version_of(store: dict, user: str, name: str) -> int | None:
    """Parses the version from a versioned name, or None."""
    ext = re.escape(store.get('ext', '')) if store['kind'] == 'file' else ''
    m = re.fullmatch(rf"{re.escape(user)}-{re.escape(store['name'])}-v(\d+){ext}", name)
    return int(m.group(1)) if m else None


def _backups(conn, store: dict, user: str, cfg: dict | None = None) -> list[int]:
    """Returns the versions of this user's backups of the store, oldest first."""
    try:
        entries = ssh_sftp.list_remote_dir(backup_dir(store, cfg), connection=conn)
    except Exception:
        return []
    want = 'file' if store['kind'] == 'file' else 'directory'
    found = []
    for e in entries:
        if e.get('type') != want:
            continue
        v = _version_of(store, user, e.get('name') or '')
        if v is not None:
            found.append(v)
    return sorted(found)


# ---- State: <root>/.margie/stores.json holds working versions; op.* is the current or last copy ----

def _state_path(root: str, name: str) -> str:
    return f'{root}/{STATE_DIR}/{name}'


def _read_text(conn, path: str) -> str:
    code, out = _run(conn, f'cat {shlex.quote(path)} 2>/dev/null')
    return out if code == 0 else ''


def read_manifest(conn, root: str) -> dict:
    try:
        data = json.loads(_read_text(conn, _state_path(root, 'stores.json')) or '{}')
        return data if isinstance(data, dict) else {}
    except json.JSONDecodeError:
        return {}


def _write_manifest(conn, root: str, manifest: dict) -> None:
    ssh_sftp.write_remote_text_file(_state_path(root, 'stores.json'), json.dumps(manifest, indent=2) + '\n', connection=conn)


def working_path(store: dict, root: str, user: str, version: int) -> str:
    return f"{root}/{store['sub']}/{versioned(store, user, version)}"


def _probe(conn, root: str, user: str, cfg: dict | None = None) -> dict:
    """Collects the manifest, working copies, backups and copy status in one cluster round-trip
    (plus a small second one for the working copies)."""
    manifest_path = _state_path(root, 'stores.json')
    status_path = _state_path(root, 'op.status')
    log_path = _state_path(root, 'op.log')
    parts = [f'echo "@@manifest"; cat {shlex.quote(manifest_path)} 2>/dev/null',
             f'echo "@@status"; cat {shlex.quote(status_path)} 2>/dev/null',
             f'echo "@@log"; tail -c 16000 {shlex.quote(log_path)} 2>/dev/null | tr "\\r" "\\n"',
             'echo "@@host"; hostname',
             f'echo "@@squeue"; j=$(cat {shlex.quote(_state_path(root, "op.jobid"))} 2>/dev/null); '
             f'[ -n "$j" ] && {{ echo "job=$j"; squeue -h -j "$j" -o "%T|%r" 2>/dev/null; }}',
             f'echo "@@slurmout"; tail -n 20 {shlex.quote(_state_path(root, "op.slurm.out"))} 2>/dev/null',
             f'echo "@@logage"; echo $(( $(date +%s) - $(stat -c %Y {shlex.quote(log_path)} 2>/dev/null || echo 0) ))']
    for st in STORES:
        parts.append(f'echo "@@depot {st["id"]}"; ls -1p {shlex.quote(backup_dir(st, cfg))} 2>/dev/null | grep "^{user}-{st["name"]}-v"')
    code, out = _run(conn, '; '.join(parts), timeout=60)
    sections: dict[str, list[str]] = {}
    current = None
    for line in out.splitlines():
        if line.startswith('@@'):
            current = line[2:].strip()
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    try:
        manifest = json.loads('\n'.join(sections.get('manifest', [])) or '{}')
    except json.JSONDecodeError:
        manifest = {}
    # Checks the working copies named by the manifest.
    wanted = {}
    for st in STORES:
        v = (manifest.get(st['id']) or {}).get('version') if isinstance(manifest, dict) else None
        if v:
            wanted[st['id']] = (int(v), working_path(st, root, user, int(v)))
    present = set()
    if wanted:
        cmd = '; '.join(f'test -e {shlex.quote(p)} && echo {sid}' for sid, (_, p) in wanted.items())
        _, out2 = _run(conn, cmd)
        present = set(out2.split())
    backups = {}
    for st in STORES:
        names = [n.rstrip('/') for n in sections.get(f'depot {st["id"]}', [])]
        is_dir = st['kind'] == 'dir'
        vs = sorted(v for n in names
                    if (v := _version_of(st, user, n)) is not None
                    and any(x == n + '/' for x in sections.get(f'depot {st["id"]}', [])) == is_dir)
        backups[st['id']] = vs
    return {
        'manifest': manifest if isinstance(manifest, dict) else {},
        'wanted': wanted,
        'present': present,
        'backups': backups,
        'op': _op_from(sections, root),
    }


def _op_from(sections: dict, root: str) -> dict:
    """Builds the current or last copy's state from its status file, SLURM state and log progress."""
    op: dict = {}
    for line in sections.get('status', []):
        k, _, v = line.partition('=')
        if k:
            op[k] = v
    if not op:
        return {}
    # A running op whose script is gone is failed: checked by squeue, by pid on the same
    # login node, or by a log idle for five minutes on another node.
    squeue = [l.strip() for l in sections.get('squeue', []) if l.strip()]
    job = next((l[4:] for l in squeue if l.startswith('job=')), '')
    slurm_line = next((l for l in squeue if not l.startswith('job=')), '')
    slurm_state, _, slurm_reason = slurm_line.partition('|')
    if job:
        op['job'] = job
        op['slurm_state'] = slurm_state
        # Pending reason (Resources, Priority, ...) shown beside the queued state.
        if slurm_state == 'PENDING' and slurm_reason and slurm_reason != 'None':
            op['slurm_reason'] = slurm_reason
    if op.get('state') in ('queued', 'running') and job:
        # A SLURM job missing from squeue has ended.
        if not slurm_state:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy job ended before it finished.'
            op['slurm_out'] = sections.get('slurmout', [])[-10:]
    elif op.get('state') == 'running':
        here = (sections.get('host') or [''])[0].strip()
        age = (sections.get('logage') or ['0'])[0].strip()
        if here and here == op.get('host') and op.get('pid', '').isdigit():
            op['_check_pid'] = op['pid']
        elif age.isdigit() and int(age) > 300:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy stopped before it finished.'
    lines = [l.rstrip() for l in sections.get('log', []) if l.strip()]
    # rsync --info=progress2 prints "   1,234,567  42%  ..." for the file in hand.
    pct = 0.0
    for l in reversed(lines):
        m = re.search(r'\s(\d{1,3})%\s', f' {l} ')
        if m:
            pct = float(m.group(1))
            break
    total = float(op.get('total_bytes') or 0)
    done = float(op.get('done_bytes') or 0)
    item = float(op.get('item_bytes') or 0)
    if op.get('state') == 'done':
        percent = 100.0
    elif total > 0:
        percent = min(100.0, 100.0 * (done + item * pct / 100.0) / total)
    else:
        percent = pct
    op['percent'] = round(percent, 1)
    # Drops rsync progress lines from the log shown.
    op['log'] = [l for l in lines if not re.match(r'^\s*[\d,]+\s+\d{1,3}%', l)][-40:]
    return op


def _read_op(conn, root: str) -> dict:
    """Returns the current or last copy's state, with its pid checked."""
    op = _probe(conn, root, '__none__')['op']
    return _settle_pid(conn, op)


def _settle_pid(conn, op: dict) -> dict:
    """Marks the op failed if its pid (checked on its own login node) is gone."""
    pid = op.pop('_check_pid', None)
    if pid:
        code, _ = _run(conn, f'kill -0 {int(pid)} 2>/dev/null')
        if code != 0:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy stopped before it finished.'
    return op


# ---- Copy script: runs a plan of tab-separated steps and writes op.status/op.log ----
#   copy<TAB>label<TAB>src<TAB>dst<TAB>bytes   rsync src -> dst (via dst.partial)
#   files<TAB>label<TAB>srcdir<TAB>dst<TAB>bytes<TAB>name,name,...
#   move<TAB>src<TAB>dst                        rename (same filesystem: instant)
#   mkdir<TAB>path
#   link<TAB>src<TAB>dst                        symlink dst -> src, unless dst exists
#   pull<TAB>label<TAB>url<TAB>dst              apptainer pull url -> dst (via dst.partial)
#   repo<TAB>folder                             margie-build's folder on the cluster, for build
#   build<TAB>label<TAB>mode<TAB>tool<TAB>sifdir<TAB>dbdir<TAB>statement
#                                               margie-build's build.sh --<mode> <tool>;
#                                               a statement accepts a gated tool's licence
# (link, pull, repo and build are used by tool_assets.py.)
COPY_SCRIPT = r'''#!/bin/bash
set -u
dir="$1"; plan="$dir/op.plan"; status="$dir/op.status"; log="$dir/op.log"
total=$(awk -F'\t' '$1=="copy"||$1=="files"{s+=$5} END{print s+0}' "$plan")
done_bytes=0
say() { { printf 'state=%s\npid=%s\nhost=%s\nop=%s\nlabel=%s\ntotal_bytes=%s\ndone_bytes=%s\nitem_bytes=%s\nmessage=%s\n' \
          "$1" "$$" "$(hostname)" "$OP" "$2" "$total" "$done_bytes" "$3" "$4"; } > "$status.tmp" && mv "$status.tmp" "$status"; }
fail() { say failed "$1" 0 "$2"; echo "FAILED: $2" >> "$log"; exit 1; }
OP="$(head -n 1 "$dir/op.kind" 2>/dev/null)"
: > "$log"
say running "Starting" 0 ""
while IFS=$'\t' read -r kind a b c d e f; do
  case "$kind" in
    copy)
      label="$a"; src="$b"; dst="$c"; bytes="$d"
      say running "$label" "$bytes" ""
      echo "== $label: $src -> $dst" >> "$log"
      mkdir -p "$(dirname "$dst")" || fail "$label" "could not make $(dirname "$dst")"
      if [ -d "$src" ]; then
        rsync -a --info=progress2 --no-inc-recursive "$src/" "$dst.partial/" >> "$log" 2>&1 || fail "$label" "copy of $src failed"
      else
        rsync -a --info=progress2 "$src" "$dst.partial" >> "$log" 2>&1 || fail "$label" "copy of $src failed"
      fi
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      done_bytes=$(( done_bytes + bytes ))
      ;;
    files)
      label="$a"; src="$b"; dst="$c"; bytes="$d"; names="$e"
      say running "$label" "$bytes" ""
      echo "== $label: $src/{$names} -> $dst" >> "$log"
      mkdir -p "$dst.partial" || fail "$label" "could not make $dst"
      list=(); IFS=',' read -r -a want <<< "$names"
      for n in "${want[@]}"; do [ -e "$src/$n" ] && list+=("$src/$n"); done
      if [ ${#list[@]} -gt 0 ]; then
        rsync -a --info=progress2 "${list[@]}" "$dst.partial/" >> "$log" 2>&1 || fail "$label" "copy from $src failed"
      fi
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      done_bytes=$(( done_bytes + bytes ))
      ;;
    move)
      echo "== rename $a -> $b" >> "$log"
      mv -T "$a" "$b" || fail "Renaming" "could not rename $a"
      ;;
    mkdir)
      mkdir -p "$a" || fail "Folders" "could not make $a"
      ;;
    link)
      [ -e "$b" ] || [ -L "$b" ] || ln -s "$a" "$b" || fail "Links" "could not link $b"
      ;;
    pull)
      label="$a"; url="$b"; dst="$c"
      say running "$label" 0 ""
      echo "== $label: $url -> $dst" >> "$log"
      app=$(command -v apptainer || command -v singularity) || fail "$label" "apptainer is not installed here"
      mkdir -p "$(dirname "$dst")" || fail "$label" "could not make $(dirname "$dst")"
      # Its cache and scratch space beside this script, not in the small home folder.
      export APPTAINER_CACHEDIR="$dir/apptainer-cache" APPTAINER_TMPDIR="$dir/apptainer-tmp"
      mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
      "$app" pull --force "$dst.partial" "$url" >> "$log" 2>&1 || fail "$label" "could not pull $url"
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      ;;
    repo)
      repo_url="$a"
      ;;
    build)
      label="$a"; mode="$b"; tool="$c"; sifdir="$d"; dbdir="$e"; statement="$f"
      say running "$label" 0 ""
      echo "== $label" >> "$log"
      # margie-build's own folder on the cluster (the lab's, or the user's clone),
      # used where it is: build.sh writes only to --sif-dir, --db-dir, its
      # --work-dir and the licence records. --local: this is already a SLURM
      # job, and build.sh must not submit another. CONFIG_FILE is an empty file of our own: whoever's
      # .build-config.sh is in that folder must not decide where this user's
      # builds go (build.sh looks for one unless CONFIG_FILE names a file).
      repo="${repo_url:-}"
      [ -f "$repo/build.sh" ] || fail "$label" "margie-build is not at ${repo:-any folder}; set it up on Install"
      : > "$dir/build-config.sh"
      accept=()
      [ -n "$statement" ] && accept=(--accept-"$tool"-licence --licence-statement "$statement")
      mkdir -p "$sifdir" "$dbdir" "$dir/apptainer-cache" "$dir/apptainer-tmp" || fail "$label" "could not make $sifdir or $dbdir"
      ( cd "$repo" && CONFIG_FILE="$dir/build-config.sh" APPTAINER_CACHEDIR="$dir/apptainer-cache" APPTAINER_TMPDIR="$dir/apptainer-tmp" \
          THREADS="${SLURM_CPUS_PER_TASK:-4}" ./build.sh --"$mode" "$tool" --local --work-dir "$dir/build-work" --sif-dir "$sifdir" --db-dir "$dbdir" \
          --audit-log "$dir/licence-acceptances.tsv" --no-color --quiet-notice ${accept[@]+"${accept[@]}"} ) >> "$log" 2>&1 \
        || fail "$label" "margie-build could not build $tool (the log above says why)"
      ;;
  esac
done < "$plan"
say done "Finished" 0 ""
echo "== finished" >> "$log"
'''


def _copy_cpus(cfg: dict | None) -> int:
    """Returns the copy job's cores: margie_sb.stores_copy_cpus, default 4, capped at 16 (the lab's limit)."""
    try:
        n = int(_cfg_get(cfg or {}, 'margie_sb.stores_copy_cpus') or 4)
    except (TypeError, ValueError):
        n = 4
    return max(1, min(16, n))


def _start(conn, root: str, kind: str, plan: list[list], after: dict, cfg: dict | None = None,
           walltime: str = '12:00:00', mem: str = '4G') -> None:
    """Writes the plan and starts the copy script: via sbatch when an account is set, else detached on the login node."""
    state = f'{root}/{STATE_DIR}'
    code, out = _run(conn, f'mkdir -p {shlex.quote(state)} && chmod 700 {shlex.quote(state)}')
    if code != 0:
        raise StoreError(f'Could not make {state}: {out.strip()}')
    ssh_sftp.write_remote_text_file(f'{state}/op.sh', COPY_SCRIPT, connection=conn)
    ssh_sftp.write_remote_text_file(
        f'{state}/op.plan', ''.join('\t'.join(str(x) for x in step) + '\n' for step in plan), connection=conn)
    ssh_sftp.write_remote_text_file(f'{state}/op.kind', kind + '\n', connection=conn)
    # Recorded by apply_finished once the copy is done.
    ssh_sftp.write_remote_text_file(f'{state}/op.after.json', json.dumps(after) + '\n', connection=conn)
    _run(conn, f'cd {shlex.quote(state)} && rm -f op.status op.log op.jobid op.slurm.out')
    account = str(_cfg_get(cfg or {}, 'compute.cluster_default.account') or '').strip()
    partition = str(_cfg_get(cfg or {}, 'compute.cluster_default.partition') or '').strip()
    q = shlex.quote
    if account:
        # Shows the job as queued until it starts.
        ssh_sftp.write_remote_text_file(
            f'{state}/op.status', f'state=queued\nop={kind}\nlabel=Waiting for a SLURM slot\n', connection=conn)
        cmd = (f'cd {q(state)} && sbatch --parsable --job-name=margie-databases --account={q(account)} '
               + (f'--partition={q(partition)} ' if partition else '')
               + f'--time={q(walltime)} --ntasks=1 --cpus-per-task={_copy_cpus(cfg)} --mem={q(mem)} --output={q(state)}/op.slurm.out '
               + f'op.sh {q(state)}')
        code, out = _run(conn, cmd)
        job = out.strip().split(';')[0].split()[-1] if out.strip() else ''
        if code != 0 or not job.isdigit():
            raise StoreError(f'SLURM would not take the copy job: {out.strip()}')
        ssh_sftp.write_remote_text_file(f'{state}/op.jobid', job + '\n', connection=conn)
        return
    code, out = _run(conn, f'cd {q(state)} && setsid nohup bash op.sh {q(state)} '
                            f'> /dev/null 2>&1 < /dev/null & echo started')
    if 'started' not in out:
        raise StoreError(f'Could not start the copy: {out.strip()}')
    # Waits briefly for the status file so the op is never seen in no state.
    for _ in range(20):
        if _read_text(conn, f'{state}/op.status'):
            break
        time.sleep(0.25)


def _size(conn, paths: list[str]) -> int:
    """Returns the total byte size of paths via du."""
    if not paths:
        return 0
    code, out = _run(conn, 'du -sbc ' + ' '.join(shlex.quote(p) for p in paths) + ' 2>/dev/null | tail -1 | cut -f1',
                     timeout=300)
    try:
        return int(out.strip() or 0)
    except ValueError:
        return 0


# ---- Page requests ----

def _config_paths(store: dict, target: str) -> dict:
    """Returns the config values for a store whose working copy is at target."""
    return {key: (posixpath.join(target, inner) if inner else target) if store['kind'] == 'dir' else target
            for key, inner in store['config'].items()}


def status(conn, cfg: dict, user: str) -> dict:
    """Returns each store's working version, latest backup and the copy under way."""
    root = scratch_root(conn, cfg)
    probe = _probe(conn, root, user, cfg)
    op = _settle_pid(conn, probe['op'])
    stores = []
    ready = True
    for s in STORES:
        v, path = probe['wanted'].get(s['id'], (None, ''))
        present = s['id'] in probe['present']
        ready = ready and present
        backups = probe['backups'].get(s['id'], [])
        stores.append({
            'id': s['id'],
            'label': s['label'],
            'note': s['note'],
            'version': v if present else None,
            'path': path if present else '',
            'base': s['base'],
            'backup': ({'version': backups[-1],
                        'path': posixpath.join(backup_dir(s, cfg), versioned(s, user, backups[-1]))}
                       if backups else None),
            'backups': len(backups),
        })
    return {'root': root, 'user': user, 'ready': ready, 'stores': stores, 'op': op or None}


def progress(conn, cfg: dict) -> dict | None:
    """Returns only the copy under way, for the progress bar's polling."""
    root = scratch_root(conn, cfg)
    return _read_op(conn, root) or None


def run_here(conn, cfg: dict) -> dict | None:
    """Cancels a queued SLURM copy job and runs the same script detached on the login node."""
    root = scratch_root(conn, cfg)
    op = _read_op(conn, root)
    if not op or op.get('state') != 'queued' or not op.get('job'):
        raise StoreError('There is no copy waiting in the queue.', 409)
    if op.get('slurm_state') != 'PENDING':
        raise StoreError('The copy has already started on a compute node.', 409)
    state = f'{root}/{STATE_DIR}'
    q = shlex.quote
    code, out = _run(conn, f'scancel {q(op["job"])} && cd {q(state)} && rm -f op.jobid op.slurm.out op.status '
                            f'&& (setsid nohup bash op.sh {q(state)} > /dev/null 2>&1 < /dev/null &) && echo started')
    if 'started' not in out:
        raise StoreError(f'Could not start the copy here: {out.strip()}')
    for _ in range(20):
        if _read_text(conn, f'{state}/op.status'):
            break
        time.sleep(0.25)
    return _read_op(conn, root) or None


def _busy(op: dict | None) -> bool:
    """Returns True while a copy is queued or running."""
    return bool(op) and op.get('state') in ('queued', 'running')


def start_setup(conn, cfg: dict, user: str) -> dict:
    """Copies every store missing on scratch from the user's newest backup, or else the base."""
    st = status(conn, cfg, user)
    if _busy(st['op']):
        raise StoreError('A copy is already under way.', 409)
    root = st['root']
    plan: list[list] = [['mkdir', f'{root}/{sub}'] for sub in OUTPUT_DIRS.values()]
    after: dict = {'kind': 'setup', 'versions': {}}
    for s, info in zip(STORES, st['stores']):
        if info['version']:
            continue
        backups = _backups(conn, s, user, cfg)
        if backups:
            src = posixpath.join(backup_dir(s, cfg), versioned(s, user, backups[-1]))
            version = backups[-1] + 1
            label = f"{s['label']} (from your backup v{backups[-1]})"
        else:
            src = s['base']
            if not src:
                raise StoreError('No base databases to copy from: set MARGIE_SHARED_ROOT on the server, '
                                 'or restore from a backup of your own.', 409)
            version = 1
            label = s['label']
        dst = working_path(s, root, user, version)
        if s['kind'] == 'file':
            plan.append(['copy', label, src, dst, _size(conn, [src])])
            for c in s['companions']:
                code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
                if code == 0:
                    plan.append(['copy', f"{s['label']} ({c.lstrip('.')})", src + c, dst + c, _size(conn, [src + c])])
        elif s.get('base_files') and src == s['base']:
            files = [posixpath.join(src, f) for f in s['base_files']]
            plan.append(['files', label, src, dst, _size(conn, files), ','.join(s['base_files'])])
        else:
            plan.append(['copy', label, src, dst, _size(conn, [src])])
        after['versions'][s['id']] = version
    _start(conn, root, 'setup', plan, after, cfg)
    return status(conn, cfg, user)


_UNITS = {'B': 1, 'KB': 1024, 'MB': 1024 ** 2, 'GB': 1024 ** 3, 'TB': 1024 ** 4, 'PB': 1024 ** 5,
          'K': 1024, 'M': 1024 ** 2, 'G': 1024 ** 3, 'T': 1024 ** 4, 'P': 1024 ** 5}


def _bytes(text: str) -> int | None:
    """Parses a size such as '1.5TB' or '20G' into bytes."""
    m = re.fullmatch(r'([\d.]+)\s*([KMGTP]?B?)', text.strip(), re.I)
    return int(float(m.group(1)) * _UNITS.get(m.group(2).upper() or 'B', 1)) if m else None


def depot_free(conn, where: str = DEPOT) -> int | None:
    """Returns the free bytes on depot: the smaller of the group's quota headroom (myquota) and df."""
    group = where.split('/')[2] if where.startswith('/depot/') and len(where.split('/')) > 2 else ''
    # Runs df on the nearest existing folder, since the backup folder may not exist yet.
    code, out = _run(conn, f'myquota 2>/dev/null; echo "@@df"; d={shlex.quote(where)}; '
                           f'while [ ! -e "$d" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; '
                           f'df -B1 --output=avail "$d" 2>/dev/null | tail -1')
    quota_free = None
    head, _, df = out.partition('@@df')
    for line in head.splitlines():
        cols = line.split()
        if group and len(cols) >= 4 and cols[0] == 'depot' and cols[1] == group:
            used, limit = _bytes(cols[2]), _bytes(cols[3])
            if used is not None and limit is not None:
                quota_free = max(0, limit - used)
    try:
        fs_free = int(df.strip())
    except ValueError:
        fs_free = None
    known = [x for x in (quota_free, fs_free) if x is not None]
    return min(known) if known else None


def backup_check(conn, cfg: dict, user: str, store_id: str) -> dict:
    """Returns a backup's size, depot's free space and whether it fits."""
    s = STORE_BY_ID.get(store_id)
    if not s:
        raise StoreError('No such database.', 404)
    st = status(conn, cfg, user)
    info = next(i for i in st['stores'] if i['id'] == store_id)
    if not info['version']:
        raise StoreError(f"{s['label']} is not set up on scratch yet.", 409)
    src = info['path']
    paths = [src] + [src + c for c in (s.get('companions') or [])]
    size = _size(conn, paths)
    free = depot_free(conn, backup_root(cfg))
    # Keeps 5 GB free on the shared depot.
    margin = 5 * 1024 ** 3
    fits = free is None or free - margin >= size
    return {
        'id': store_id,
        'label': s['label'],
        'version': info['version'],
        'size': size,
        'free': free,
        'fits': fits,
        'target': posixpath.join(backup_dir(s, cfg), versioned(s, user, info['version'])),
    }


def start_backup(conn, cfg: dict, user: str, store_id: str) -> dict:
    """Copies the working vN to depot as <user>-<name>-vN, then renames the working copy to vN+1."""
    s = STORE_BY_ID.get(store_id)
    if not s:
        raise StoreError('No such database.', 404)
    st = status(conn, cfg, user)
    if _busy(st['op']):
        raise StoreError('A copy is already under way.', 409)
    if not backup_root(cfg):
        raise StoreError('No backup folder: set margie_sb.backup_root in Settings.', 409)
    check = backup_check(conn, cfg, user, store_id)
    if not check['fits']:
        gb = lambda n: f'{n / 1024 ** 3:.1f} GB'
        raise StoreError(f"Depot has {gb(check['free'])} left and this backup needs {gb(check['size'])} "
                         '(with 5 GB kept free): it would not fit.', 409)
    info = next(i for i in st['stores'] if i['id'] == store_id)
    v = info['version']
    if not v:
        raise StoreError(f"{s['label']} is not set up on scratch yet.", 409)
    root = st['root']
    src = working_path(s, root, user, v)
    dst = posixpath.join(backup_dir(s, cfg), versioned(s, user, v))
    code, _ = _run(conn, f'test -e {shlex.quote(dst)}')
    if code == 0:
        raise StoreError(f'{dst} already exists; nothing is overwritten.', 409)
    nxt = working_path(s, root, user, v + 1)
    plan: list[list] = [['copy', f"{s['label']} v{v} to depot", src, dst, _size(conn, [src])]]
    companions = s.get('companions', []) if s['kind'] == 'file' else []
    for c in companions:
        code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
        if code == 0:
            plan.append(['copy', f"{s['label']} ({c.lstrip('.')})", src + c, dst + c, _size(conn, [src + c])])
    # The rename runs only after the backup copy succeeds.
    plan.append(['move', src, nxt])
    for c in companions:
        code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
        if code == 0:
            plan.append(['move', src + c, nxt + c])
    _start(conn, root, f'backup:{store_id}', plan, {'kind': 'backup', 'versions': {store_id: v + 1}}, cfg)
    return status(conn, cfg, user)


def apply_finished(conn, cfg: dict, user: str) -> bool:
    """Records a finished copy's new versions and config changes; idempotent.
    Returns True when the config changed (the caller saves it)."""
    root = scratch_root(conn, cfg)
    op = _read_op(conn, root)
    if op.get('state') != 'done':
        return False
    state = f'{root}/{STATE_DIR}'
    after_text = _read_text(conn, f'{state}/op.after.json')
    changed = False
    if after_text:
        try:
            after = json.loads(after_text)
        except json.JSONDecodeError:
            after = {}
        manifest = read_manifest(conn, root)
        for sid, version in (after.get('versions') or {}).items():
            manifest[sid] = {'version': int(version)}
        _write_manifest(conn, root, manifest)
        # Config keys set by a tool-assets copy.
        for key, value in (after.get('config') or {}).items():
            if _cfg_get(cfg, key) != value:
                _cfg_set(cfg, key, value)
                changed = True
        _run(conn, f'mv -f {shlex.quote(state)}/op.after.json {shlex.quote(state)}/op.applied.json')
    if op.get('op') == 'assets':
        # A tool-assets copy leaves the store settings alone.
        return changed
    return point_config(conn, cfg, user, root) or changed


def point_config(conn, cfg: dict, user: str, root: str | None = None) -> bool:
    """Sets every store's and output folder's config key to its scratch path."""
    root = root or scratch_root(conn, cfg)
    manifest = read_manifest(conn, root)
    changed = False
    wanted: dict = {}
    for s in STORES:
        v = (manifest.get(s['id']) or {}).get('version')
        if v:
            wanted.update(_config_paths(s, working_path(s, root, user, v)))
    for key, sub in OUTPUT_DIRS.items():
        wanted[key] = f'{root}/{sub}'
    # Cached so later requests skip the cluster lookup.
    wanted['margie_sb.stores_root'] = root
    for key, value in wanted.items():
        if _cfg_get(cfg, key) != value:
            _cfg_set(cfg, key, value)
            changed = True
    return changed


def is_ready(conn, cfg: dict, user: str) -> bool:
    """Returns True when every store has a working copy on scratch."""
    return status(conn, cfg, user)['ready']
