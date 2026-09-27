"""
Per-protein annotation cache in the job's SQLite database, keyed by the protein's
sequence hash and a tool_key (container, database and settings).

split() writes novel.faa plus the cached rows before a tool runs, MERGE_SH joins
both into the tool's results file, and store() adds the genome's new proteins.
Only tools whose rows depend on the protein alone are listed in TOOLS.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

try:
    from bioinformatics_tools.workflow_tools.genome_identity import protein_hash, read_fasta
    from bioinformatics_tools.workflow_tools.workflow_helpers import db_path as tool_db_path, rc, sif_path
except ImportError:  # imported from margie_sb.smk, with this folder on sys.path
    from genome_identity import protein_hash, read_fasta
    from workflow_helpers import db_path as tool_db_path, rc, sif_path

LOGGER = logging.getLogger(__name__)

SCHEMA = 1
WORKFLOW = 'margie_sb'
PC_DIR = 'protein-cache'

CREATE_SQL = """
CREATE TABLE IF NOT EXISTS protein_cache (
    aa_hash TEXT NOT NULL,
    tool TEXT NOT NULL,
    tool_key TEXT NOT NULL,
    filename TEXT NOT NULL,
    lines TEXT NOT NULL,
    cached_at TEXT NOT NULL,
    PRIMARY KEY (tool, tool_key, filename, aa_hash)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS protein_cache_header (
    tool TEXT NOT NULL,
    tool_key TEXT NOT NULL,
    filename TEXT NOT NULL,
    header TEXT NOT NULL,
    PRIMARY KEY (tool, tool_key, filename)
) WITHOUT ROWID;
"""

# InterPro's member databases, as margie_sb.smk names them.
INTERPRO_ALL_ANALYSES = {
    "AntiFam": "antifam", "CDD": "cdd", "Coils": "coils", "FunFam": "funfam",
    "Gene3D": "gene3d", "Hamap": "hamap", "MobiDBLite": "mobidb", "NCBIfam": "ncbifam",
    "PANTHER": "panther", "Pfam": "pfam", "PIRSF": "pirsf", "PIRSR": "pirsr",
    "PRINTS": "prints", "ProSitePatterns": "prosite_patterns", "ProSiteProfiles": "prosite_profiles",
    "SFLD": "sfld", "SMART": "smart", "SUPERFAMILY": "superfamily",
}
INTERPRO_DEFAULT_ANALYSES = ["Hamap", "NCBIfam", "CDD", "PIRSF"]

# Tools whose rows depend only on the protein; each tool's first file marks a protein as cached.
_SIMPLE = ['cog', 'kegg', 'eggnog', 'uniprot', 'pfam', 'merops', 'tcdb', 'dbcan', 'pgap', 'tmbed', 'signalp6']
TOOLS = _SIMPLE + ['tigrfam', 'phobius', 'interpro']

# Rule params that change a tool's answer, with margie_sb.smk's defaults.
TOOL_PARAMS = {
    'cog': [('cog.evalue', '1e-2')],
    'merops': [('merops.evalue', '1e-5')],
    'tcdb': [('tcdb.evalue', '1e-5'), ('tcdb.pct_id', '30')],
    'uniprot': [('uniprot.evalue', '1e-5'), ('uniprot.pct_id', '30')],
}
# Tools run from cluster modules rather than a container.
MODULE_TOOLS = {'signalp6': 'biocontainers/default + signalp6/6.0-fast --organism other --mode fast'}
NO_DATABASE = {'phobius', 'signalp6'}


def interpro_analyses(cfg: dict) -> list[str]:
    return list(rc('interpro.analyses', INTERPRO_DEFAULT_ANALYSES, config=cfg))


def tool_files(tool: str, cfg: dict) -> list[tuple[str, str]]:
    """Returns (file name, kind) for each results file of the tool, marker file first.

    kind is tsv, tsv? (may be missing or empty) or domtbl (HMMER --domtblout)."""
    if tool in _SIMPLE:
        return [(f'{tool}_results.tsv', 'tsv')]
    if tool == 'tigrfam':
        return [('tigrfam_results.tsv', 'tsv'), ('tigrfam_domtbl.out', 'domtbl')]
    if tool == 'phobius':
        return [('phobius_results.tsv', 'tsv'), ('phobius_top1.tsv', 'tsv')]
    if tool == 'interpro':
        return [('interpro_results.tsv', 'tsv')] + [
            (f'interpro_{INTERPRO_ALL_ANALYSES[a]}_results.tsv', 'tsv?')
            for a in interpro_analyses(cfg) if a in INTERPRO_ALL_ANALYSES]
    raise KeyError(tool)


def _stat(path: str) -> list:
    try:
        st = os.stat(path)
        return [os.path.basename(path.rstrip('/')), st.st_size if not os.path.isdir(path) else 0, int(st.st_mtime)]
    except OSError:
        return [path, None, None]


def tool_key(tool: str, cfg: dict) -> str:
    """Hashes what the tool's answer depends on besides the protein: container, database, settings."""
    parts: dict = {'schema': SCHEMA, 'tool': tool}
    parts['runs'] = MODULE_TOOLS.get(tool) or _stat(sif_path(f'{tool}.sif', config=cfg, workflow_id=WORKFLOW))
    if tool not in NO_DATABASE:
        db = tool_db_path(tool, config=cfg, workflow_id=WORKFLOW)
        parts['db'] = [db] + _stat(db)[1:]
    parts['params'] = {k: str(rc(k, d, config=cfg)) for k, d in TOOL_PARAMS.get(tool, [])}
    if tool == 'interpro':
        analyses = interpro_analyses(cfg)
        parts['params']['analyses'] = analyses
        parts['params']['applications'] = str(rc('interpro.applications', ','.join(analyses), config=cfg))
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:32]


def enabled(cfg: dict) -> bool:
    v = rc(f'{WORKFLOW}.protein_cache', True, config=cfg)
    return str(v).strip().lower() not in ('false', '0', 'no', 'off')


def pc_dir(prefix: str, tool: str) -> str:
    return f'{prefix}{tool}/{PC_DIR}'


def clear(prefix: str, tool: str) -> None:
    """Removes a split so the tool's rule reads rast.faa again."""
    shutil.rmtree(pc_dir(prefix, tool), ignore_errors=True)


# ---- database ----

def _connect(db: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db, timeout=120)
    conn.execute('PRAGMA busy_timeout=120000')
    conn.executescript(CREATE_SQL)
    return conn


def _chunks(items: list, n: int = 500):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def _cached_hashes(conn, tool: str, key: str, filename: str, hashes: list[str]) -> set[str]:
    found: set[str] = set()
    for chunk in _chunks(hashes):
        found.update(r[0] for r in conn.execute(
            'SELECT aa_hash FROM protein_cache WHERE tool = ? AND tool_key = ? AND filename = ? '
            f'AND aa_hash IN ({",".join("?" * len(chunk))})', (tool, key, filename, *chunk)))
    return found


# ---- rows: remapped to this genome's feature ids and name ----

_DOMTBL_ID = re.compile(r'^(\S+\s+\S+\s+\S+\s+)(\S+)')


def _columns(header: str) -> tuple[int | None, int | None]:
    cols = header.rstrip('\n').split('\t')
    return (cols.index('feature_id') if 'feature_id' in cols else None,
            cols.index('organism_name') if 'organism_name' in cols else None)


def _row_id(line: str, kind: str, fid_col: int | None) -> str | None:
    if kind == 'domtbl':
        if line.startswith('#'):
            return None
        m = _DOMTBL_ID.match(line)
        return m.group(2) if m else None
    if fid_col is None:
        return None
    cols = line.split('\t')
    return cols[fid_col] if fid_col < len(cols) else None


def _remap(line: str, kind: str, cols: tuple[int | None, int | None], fid: str, genome: str) -> str:
    if kind == 'domtbl':
        return _DOMTBL_ID.sub(lambda m: m.group(1) + fid, line, count=1)
    fid_col, org_col = cols
    parts = line.split('\t')
    if fid_col is not None and fid_col < len(parts):
        parts[fid_col] = fid
    if org_col is not None and org_col < len(parts):
        parts[org_col] = genome
    return '\t'.join(parts)


def _read_rows(path: Path, kind: str) -> tuple[str | None, dict[str, list[str]]]:
    """Returns (header, {feature id: lines}) for one results file."""
    if not path.exists() or path.stat().st_size == 0:
        return None, {}
    lines = path.read_text(encoding='utf-8', errors='replace').splitlines()
    rows: dict[str, list[str]] = {}
    if kind == 'domtbl':
        head = []
        for line in lines:
            if line.startswith('#'):
                if not rows:
                    head.append(line)
                continue
            fid = _row_id(line, kind, None)
            if fid:
                rows.setdefault(fid, []).append(line)
        return '\n'.join(head[:3]), rows
    header, data = lines[0], lines[1:]
    fid_col, _ = _columns(header)
    if fid_col is None:
        return header, {}
    for line in data:
        fid = _row_id(line, kind, fid_col)
        if fid:
            rows.setdefault(fid, []).append(line)
    return header, rows


# ---- before the annotation stage ----

def _proteins(faa: str) -> list[tuple[str, str, str, str]]:
    """Returns (feature id, hash, header, sequence) for each protein of the genome."""
    return [(fid, protein_hash(seq), header, seq) for fid, header, seq in read_fasta(faa) if fid and seq]


def split(db: str, faa: str, prefix: str, genome: str, tools: list[str], cfg: dict) -> dict[str, tuple[int, int]]:
    """Writes novel.faa and the cached proteins' rows per tool; returns {tool: (cached, to run)}."""
    proteins = _proteins(faa)
    hashes = sorted({h for _, h, _, _ in proteins})
    now = {}
    conn = _connect(db)
    try:
        for tool in tools:
            key = tool_key(tool, cfg)
            files = tool_files(tool, cfg)
            cached = _cached_hashes(conn, tool, key, files[0][0], hashes)
            out = Path(pc_dir(prefix, tool))
            if out.exists():
                shutil.rmtree(out)  # working files from an earlier attempt
            (out / 'cached').mkdir(parents=True)
            with open(out / 'novel.faa', 'w') as fh:
                for fid, h, header, seq in proteins:
                    if h not in cached:
                        fh.write(f'>{header}\n')
                        for i in range(0, len(seq), 60):
                            fh.write(seq[i:i + 60] + '\n')
            want = sorted(cached)
            for filename, kind in files:
                row = conn.execute('SELECT header FROM protein_cache_header WHERE tool = ? AND tool_key = ? AND filename = ?',
                                   (tool, key, filename)).fetchone()
                header = row[0] if row else None
                stored: dict[str, str] = {}
                for chunk in _chunks(want):
                    stored.update(conn.execute(
                        'SELECT aa_hash, lines FROM protein_cache WHERE tool = ? AND tool_key = ? AND filename = ? '
                        f'AND aa_hash IN ({",".join("?" * len(chunk))})', (tool, key, filename, *chunk)))
                if header is None and not any(stored.values()):
                    continue
                cols = _columns(header) if header and kind != 'domtbl' else (None, None)
                with open(out / 'cached' / filename, 'w') as fh:
                    if header:
                        fh.write(header + '\n')
                    for fid, h, _, _ in proteins:
                        for line in (stored.get(h) or '').splitlines():
                            fh.write(_remap(line, kind, cols, fid, genome) + '\n')
            (out / 'plan.json').write_text(json.dumps({
                'tool': tool, 'tool_key': key, 'cached': sum(1 for p in proteins if p[1] in cached),
                'novel': sum(1 for p in proteins if p[1] not in cached)}) + '\n')
            now[tool] = (sum(1 for p in proteins if p[1] in cached), sum(1 for p in proteins if p[1] not in cached))
    finally:
        conn.close()
    return now


# ---- after a tool's results are in ----

def store(db: str, faa: str, prefix: str, tool: str, cfg: dict) -> int:
    """Adds the genome's uncached proteins from the tool's final results; returns the count."""
    key = tool_key(tool, cfg)
    files = tool_files(tool, cfg)
    folder = Path(f'{prefix}{tool}')
    if not (folder / files[0][0]).exists():
        return 0
    proteins = _proteins(faa)
    first: dict[str, str] = {}
    for fid, h, _, _ in proteins:
        first.setdefault(h, fid)
    conn = _connect(db)
    try:
        new = sorted(set(first) - _cached_hashes(conn, tool, key, files[0][0], sorted(first)))
        if not new:
            return 0
        tables = {}
        for filename, kind in files:
            header, rows = _read_rows(folder / filename, kind)
            if header is None and kind == 'tsv':
                LOGGER.warning('protein cache: %s has no header; not caching %s', folder / filename, tool)
                return 0
            tables[filename] = (header, rows)
        # Skips caching when the results' ids are not this genome's proteins.
        ids = {fid for fid, _, _, _ in proteins}
        rows = tables[files[0][0]][1]
        strangers = sum(1 for fid in rows if fid not in ids)
        if rows and strangers > 0.01 * len(rows):
            LOGGER.warning('protein cache: %d of %d feature ids in %s are not proteins of %s; not caching %s',
                           strangers, len(rows), folder / files[0][0], faa, tool)
            return 0
        stamp = datetime.now(timezone.utc).isoformat()
        with conn:
            for filename, (header, _) in tables.items():
                if header:
                    conn.execute('INSERT OR IGNORE INTO protein_cache_header VALUES (?, ?, ?, ?)',
                                 (tool, key, filename, header))
            for filename, (_, rows) in tables.items():
                conn.executemany('INSERT OR IGNORE INTO protein_cache VALUES (?, ?, ?, ?, ?, ?)',
                                 [(h, tool, key, filename, '\n'.join(rows.get(first[h], [])), stamp) for h in new])
        return len(new)
    finally:
        conn.close()


# ---- inside the rule: merging tool rows with cached rows ----

# POSIX sh, since it runs inside each tool's container. Usage: merge.sh PRODUCED DEST PC_DIR FILENAME KIND
MERGE_SH = r'''#!/bin/sh
# Written by protein_cache.py. The tool's own rows (it ran on the proteins
# not in the cache), then the cached proteins' rows.
produced="$1"; dest="$2"; cached="$3/cached/$4"; kind="$5"
if [ -s "$produced" ]; then
  [ "$produced" = "$dest" ] || cp "$produced" "$dest" || exit 1
  if [ -s "$cached" ]; then
    if [ "$kind" = domtbl ]; then grep -v '^#' "$cached" >> "$dest" || true
    else tail -n +2 "$cached" >> "$dest" || exit 1; fi
  fi
elif [ -s "$cached" ]; then
  cp "$cached" "$dest" || exit 1
elif [ -e "$produced" ] || [ "$kind" = "tsv?" ]; then
  : > "$dest"
elif [ -e "$3/novel.faa" ] && [ ! -s "$3/novel.faa" ]; then
  : > "$dest"   # every protein was cached, and none has a row here
else
  echo "protein cache: $produced was not written" >&2; exit 1
fi
'''


def install_merge_script(folder: str) -> str:
    """Writes MERGE_SH into the run's output folder, where every container can read it."""
    path = Path(folder) / '.protein-cache' / 'merge.sh'
    try:
        if not path.exists() or path.read_text() != MERGE_SH:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix('.tmp')
            tmp.write_text(MERGE_SH)
            tmp.replace(path)
    except OSError as exc:
        LOGGER.warning('protein cache: could not write %s: %s', path, exc)
    return str(path)
