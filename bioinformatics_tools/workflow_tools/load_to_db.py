"""
Load annotation tool output into a SQLite database.

Usage:
    python load_to_db.py gff  <input> <db_path> <source_tool> [--token <token_file>]
    python load_to_db.py csv  <input> <db_path> <table_name>  [--token <token_file>]
    python load_to_db.py tsv  <input> <db_path> <table_name>  [--token <token_file>]

Subcommands:
    gff  - Load GFF3 output (prodigal, etc.) into the `annotations` table
    csv  - Load any CSV with headers into a table named after the tool.
           Columns and types are inferred from the headers and data.
    tsv  - Load any TSV with headers into a table named after the tool.
           Same as csv but tab-delimited (e.g. COGclassifier output).

All subcommands write to the same .db file so all results live together.
"""

# ─────────────────────────── Pipeline version ────────────────────────── #
# Bump this string whenever any scoring/labeling/consolidation script
# changes its output in a way that makes existing DB rows stale.
# Existing rows with a different version are deleted and reloaded.
# is_already_processed() in workflow.py uses this to skip Stage 2
# entirely for organisms that are already complete at the current version.
PIPELINE_VERSION = "1.7.2027"
import argparse
import csv
import hashlib
import sqlite3
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

# consolidated-merged-all-columns.tsv carries columns like na_seq/aa_seq and
# concatenated multi-tool command_used strings well past Python's csv
# module's 131072-byte default field limit -- without raising it, loading
# that table raises _csv.Error: field larger than field limit.
csv.field_size_limit(10_000_000)


# Hardened connections, same as output_cache.py's _get_connection/
# _retry_operation: bare sqlite3.connect() doesn't retry on lock contention,
# and InterPro's 18-way concurrent load burst needs it.

def _get_connection(db_path: str, timeout: float = 120.0) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path, timeout=timeout)
    conn.execute(f"PRAGMA busy_timeout={int(timeout * 1000)}")
    return conn


def _retry_operation(func, max_retries: int = 12, initial_delay: float = 2.0):
    """Retry func() on transient lock/IO errors with exponential backoff.

    func is expected to open its own connection and close it on every call
    (including retries) rather than reuse one across attempts, since a
    connection that errored mid-transaction shouldn't be trusted afterward."""
    delay = initial_delay
    last_error: sqlite3.OperationalError | None = None
    for attempt in range(max_retries):
        try:
            return func()
        except sqlite3.OperationalError as e:
            last_error = e
            error_str = str(e).lower()
            if any(err in error_str for err in ("disk i/o error", "database is locked", "unable to open")):
                if attempt < max_retries - 1:
                    print(f"[load_to_db] database busy (attempt {attempt + 1}/{max_retries}): {e}. "
                          f"Retrying in {delay:.1f}s...", file=sys.stderr)
                    time.sleep(delay)
                    delay *= 2
            else:
                raise
    raise last_error


# ─────────────────────────── Provenance ─────────────────────────── #

CREATE_RUN_LOG_SQL = """
CREATE TABLE IF NOT EXISTS run_log (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL,
    input_hash TEXT NOT NULL,
    tool TEXT NOT NULL,
    input_path TEXT,
    row_count INTEGER,
    rules_completed INTEGER,
    status TEXT NOT NULL,
    loaded_at TEXT NOT NULL,
    fasta_hash TEXT,
    pipeline_version TEXT,
    UNIQUE(input_hash, tool)
);
"""

_RUN_LOG_NEW_COLS = [("fasta_hash", "TEXT"), ("pipeline_version", "TEXT"),
                     # organism_name is what the DATA tables are keyed by, but a
                     # genome's identity is its fasta_hash -- the same sequence can
                     # be submitted under a new display name. Recording the name
                     # alongside the hash is what lets a reload find and clear the
                     # rows it wrote under any EARLIER name for the same genome
                     # (see _organism_names_for_fasta).
                     ("organism_name", "TEXT")]


def _compute_file_hash(file_path: str) -> str:
    """Compute full SHA-256 hash of a file's contents."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


# Public alias used by workflow.py to compute the genome FASTA hash that
# keys is_already_processed().  Returns the full 64-char SHA-256 — distinct
# from output_cache.compute_file_hash which truncates to 16 chars for
# the tool output cache.
compute_fasta_hash = _compute_file_hash


def _ensure_run_log(conn: sqlite3.Connection) -> None:
    """Create run_log if absent; add new columns if the table predates them."""
    conn.execute(CREATE_RUN_LOG_SQL)
    existing = {row[1] for row in conn.execute("PRAGMA table_info(run_log)")}
    for col, col_type in _RUN_LOG_NEW_COLS:
        if col not in existing:
            conn.execute(f"ALTER TABLE run_log ADD COLUMN {col} {col_type}")


def _already_loaded(db_path: str, input_hash: str, tool: str,
                    fasta_hash: str | None = None) -> bool:
    """Return True if this tool's data is already current in the DB.

    When fasta_hash is provided, the check is (fasta_hash, PIPELINE_VERSION,
    tool) — version-aware.  Without it, falls back to (input_hash, tool)
    for backward-compat with non-genome tables (e.g. prodigal GFF).
    """
    if not Path(db_path).exists():
        return False

    def _check() -> bool:
        conn = _get_connection(db_path)
        try:
            _ensure_run_log(conn)
            if fasta_hash:
                row = conn.execute(
                    "SELECT id FROM run_log "
                    "WHERE fasta_hash = ? AND pipeline_version = ? AND tool = ? AND status = 'success'",
                    (fasta_hash, PIPELINE_VERSION, tool),
                ).fetchone()
            else:
                row = conn.execute(
                    "SELECT id FROM run_log WHERE input_hash = ? AND tool = ?",
                    (input_hash, tool),
                ).fetchone()
            return row is not None
        finally:
            conn.close()

    return _retry_operation(_check)


def is_already_processed(db_path: str, fasta_hash: str,
                         pipeline_version: str = PIPELINE_VERSION,
                         tool: str = "scoring_confidence_final") -> bool:
    """Return True when this genome FASTA was fully processed at pipeline_version.

    workflow.py calls this before launching Stage 2 so it can skip the
    entire Snakemake run for organisms already at the current version.
    """
    if not Path(db_path).exists():
        return False

    def _check() -> bool:
        conn = _get_connection(db_path)
        try:
            _ensure_run_log(conn)
            row = conn.execute(
                "SELECT id FROM run_log "
                "WHERE fasta_hash = ? AND pipeline_version = ? AND tool = ? AND status = 'success'",
                (fasta_hash, pipeline_version, tool),
            ).fetchone()
            return row is not None
        finally:
            conn.close()

    return _retry_operation(_check)


def _organism_names_for_fasta(db_path: str, fasta_hash: str) -> set[str]:
    """Every organism_name this FASTA has ever been loaded under.

    A genome's identity is its sequence, not its label. When the same FASTA is
    re-run under a new name, the data tables still hold the rows written under
    the OLD name -- and those rows are keyed only by organism_name, so a delete
    scoped to the new name cannot see them.
    """
    if not fasta_hash:
        return set()

    def _query() -> set[str]:
        conn = _get_connection(db_path)
        try:
            _ensure_run_log(conn)
            rows = conn.execute(
                "SELECT DISTINCT organism_name FROM run_log "
                "WHERE fasta_hash = ? AND organism_name IS NOT NULL AND organism_name != ''",
                (fasta_hash,),
            ).fetchall()
            return {r[0] for r in rows}
        finally:
            conn.close()

    try:
        return _retry_operation(_query)
    except sqlite3.OperationalError:
        return set()


def _delete_stale_organism_rows(db_path: str, table_name: str,
                                organism_name: str,
                                fasta_hash: str | None = None) -> int:
    """Delete existing rows for THIS GENOME from table_name before reloading.

    Called only when _already_loaded() returned False, so this is always safe:
    either the organism is new (0 rows deleted) or it exists at a stale/unversioned
    version (old rows cleared before fresh insert).  Handles pre-versioning data
    in margie.db that has no fasta_hash in run_log.

    Scoped by genome IDENTITY, not display name. Deleting only *organism_name*
    silently leaves a full duplicate copy of the same genome behind whenever it
    is re-run under a different label: the delete matches nothing, the insert
    adds a second set, and the DB ends up holding one genome twice under two
    names. That is reachable on every --force load (scoring does one each run),
    so the alias sweep below is what actually keeps the tables unique.
    """
    targets = {organism_name} | _organism_names_for_fasta(db_path, fasta_hash)
    targets = {t for t in targets if t}

    def _do_delete() -> int:
        conn = _get_connection(db_path)
        try:
            exists = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (table_name,),
            ).fetchone()
            if not exists:
                return 0
            cols = {row[1] for row in conn.execute(f'PRAGMA table_info("{table_name}")')}
            if "organism_name" not in cols:
                return 0
            placeholders = ",".join("?" * len(targets))
            cursor = conn.execute(
                f'DELETE FROM "{table_name}" WHERE organism_name IN ({placeholders})',
                tuple(sorted(targets)),
            )
            deleted = cursor.rowcount
            if deleted:
                conn.commit()
            return deleted
        finally:
            conn.close()

    if not targets:
        return 0
    if len(targets) > 1:
        print(f"Note: {table_name} — this FASTA was previously loaded as "
              f"{sorted(targets - {organism_name})}; clearing those rows too so the "
              f"same genome is not stored twice under different names.")
    return _retry_operation(_do_delete)


def _record_load(db_path: str, input_hash: str, tool: str,
                 input_path: str, row_count: int,
                 fasta_hash: str | None = None,
                 organism_name: str | None = None) -> None:
    """Record a successful annotation load in the run_log table."""
    def _insert() -> None:
        conn = _get_connection(db_path)
        try:
            _ensure_run_log(conn)
            conn.execute(
                "INSERT OR IGNORE INTO run_log "
                "(run_id, input_hash, tool, input_path, row_count, rules_completed, "
                " status, loaded_at, fasta_hash, pipeline_version, organism_name) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (str(uuid.uuid4()), input_hash, tool, input_path, row_count, 0, 'success',
                 datetime.now(timezone.utc).isoformat(),
                 fasta_hash, PIPELINE_VERSION if fasta_hash else None, organism_name),
            )
            conn.commit()
        finally:
            conn.close()

    _retry_operation(_insert)


# ──────────────────────────── GFF loader ──────────────────────────── #

CREATE_GFF_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS annotations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    seqid TEXT NOT NULL,
    source TEXT NOT NULL,
    type TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    score REAL,
    strand TEXT,
    phase TEXT,
    attributes TEXT,
    gene_id TEXT,
    partial TEXT,
    start_type TEXT,
    rbs_motif TEXT,
    gc_content REAL,
    confidence REAL
);
"""

INSERT_GFF_SQL = """
INSERT INTO annotations
    (seqid, source, type, start, end, score, strand, phase, attributes,
     gene_id, partial, start_type, rbs_motif, gc_content, confidence)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
"""


def parse_attributes(attr_string: str) -> dict:
    """Parse GFF3 attribute column (key=value;key=value) into a dict."""
    attrs = {}
    for pair in attr_string.strip().rstrip(";").split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def safe_float(value: str | None) -> float | None:
    if value is None or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_gff_to_db(gff_path: str, db_path: str, source_tool: str) -> int:
    """Parse a GFF3 file and insert rows into the annotations table."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue

            seqid, source, type_, start, end, score, strand, phase, attributes = cols
            attrs = parse_attributes(attributes)

            rows.append((
                seqid,
                source_tool,
                type_,
                int(start),
                int(end),
                safe_float(score),
                strand,
                phase,
                attributes,
                attrs.get("ID"),
                attrs.get("partial"),
                attrs.get("start_type"),
                attrs.get("rbs_motif"),
                safe_float(attrs.get("gc_cont")),
                safe_float(attrs.get("conf")),
            ))

    def _insert() -> None:
        conn = _get_connection(db_path)
        try:
            conn.execute(CREATE_GFF_TABLE_SQL)
            conn.executemany(INSERT_GFF_SQL, rows)
            conn.commit()
        finally:
            conn.close()

    _retry_operation(_insert)
    return len(rows)


# ──────────────────────────── CSV loader ──────────────────────────── #

def _infer_type(value: str) -> str:
    """Guess SQLite column type from a sample value."""
    # If value contains a comma, it's likely a list - treat as TEXT
    if ',' in value:
        return "TEXT"
    # If value contains special characters or ranges, treat as TEXT
    if any(c in value for c in ['-', '_', '/', ':', ';', ' ']):
        # Exception: negative numbers
        if value.startswith('-'):
            try:
                int(value)
                return "INTEGER"
            except ValueError:
                try:
                    float(value)
                    return "REAL"
                except ValueError:
                    return "TEXT"
        return "TEXT"
    try:
        int(value)
        return "INTEGER"
    except ValueError:
        pass
    try:
        float(value)
        return "REAL"
    except ValueError:
        pass
    return "TEXT"


def load_csv_to_db(csv_path: str, db_path: str, table_name: str,
                   delimiter: str = ",") -> int:
    """Load a delimited file with headers into a table named `table_name`.

    - Creates the table from headers if it doesn't exist.
    - Infers column types (INTEGER/REAL/TEXT) from sample rows.
    - Adds an autoincrement `id` primary key.
    """
    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        headers = [h.strip() for h in next(reader)]
        data_rows = list(reader)

    if not data_rows:
        return 0

    # Infer types from multiple sample rows (up to first 10) for better accuracy
    sample_size = min(10, len(data_rows))
    col_types = []
    for col_idx in range(len(headers)):
        # Check all sample values for this column
        inferred = "INTEGER"
        for row_idx in range(sample_size):
            if col_idx < len(data_rows[row_idx]):
                val = data_rows[row_idx][col_idx].strip()
                if val:  # Skip empty values
                    col_type = _infer_type(val)
                    # Use most permissive type seen: TEXT > REAL > INTEGER
                    if col_type == "TEXT":
                        inferred = "TEXT"
                        break
                    elif col_type == "REAL" and inferred == "INTEGER":
                        inferred = "REAL"
        col_types.append(inferred)

    quoted_headers = [f'"{h}"' for h in headers]
    col_defs = ",\n    ".join(
        f'"{h}" {t}' for h, t in zip(headers, col_types)
    )
    create_sql = f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        {col_defs}
    );
    """

    placeholders = ", ".join("?" for _ in headers)
    insert_sql = f"INSERT INTO {table_name} ({', '.join(quoted_headers)}) VALUES ({placeholders});"

    # Cast values to match inferred types. Falls back to a looser type (or
    # the raw string) if a later row doesn't fit the sampled type -- e.g.
    # quast's report.tsv mixes integer and float metric rows.
    def cast_row(row):
        result = []
        for val, typ in zip(row, col_types):
            val = val.strip()
            if not val:
                result.append(None)
            elif typ == "INTEGER":
                try:
                    result.append(int(val))
                except ValueError:
                    try:
                        result.append(float(val))
                    except ValueError:
                        result.append(val)
            elif typ == "REAL":
                try:
                    result.append(float(val))
                except ValueError:
                    result.append(val)
            else:
                result.append(val)
        return tuple(result)

    rows = [cast_row(r) for r in data_rows]

    def _insert() -> None:
        conn = _get_connection(db_path)
        try:
            conn.execute(create_sql)
            # Add any columns that are in the TSV but missing from the table
            # (handles schema evolution when new columns are added to output files)
            existing = {row[1] for row in conn.execute(f"PRAGMA table_info({table_name})")}
            for h, t in zip(headers, col_types):
                if h not in existing:
                    conn.execute(f'ALTER TABLE {table_name} ADD COLUMN "{h}" {t}')
            conn.executemany(insert_sql, rows)
            conn.commit()
        finally:
            conn.close()

    _retry_operation(_insert)
    return len(rows)


# ──────────────────────────── CLI ──────────────────────────── #

def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--token", help="Write a token file on success")
    p.add_argument("--fasta", metavar="PATH",
                   help="Genome FASTA path used to key the pipeline result cache")
    p.add_argument("--delete-organism", metavar="NAME",
                   help="Delete existing rows for this organism before loading "
                        "(used when pipeline_version changed)")
    p.add_argument("--force", action="store_true",
                   help="Reload even if this input was already loaded at the "
                        "current pipeline version (bypasses the version-aware "
                        "skip). Used for scoring, which must overwrite the DB on "
                        "every run because its OCC operon reference is a moving "
                        "target -- see rule load_scoring_to_db in margie_sb.smk.")


def main():
    parser = argparse.ArgumentParser(description="Load annotation output into SQLite")
    sub = parser.add_subparsers(dest="format", required=True)

    gff_p = sub.add_parser("gff", help="Load GFF3 file into the annotations table")
    gff_p.add_argument("input_file", help="Path to GFF3 file")
    gff_p.add_argument("db_path", help="Path to SQLite database")
    gff_p.add_argument("source_tool", help="Tool name (e.g. prodigal)")
    _add_common_args(gff_p)

    csv_p = sub.add_parser("csv", help="Load CSV file into a named table")
    csv_p.add_argument("input_file", help="Path to CSV file with headers")
    csv_p.add_argument("db_path", help="Path to SQLite database")
    csv_p.add_argument("table_name", help="Table name (e.g. pfam)")
    _add_common_args(csv_p)

    tsv_p = sub.add_parser("tsv", help="Load TSV file into a named table")
    tsv_p.add_argument("input_file", help="Path to TSV file with headers")
    tsv_p.add_argument("db_path", help="Path to SQLite database")
    tsv_p.add_argument("table_name", help="Table name (e.g. cog)")
    _add_common_args(tsv_p)

    args = parser.parse_args()

    if not Path(args.input_file).exists():
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    label = args.source_tool if args.format == "gff" else args.table_name
    fasta_hash = _compute_file_hash(args.fasta) if getattr(args, "fasta", None) else None
    organism_name = getattr(args, "delete_organism", None)

    # Version-aware skip: if this FASTA was already loaded at PIPELINE_VERSION, skip.
    # --force bypasses this so scoring reloads (overwrites) on every run -- its OCC
    # operon reference grows over time, so the same genome must be re-scored and the
    # DB refreshed each run (see rule load_scoring_to_db in margie_sb.smk).
    input_hash = _compute_file_hash(args.input_file)
    if not getattr(args, "force", False) and _already_loaded(args.db_path, input_hash, label, fasta_hash=fasta_hash):
        print(f"Skipped {label}: already at pipeline version {PIPELINE_VERSION}")
        if args.token:
            Path(args.token).parent.mkdir(parents=True, exist_ok=True)
            Path(args.token).write_text(f"0 rows loaded from {label} (already current)\n")
        return

    # Pre-delete any existing rows for this organism before loading fresh data.
    # Safe because _already_loaded() above would have returned True and exited
    # if current-version rows were already present.
    if organism_name:
        deleted = _delete_stale_organism_rows(args.db_path, label, organism_name,
                                              fasta_hash=fasta_hash)
        if deleted:
            print(f"Deleted {deleted} stale rows from {label} (pre-delete before reload)")

    if args.format == "gff":
        n = load_gff_to_db(args.input_file, args.db_path, args.source_tool)
    elif args.format == "tsv":
        n = load_csv_to_db(args.input_file, args.db_path, args.table_name,
                           delimiter="\t")
    else:
        n = load_csv_to_db(args.input_file, args.db_path, args.table_name)

    _record_load(args.db_path, input_hash, label, args.input_file, n,
                 fasta_hash=fasta_hash, organism_name=organism_name)
    print(f"Loaded {n} rows from {label} into {args.db_path}")

    if args.token:
        Path(args.token).parent.mkdir(parents=True, exist_ok=True)
        Path(args.token).write_text(f"{n} rows loaded from {label}\n")


if __name__ == "__main__":
    main()
