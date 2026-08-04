#!/usr/bin/env python3
"""Admin utility scaffold for merging per-user SQLite databases into a new
shared versioned database.

This intentionally starts as a safe planner/validator and does not mutate any
inputs by default. It prints a deterministic merge plan and required checks.

Example:
  python scripts/merge_user_databases.py \
    --base-shared /depot/lindems/data/margie/sqlite/shared-v2.db \
    --output-shared /depot/lindems/data/margie/sqlite/shared-v3.db \
    --source-db /depot/lindems/data/margie/sqlite/bhattar3-shared-v2-v1.db \
    --source-db /depot/lindems/data/margie/sqlite/rraghun-shared-v2-v1.db \
    --dedupe-key fasta_hash
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path


SUPPORTED_DEDUPE_KEYS = ("fasta_hash", "genome_hash")


def _must_exist_file(path_text: str, label: str) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_file():
        raise SystemExit(f"{label} not found or not a file: {path}")
    return path


def _open_readonly(path: Path) -> sqlite3.Connection:
    # URI mode keeps this planner read-only for source/base checks.
    return sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)


def _list_tables(conn: sqlite3.Connection) -> set[str]:
    rows = conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    return {r[0] for r in rows}


def _validate_sources(base_shared: Path, source_dbs: list[Path], dedupe_key: str) -> dict:
    report: dict = {
        "base_shared": str(base_shared),
        "source_count": len(source_dbs),
        "dedupe_key": dedupe_key,
        "checks": [],
    }

    with _open_readonly(base_shared) as base_conn:
        base_tables = _list_tables(base_conn)
        report["checks"].append({
            "name": "base_tables_present",
            "ok": bool(base_tables),
            "detail": sorted(base_tables),
        })

    for source in source_dbs:
        with _open_readonly(source) as source_conn:
            tables = _list_tables(source_conn)
            report["checks"].append({
                "name": "source_tables_present",
                "source": str(source),
                "ok": bool(tables),
                "detail": sorted(tables),
            })

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Plan and validate a per-user DB merge into shared vN+1")
    parser.add_argument("--base-shared", required=True, help="Existing shared database version (read-only input)")
    parser.add_argument("--output-shared", required=True, help="Target shared database version to create")
    parser.add_argument("--source-db", action="append", default=[], help="Per-user source database (repeatable)")
    parser.add_argument("--dedupe-key", default="fasta_hash", help="Organism dedupe key (default: fasta_hash)")
    parser.add_argument("--execute", action="store_true", help="Execute merge (not yet implemented)")

    args = parser.parse_args()

    if args.dedupe_key not in SUPPORTED_DEDUPE_KEYS:
        raise SystemExit(
            f"Unsupported --dedupe-key '{args.dedupe_key}'. Supported: {', '.join(SUPPORTED_DEDUPE_KEYS)}"
        )

    if not args.source_db:
        raise SystemExit("At least one --source-db is required")

    base_shared = _must_exist_file(args.base_shared, "--base-shared")
    source_dbs = [_must_exist_file(v, "--source-db") for v in args.source_db]
    output_shared = Path(args.output_shared).expanduser()

    if output_shared.exists():
        raise SystemExit(f"--output-shared already exists: {output_shared}")

    report = _validate_sources(base_shared, source_dbs, args.dedupe_key)

    plan = {
        "mode": "execute" if args.execute else "plan",
        "base_shared": str(base_shared),
        "output_shared": str(output_shared),
        "sources": [str(p) for p in source_dbs],
        "dedupe_key": args.dedupe_key,
        "next_steps": [
            "copy base_shared to output_shared atomically",
            "attach each source db and upsert rows by dedupe_key",
            "record merge audit (inserted/duplicate/conflict counts)",
            "validate output_shared integrity and publish alias switch",
        ],
        "validation": report,
    }

    print(json.dumps(plan, indent=2, sort_keys=True))

    if args.execute:
        raise SystemExit("--execute requested, but merge execution is not implemented yet")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
