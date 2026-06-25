#!/usr/bin/env python3
"""detect-columns.py — Stage 1 of margie_sb's consolidation pipeline.

A discovery/diagnostic step: scans output/<genome>/ and reports what's
actually there before merge-all-columns.py commits to a join -- which
tool directories exist, their real (raw) column names vs. how they'll be
normalised/prefixed, a real sample value per column, whether a gene can
ever have more than one row for that tool, and (when it can) exactly what
row-key strategy merge-all-columns.py will use to join those rows.

Pure inspection -- writes nothing back, never joins anything. Safe to run
any time as a sanity check, including against a genome that hasn't been
merged yet.

Output is JSON (detected-columns.json), not TSV -- this is inherently
nested per-tool metadata (a list of columns per tool, a key spec per
tool), which doesn't flatten cleanly into one TSV row shape.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from _shared import (
    ToolTable,
    discover_tool_tables,
    load_tool_table,
    resolve_key_spec,
)


def first_real_sample(table: ToolTable, col: str) -> str:
    for rows in table.rows_by_feature.values():
        for row in rows:
            val = row.get(col, "")
            if val:
                return val
    return ""


def max_rows_per_gene(table: ToolTable) -> int:
    if not table.rows_by_feature:
        return 0
    return max(len(rows) for rows in table.rows_by_feature.values())


def describe_tool(table: ToolTable) -> dict:
    max_rows = max_rows_per_gene(table)
    multi_row_capable = max_rows > 1
    spec = resolve_key_spec(table.tool_name)

    columns = []
    for col in table.tool_columns:
        columns.append({
            "name": col,
            "sample_value": first_real_sample(table, col),
        })

    return {
        "tool_name": table.tool_name,
        "source_file": str(table.source_path),
        "n_feature_ids": len(table.rows_by_feature),
        "n_raw_columns": len(table.raw_columns),
        "raw_columns": table.raw_columns,
        "max_rows_per_gene": max_rows,
        "multi_row_capable": multi_row_capable,
        "row_key_spec": (
            {"kind": spec.kind, "columns": list(spec.columns)}
            if multi_row_capable else None
        ),
        "columns": columns,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-root", required=True,
                   help="Genome output root, e.g. output/<genome>.")
    p.add_argument("--organism-name", default="")
    p.add_argument("--exclude-tool", action="append", default=[])
    p.add_argument("--output", required=True, help="Path for detected-columns.json")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    output_root = Path(args.input_root)
    if not output_root.is_dir():
        print(f"[detect-columns] ERROR: input root not a directory: {output_root}", file=sys.stderr)
        raise SystemExit(1)

    tool_table_pairs = discover_tool_tables(output_root, set(args.exclude_tool))
    if not tool_table_pairs:
        print("[detect-columns] ERROR: no per-tool tables found", file=sys.stderr)
        raise SystemExit(1)

    print(f"[detect-columns] {len(tool_table_pairs)} tool tables found under {output_root}:")
    tool_reports = []
    for tool_name, src in tool_table_pairs:
        table = load_tool_table(tool_name, src)
        report = describe_tool(table)
        tool_reports.append(report)
        flag = f"multi-row (max {report['max_rows_per_gene']})" if report["multi_row_capable"] else "single-row"
        print(f"    {tool_name:14s} ← {src}  ({report['n_feature_ids']} feature_ids, "
              f"{len(report['columns'])} cols, {flag})")

    result = {
        "organism_name": args.organism_name,
        "input_root": str(output_root),
        "n_tools": len(tool_reports),
        "tools": tool_reports,
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"[detect-columns] Wrote {len(tool_reports)} tool reports → {output_path}")
    print("[detect-columns] Done.")


if __name__ == "__main__":
    main()
