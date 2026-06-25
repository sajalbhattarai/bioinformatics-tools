#!/usr/bin/env python3
"""extract-database-coverage.py — derives
consolidated-database-coverage-hits-boolean.tsv directly from
merge-all-columns.py's output: basic identity columns plus
{TOOL}_has_hit for every discovered tool, 1 if that tool found anything
for this gene, 0 otherwise.

Same tool categorization as extract-hit-counts.py (see _shared.py's
ID_BASED_TOOLS/SEGMENT_BASED_TOOLS/SINGLE_ROW_TOOLS) -- this is just that
script's per-tool count collapsed to a boolean (count >= 1 -> 1). Use
extract-hit-counts.py when the actual number of hits matters (e.g. "how
many PFAM domains"); use this one for a quick yes/no coverage matrix
across every tool at once (e.g. "which databases annotated this gene at
all").

Every listed tool gets its own {TOOL}_has_hit column regardless of
whether this genome happened to have any hits -- written as 0 rather
than omitted, so the schema is stable across genomes.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from _shared import ALL_HIT_COUNT_TOOLS, IDENTITY_COLUMNS, tool_hit_count

csv.field_size_limit(10_000_000)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[extract-database-coverage] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = list(IDENTITY_COLUMNS) + [f"{t}_has_hit" for t in ALL_HIT_COUNT_TOOLS]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        totals = {t: 0 for t in ALL_HIT_COUNT_TOOLS}
        with open(output_path, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                out_row = {col: row.get(col, "") for col in IDENTITY_COLUMNS}
                for tool in ALL_HIT_COUNT_TOOLS:
                    has_hit = 1 if tool_hit_count(row, tool) else 0
                    out_row[f"{tool}_has_hit"] = has_hit
                    totals[tool] += has_hit
                writer.writerow(out_row)
                n += 1

    print(f"[extract-database-coverage] Wrote {n} rows × {len(out_columns)} cols → {output_path}")
    for tool, count in totals.items():
        print(f"    {tool}_has_hit: {count} genes covered")


if __name__ == "__main__":
    main()
