#!/usr/bin/env python3
"""extract-hit-counts.py — derives consolidated-hits-counts-only.tsv
directly from merge-all-columns.py's output: basic identity columns plus
{TOOL}_hit_count for every discovered tool, recording how many distinct
hits that tool found for each gene.

"Hit count" means a different count depending on the tool's shape (see
_shared.py's ID_BASED_TOOLS/SEGMENT_BASED_TOOLS/SINGLE_ROW_TOOLS for the
full categorization, shared with extract-database-coverage.py):
  - Accession-ID tools (COG, KEGG, eggNOG, PFAM, TIGRFAM, PGAP, MEROPS,
    TCDB, dbCAN, UniProt, GeneProp, the shared InterPro reference, and
    every interpro_<db> table): count of the ";"-joined accessions in
    their bare {PREFIX}_id column, 0 if empty.
  - Segment/region tools with no accession (Phobius, TMbed, DeepSig):
    counted from their own representative raw column -- a segment/region
    count, not a database-hit count, but the closest equivalent these
    tools have.
  - Always-single-row tools (PSORTb, SignalP4, SignalP6, Operon): 1 if
    that tool made a positive call for the gene, 0 otherwise.
  - rasttk is excluded -- every gene has exactly one rasttk entry by
    definition, not a "hit" against an external reference.

Every listed tool gets its own {TOOL}_hit_count column regardless of
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
        print(f"[extract-hit-counts] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = list(IDENTITY_COLUMNS) + [f"{t}_hit_count" for t in ALL_HIT_COUNT_TOOLS]

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
                    c = tool_hit_count(row, tool)
                    out_row[f"{tool}_hit_count"] = c
                    totals[tool] += 1 if c else 0
                writer.writerow(out_row)
                n += 1

    print(f"[extract-hit-counts] Wrote {n} rows × {len(out_columns)} cols → {output_path}")
    for tool, count in totals.items():
        print(f"    {tool}_hit_count: {count} genes with >=1 hit")


if __name__ == "__main__":
    main()
