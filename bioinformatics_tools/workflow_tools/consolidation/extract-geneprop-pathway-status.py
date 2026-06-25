#!/usr/bin/env python3
"""extract-geneprop-pathway-status.py — derives
consolidated-only-geneprop-pathway-status.tsv directly from
merge-all-columns.py's output: for each gene, which GeneProp pathway(s)
it belongs to and that pathway's completion status (YES/PARTIAL/NO),
extracted into a clean per-pathway breakdown.

Pure extraction -- no cross-gene comparison logic here (that's a separate
question, "do operon-mates share a pathway", which belongs in scoring/ as
its own step once this extraction is confirmed).

GeneProp is genuinely multi-valued per gene -- a single gene can match
several independent pathways, each with its own status/step-completion.
The merge table carries this as "id1: val1; id2: val2" (the same
multi-hit convention used throughout this pipeline) for genes with more
than one pathway match, or a bare single value for genes with exactly
one. This script normalises both shapes into one consistent
"GenPropNNNN:STATUS" pairing per pathway, so a gene with one match and a
gene with several matches are equally easy to parse downstream.

GENEPROP_status values: YES (complete), PARTIAL (some but not all
required steps present), NO (present in the genome scan but doesn't meet
the completion threshold).

Output columns: identity columns, geneprop_pathway_count (how many
distinct pathways this gene matches, 0 if none), geneprop_pathway_status
("GenProp0263:PARTIAL;GenProp1162:PARTIAL" -- one ID:STATUS pair per
pathway, ";"-joined), plus the raw steps-matched/steps-total columns
passed through unchanged for full traceability.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from _shared import IDENTITY_COLUMNS

csv.field_size_limit(10_000_000)

_RAW_COLUMNS = [
    "GENEPROP_id", "GENEPROP_status", "GENEPROP_steps_matched",
    "GENEPROP_steps_matched_required", "GENEPROP_core_steps_required",
    "GENEPROP_steps_total_including_optional",
]


def parse_keyed(value: str) -> dict[str, str]:
    """Parse 'id1: val1; id2: val2' into {id1: val1, id2: val2}."""
    result: dict[str, str] = {}
    if not value:
        return result
    for part in value.split("; "):
        if ": " in part:
            k, v = part.split(": ", 1)
            result[k] = v
    return result


def normalize_pathway_status(geneprop_id: str, geneprop_status: str) -> list[tuple[str, str]]:
    """Returns [(pathway_id, status), ...] regardless of whether the
    source row was single-valued (bare) or multi-valued (keyed)."""
    if not geneprop_id:
        return []
    if ";" not in geneprop_id:
        return [(geneprop_id, geneprop_status)]
    ids = geneprop_id.split(";")
    status_map = parse_keyed(geneprop_status)
    return [(pid, status_map.get(pid, "")) for pid in ids]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[extract-geneprop-pathway-status] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = list(IDENTITY_COLUMNS) + [
        "geneprop_pathway_count", "geneprop_pathway_status",
    ] + [c.lower() for c in _RAW_COLUMNS]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    status_totals: dict[str, int] = {}
    n = 0
    with open(input_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            pairs = normalize_pathway_status(row.get("GENEPROP_id", ""), row.get("GENEPROP_status", ""))

            out_row = {col: row.get(col, "") for col in IDENTITY_COLUMNS}
            out_row["geneprop_pathway_count"] = str(len(pairs))
            out_row["geneprop_pathway_status"] = ";".join(f"{pid}:{status}" for pid, status in pairs)
            for raw_col in _RAW_COLUMNS:
                out_row[raw_col.lower()] = row.get(raw_col, "")

            writer.writerow(out_row)
            for _, status in pairs:
                status_totals[status or "(blank)"] = status_totals.get(status or "(blank)", 0) + 1
            n += 1

    total_pathway_mentions = sum(status_totals.values())
    print(f"[extract-geneprop-pathway-status] Wrote {n} genes → {output_path}")
    print(f"    {total_pathway_mentions} total pathway mentions across all genes:")
    for status, count in sorted(status_totals.items(), key=lambda kv: -kv[1]):
        print(f"        {status:10s} {count:6d}")


if __name__ == "__main__":
    main()
