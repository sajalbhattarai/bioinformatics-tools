#!/usr/bin/env python3
"""add-operon-info.py — labeling stage: operon context.

Joins labeled-genes.tsv with the merged table's OPERON_* columns (UniOP) and
writes identity columns plus operon id, size, position and probability.
operon_id is "operon_XXXX", "NOT_IN_AN_OPERON" (predicted singleton) or
"NOT_APPLICABLE_NON_CODING" (RNA features, which UniOP never sees).
operon_probability_geometric_mean is the operon-level geometric mean of
adjacent-pair probabilities, parsed from the raw string, which is kept too.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_GEOMETRIC_MEAN_PATTERN = re.compile(r"^([\d.]+)\s*\(")

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]


def normalize_operon_id(raw: str) -> str:
    """Returns the operon id, mapping an empty value to NOT_APPLICABLE_NON_CODING."""
    if raw == "":
        return "NOT_APPLICABLE_NON_CODING"
    return raw


def parse_operon_probability(raw: str) -> str:
    """Extracts the leading number from "0.810748 (geometric_mean_of_adjacent_pairs, ...)" by regex; empty stays empty."""
    if not raw:
        return ""
    match = _GEOMETRIC_MEAN_PATTERN.match(raw)
    return match.group(1) if match else ""


def main() -> None:
    """Joins labeled genes to OPERON_* columns with csv and writes the operon table."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True,
                        help="assign-canonical-label.py's output TSV (labeled-genes.tsv)")
    parser.add_argument("--merged-input", required=True,
                        help="merge-all-columns.py's output TSV (consolidated-merged-all-columns.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    merged_path = Path(args.merged_input)
    if not labeled_path.is_file():
        print(f"[add-operon-info] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)
    if not merged_path.is_file():
        print(f"[add-operon-info] ERROR: input not found: {merged_path}", file=sys.stderr)
        raise SystemExit(1)

    operon_by_gene: dict[str, dict[str, str]] = {}
    with open(merged_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                operon_by_gene[fid] = {
                    "operon_id": normalize_operon_id(row.get("OPERON_id", "")),
                    "operon_member_count": row.get("OPERON_member_count", ""),
                    "operon_gene_position_in_operon": row.get("OPERON_gene_position_in_operon", ""),
                    "operon_probability_geometric_mean": parse_operon_probability(
                        row.get("OPERON_probability", "")),
                    "operon_probability_raw": row.get("OPERON_probability", ""),
                }

    out_columns = _IDENTITY_COLUMNS + [
        "operon_id", "operon_member_count", "operon_gene_position_in_operon",
        "operon_probability_geometric_mean", "operon_probability_raw",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    id_state_counts: dict[str, int] = {}
    n = 0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            operon_info = operon_by_gene.get(fid, {
                "operon_id": "NOT_APPLICABLE_NON_CODING",
                "operon_member_count": "",
                "operon_gene_position_in_operon": "",
                "operon_probability_geometric_mean": "",
                "operon_probability_raw": "",
            })

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row.update(operon_info)

            writer.writerow(out_row)
            state = "in_multi_gene_operon" if operon_info["operon_id"].startswith("operon_") \
                else operon_info["operon_id"]
            id_state_counts[state] = id_state_counts.get(state, 0) + 1
            n += 1

    print(f"[add-operon-info] Wrote {n} genes → {output_path}")
    for state, count in sorted(id_state_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {state:25s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
