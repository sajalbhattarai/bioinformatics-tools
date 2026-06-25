#!/usr/bin/env python3
"""add-operon-info.py — margie_sb phase10 (labeling): operon context.

Reads labeled-genes.tsv (assign-canonical-label.py's output, READ-ONLY)
and consolidated-merged-all-columns.tsv (for the OPERON_* columns
produced by phase5's UniOP-based operon layer), joins on feature_id, and
writes its OWN standalone table: identity/decision columns for
cross-reference, plus operon identity and the operon-level probability.
Same "one core table, several focused derived views" pattern as
add-ec-consensus.py -- doesn't modify labeled-genes.tsv.

OPERON IDENTITY (operon_id): one of three states, kept distinct rather
than collapsed into a single blank/empty meaning:
  "operon_XXXX"              -- predicted member of a real multi-gene
                                 operon (>=2 genes).
  "NOT_IN_AN_OPERON"          -- UniOP ran on this gene and predicted it
                                 is a standalone singleton, not part of
                                 any multi-gene operon. A real, negative
                                 prediction.
  "NOT_APPLICABLE_NON_CODING" -- UniOP never saw this feature at all,
                                 since its input is a protein FAA and RNA
                                 genes (rRNA/tRNA) are never part of that
                                 input -- every gene with a blank
                                 OPERON_id in the merged table is an
                                 rna.* feature. Relabeled here from a
                                 bare blank specifically so it isn't
                                 mistaken for "predicted, not in an
                                 operon" -- this is "prediction doesn't
                                 apply to this feature type," a different
                                 fact entirely.

OPERON PROBABILITY: process_operon_raw_results.py (phase5) computes this
as the geometric mean of all N-1 adjacent-gene pairwise probabilities
within the operon (not an arithmetic mean -- geometric mean is the right
choice here because operon membership is a chain of independent-ish
pairwise calls multiplied together conceptually, so one weak link should
pull the summary down hard, the way a product would, rather than being
diluted by averaging). Stored upstream as a descriptive string, e.g.
"0.810748 (geometric_mean_of_adjacent_pairs, operon_genes=10,
n_adjacent_pairs=9)" -- parsed here into a clean numeric
operon_probability_geometric_mean column (the raw string is kept
alongside for full provenance, e.g. confirming exactly how many adjacent
pairs went into that number). Empty for singletons (no internal pairs to
average) and for NOT_APPLICABLE_NON_CODING genes.

Every gene in an operon shares the SAME operon_probability_geometric_mean
(it's an operon-level summary, not a per-gene score) -- this script does
not invent a per-gene operon confidence; that distinction matters for
anyone using this as a confidence signal later (e.g. scoring/phase11):
the number reflects how confident UniOP is that the WHOLE predicted gene
cluster is a real operon, not how confident it is about any one member.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_GEOMETRIC_MEAN_PATTERN = re.compile(r"^([\d.]+)\s*\(")

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def normalize_operon_id(raw: str) -> str:
    if raw == "":
        return "NOT_APPLICABLE_NON_CODING"
    return raw


def parse_operon_probability(raw: str) -> str:
    """"0.810748 (geometric_mean_of_adjacent_pairs, ...)" -> "0.810748".
    Empty string passes through unchanged (singletons / non-coding)."""
    if not raw:
        return ""
    match = _GEOMETRIC_MEAN_PATTERN.match(raw)
    return match.group(1) if match else ""


def main() -> None:
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
