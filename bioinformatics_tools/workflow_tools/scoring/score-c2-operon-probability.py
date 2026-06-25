#!/usr/bin/env python3
"""score-c2-operon-probability.py — margie_sb phase11 (scoring), metric
C2: operon probability.

Reads labeled-genes-operon-info.tsv (phase10, add-operon-info.py's
output, READ-ONLY) and converts the operon-level geometric-mean
probability into the C2 score.

C2 = operon_probability_geometric_mean, EXCEPT:
  - "NOT_IN_AN_OPERON" -> neutral 0.5, not 0. Being a singleton isn't
    evidence against the gene's label, just a different (extremely
    common) genomic context -- only a genuinely WEAK operon prediction
    should pull this down, not the mere absence of one.
  - "NOT_APPLICABLE_NON_CODING" (RNA features) -> no score at all.

Output (labeled-genes-c2-operon-probability.tsv): identity columns,
c2_score_from_operon_probability, c2_operon_id (kept for traceability -- which operon, or which
of the two non-multi-gene states, produced this score).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_NEUTRAL = 0.5

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def compute_c2(operon_id: str, geometric_mean: str) -> tuple[str, str]:
    """Returns (score, formula_text)."""
    if operon_id == "NOT_APPLICABLE_NON_CODING":
        return "", "non-coding feature, C2 not applicable"
    if operon_id == "NOT_IN_AN_OPERON":
        return f"{_NEUTRAL:.4f}", f"singleton (NOT_IN_AN_OPERON) -> c2_score={_NEUTRAL:.4f} (neutral)"
    if not geometric_mean:
        return f"{_NEUTRAL:.4f}", f"in {operon_id} but no geometric-mean probability available -> c2_score={_NEUTRAL:.4f} (neutral)"
    try:
        value = float(geometric_mean)
        return f"{value:.4f}", f"{operon_id} geometric-mean probability = {value:.4f} -> c2_score={value:.4f}"
    except ValueError:
        return f"{_NEUTRAL:.4f}", f"unparseable probability value -> c2_score={_NEUTRAL:.4f} (neutral)"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--operon-input", required=True, help="labeled-genes-operon-info.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    operon_path = Path(args.operon_input)
    if not operon_path.is_file():
        print(f"[score-c2-operon-probability] ERROR: input not found: {operon_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = _IDENTITY_COLUMNS + ["c2_score_from_operon_probability", "c2_operon_id", "c2_formula"]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    state_counts: dict[str, int] = {}
    n = 0
    with open(operon_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            operon_id = row.get("operon_id", "")
            c2, formula = compute_c2(operon_id, row.get("operon_probability_geometric_mean", ""))
            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c2_score_from_operon_probability"] = c2
            out_row["c2_operon_id"] = operon_id
            out_row["c2_formula"] = formula
            writer.writerow(out_row)
            state = "in_multi_gene_operon" if operon_id.startswith("operon_") else operon_id
            state_counts[state] = state_counts.get(state, 0) + 1
            n += 1

    print(f"[score-c2-operon-probability] Wrote {n} genes → {output_path}")
    for state, count in sorted(state_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {state:25s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
