#!/usr/bin/env python3
"""score-c2-operon-probability.py — margie_sb phase11 (scoring), metric
C2: raw pairwise operon probability (per gene).

C2 is the gene's own raw UniOP pairwise operon probability -- its local
operon-bond strength. (The operon-level geometric-mean aggregation is C3's
job, Operon Context Confidence.)

A gene in the middle of an operon has TWO adjacency probabilities (to its
upstream and its downstream neighbour); a gene at either end of the operon has
one.  C2 = the MEAN of the gene's available raw pairwise probabilities.  These
per-pair values come from operon/operon_results.tsv
(OPERON_upstream_pairwise_probability / OPERON_downstream_pairwise_probability);
the operon-info file only carries the operon-level geometric mean, which C2
does not use.

C2 rules:
  - "NOT_APPLICABLE_NON_CODING" (RNA features)  -> no score at all ("").
  - "NOT_IN_AN_OPERON" (singleton)              -> neutral 0.5 (being a
    singleton isn't evidence against the label, just a different context).
  - in an operon, with >=1 raw pairwise prob    -> mean of those probs.
  - in an operon but no pairwise prob available -> neutral 0.5.

Output (labeled-genes-c2-operon-probability.tsv): identity columns,
c2_score_from_operon_probability, c2_operon_id, c2_formula.
"""
import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_NEUTRAL = 0.5

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]

_EMPTY_TOKENS = {"", "NA", "N/A", "None", "none", "nan", "NaN"}


def load_pairwise_probs(operon_results_path):
    """feature_id -> [raw pairwise probabilities] (upstream and/or downstream)."""
    probs = {}
    p = Path(operon_results_path)
    if not p.is_file():
        print(f"[score-c2-operon-probability] WARNING: operon_results not found: "
              f"{p} (operon members will fall back to neutral 0.5)", file=sys.stderr)
        return probs
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = (row.get("feature_id") or "").strip()
            if not fid:
                continue
            vals = []
            for col in ("OPERON_upstream_pairwise_probability",
                        "OPERON_downstream_pairwise_probability"):
                raw = (row.get(col) or "").strip()
                if raw in _EMPTY_TOKENS:
                    continue
                try:
                    vals.append(float(raw))
                except ValueError:
                    pass
            if vals:
                probs[fid] = vals
    return probs


def compute_c2(operon_id, pairwise_probs):
    """Returns (c2_text, raw_uniop_text, formula_text) for one gene.

    raw_uniop_text is the EXACT UniOP operon probability (the mean of the gene's
    raw pairwise probabilities) when one exists, and is BLANK when the gene has no
    UniOP probability at all -- a singleton, an operon member with no pairwise
    probability, or a non-coding feature.  C2 itself substitutes the neutral value
    0.5 in those cases, so reporting the raw value separately means a reader can
    always tell a genuine UniOP probability from the 0.5 neutral fallback."""
    if operon_id == "NOT_APPLICABLE_NON_CODING":
        return "", "", "non-coding feature, C2 not applicable"
    if operon_id == "NOT_IN_AN_OPERON":
        return (f"{_NEUTRAL:.4f}", "",
                f"singleton (NOT_IN_AN_OPERON): a gene not in an operon has no UniOP "
                f"pairwise probability, so C2 uses the neutral value {_NEUTRAL:.4f} as "
                f"a fallback (being a singleton is not evidence against the label, only "
                f"a context in which operon corroboration is unavailable). The raw "
                f"UniOP probability column is left blank.")
    if not pairwise_probs:
        return (f"{_NEUTRAL:.4f}", "",
                f"in {operon_id} but no raw pairwise probability available: C2 uses the "
                f"neutral value {_NEUTRAL:.4f} as a fallback; the raw UniOP probability "
                f"column is left blank.")
    value = sum(pairwise_probs) / len(pairwise_probs)
    shown = ", ".join(f"{v:.4f}" for v in pairwise_probs)
    detail = "raw pairwise prob" if len(pairwise_probs) == 1 else "mean of raw pairwise probs"
    return (f"{value:.4f}", f"{value:.4f}",
            f"{operon_id} {detail} [{shown}] = {value:.4f}: here C2 IS the raw UniOP "
            f"operon probability (reported also in the UniOP_operon_probability column).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--operon-input", required=True,
                        help="labeled-genes-operon-info.tsv (drives output rows + operon membership)")
    parser.add_argument("--operon-results", required=True,
                        help="operon/operon_results.tsv (UniOP per-pair probabilities)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    operon_path = Path(args.operon_input)
    if not operon_path.is_file():
        print(f"[score-c2-operon-probability] ERROR: input not found: {operon_path}", file=sys.stderr)
        raise SystemExit(1)

    pairwise = load_pairwise_probs(args.operon_results)

    out_columns = _IDENTITY_COLUMNS + ["c2_uniop_probability_raw",
                                       "c2_score_from_operon_probability",
                                       "c2_operon_id", "c2_formula"]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    state_counts = {}
    n = 0
    with open(operon_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            operon_id = row.get("operon_id", "")
            fid = (row.get("feature_id") or "").strip()
            c2, raw_uniop, formula = compute_c2(operon_id, pairwise.get(fid, []))
            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c2_uniop_probability_raw"] = raw_uniop
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
