#!/usr/bin/env python3
"""score-confidence-final.py — margie_sb phase11 (scoring), final step:
combine C1/C2/C3/C4 into one confidence_score + confidence_flag.

Reads all four individual metric tables (score-c1-tool-coverage.py,
score-c2-operon-probability.py, score-c3-pathway-coherence.py,
score-c4-ec-agreement.py) -- ALL READ-ONLY -- joins on feature_id, and
writes ONE file carrying every component score, every formula string,
the combined confidence_score, and confidence_flag together. This is
the single reference table for "how confident are we in this gene's
canonical_label, and why" -- every number on it traces back to a
human-readable formula on the same row, no need to open any of the four
upstream files to audit a result.

    confidence_score = C1*0.40 + C2*0.05 + C3*0.15 + C4*0.40

WEIGHTS, and why they're not 0.4/0.1/0.3/0.2 (the original proposal):
  C1 and C4 (0.40 each) are the ALWAYS-APPLICABLE, DIRECT evidence
  signals -- tool coverage and cross-tool EC agreement -- so they
  dominate the weight equally (no more 2:1 C1-over-C4 imbalance).
  C2 and C3 (0.05 / 0.15) are CONTEXT-DEPENDENT bonuses that only mean
  anything when an operon exists -- their combined weight (0.20) is
  deliberately small so that sitting at their neutral default (0.5,
  for the many genes that are singletons) doesn't structurally cap an
  otherwise-perfectly-annotated gene below a high score. Verified
  worked example: a singleton with C1=1.0 and C4=1.0 (full EC
  consensus) now reaches confidence_score=0.90 exactly, where the
  original weights capped the same gene at 0.80 -- pure singleton
  status no longer prevents "highly confident."
  C3 outweighs C2 within that 0.20 (0.15 vs 0.05) because pathway
  coherence is a more direct functional cross-check than the operon's
  raw clustering probability.

confidence_flag: c4_ec_agreement_status == "conflicting" always forces
confidence_flag = "needs_review", regardless of the computed
confidence_score -- carried straight through from
score-c4-ec-agreement.py, not recomputed. A real, active disagreement
between independent EC sources deserves a reviewer's attention even when
three other components look fine; averaging it away would hide exactly
the thing most worth flagging.

This is explicitly a HEURISTIC, first-pass formula -- intended to be
validated against a broader set of genomes to check whether the
resulting confidence_score distribution actually tracks biological
plausibility, not assumed correct on the basis of a single genome.
Future signals (SEED-subsystem agreement for C3, synteny analysis,
LLM-based review including MEROPS/TCDB/DBCAN cross-checking) are
explicitly out of scope here -- to be layered on once this baseline is
validated.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

WEIGHTS = {"c1": 0.40, "c2": 0.05, "c3": 0.15, "c4": 0.40}

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def safe_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


# Named bands over the continuous confidence_score -- distinct from
# confidence_tier (hierarchy + EC agreement only). Kept as a pure function
# of the score; confidence_flag stays its own separate column rather than
# being folded into the tier name, so "what does the score say" and "is
# there an active conflict" remain independently visible.
def confidence_score_tier(score: float) -> str:
    if score >= 0.9:
        return "highest"
    if score >= 0.7:
        return "high"
    if score >= 0.5:
        return "moderate"
    if score >= 0.3:
        return "fair"
    return "low"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--c1-input", required=True, help="labeled-genes-c1-tool-coverage.tsv")
    parser.add_argument("--c2-input", required=True, help="labeled-genes-c2-operon-probability.tsv")
    parser.add_argument("--c3-input", required=True, help="labeled-genes-c3-pathway-coherence.tsv")
    parser.add_argument("--c4-input", required=True, help="labeled-genes-c4-ec-agreement.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    paths = {
        "c1": Path(args.c1_input), "c2": Path(args.c2_input),
        "c3": Path(args.c3_input), "c4": Path(args.c4_input),
    }
    for name, path in paths.items():
        if not path.is_file():
            print(f"[score-confidence-final] ERROR: {name} input not found: {path}", file=sys.stderr)
            raise SystemExit(1)

    def load(path: Path) -> dict[str, dict[str, str]]:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            return {row["feature_id"]: row for row in reader if row.get("feature_id")}

    c1_by_gene = load(paths["c1"])
    c2_by_gene = load(paths["c2"])
    c3_by_gene = load(paths["c3"])
    c4_by_gene = load(paths["c4"])

    out_columns = _IDENTITY_COLUMNS + [
        "c1_score", "c1_formula",
        "c2_score_from_operon_probability", "c2_formula",
        "c3_score", "c3_formula",
        "c4_score", "c4_formula", "c4_ec_agreement_status",
        "confidence_score", "confidence_score_formula", "confidence_score_tier", "confidence_flag",
        "hierarchy_tier_name", "hierarchy_tier_score", "combined_score", "confidence_tier",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    flag_counts: dict[str, int] = {}
    tier_counts: dict[str, int] = {}
    n = 0
    n_skipped = 0
    # Iterate over C4 (always populated for every gene, including non-coding
    # ones, since add-ec-consensus.py reads from labeled-genes.tsv directly).
    with open(paths["c4"], newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            c1_row = c1_by_gene.get(fid, {})
            c2_row = c2_by_gene.get(fid, {})
            c3_row = c3_by_gene.get(fid, {})

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c1_score"] = c1_row.get("c1_score", "")
            out_row["c1_formula"] = c1_row.get("c1_formula", "")
            out_row["c2_score_from_operon_probability"] = c2_row.get("c2_score_from_operon_probability", "")
            out_row["c2_formula"] = c2_row.get("c2_formula", "")
            out_row["c3_score"] = c3_row.get("c3_score", "")
            out_row["c3_formula"] = c3_row.get("c3_formula", "")
            out_row["c4_score"] = row.get("c4_score", "")
            out_row["c4_formula"] = row.get("c4_formula", "")
            out_row["c4_ec_agreement_status"] = row.get("c4_ec_agreement_status", "")
            out_row["hierarchy_tier_name"] = row.get("hierarchy_tier_name", "")
            out_row["hierarchy_tier_score"] = row.get("hierarchy_tier_score", "")
            out_row["combined_score"] = row.get("combined_score", "")
            out_row["confidence_tier"] = row.get("confidence_tier", "")

            c2_value = c2_row.get("c2_score_from_operon_probability", "")
            if not c2_value:
                # Non-coding feature -- C2 (and the rest of the formula) doesn't apply.
                out_row["confidence_score"] = ""
                out_row["confidence_score_formula"] = ""
                out_row["confidence_score_tier"] = "NOT_APPLICABLE_NON_CODING"
                out_row["confidence_flag"] = "NOT_APPLICABLE_NON_CODING"
                writer.writerow(out_row)
                n_skipped += 1
                continue

            c1 = safe_float(c1_row.get("c1_score", ""), 0.0)
            c2 = safe_float(c2_value, 0.5)
            c3 = safe_float(c3_row.get("c3_score", ""), 0.5)
            c4 = safe_float(row.get("c4_score", ""), 0.5)

            confidence_score = (c1 * WEIGHTS["c1"] + c2 * WEIGHTS["c2"]
                               + c3 * WEIGHTS["c3"] + c4 * WEIGHTS["c4"])
            ec_status = row.get("c4_ec_agreement_status", "no_evidence")
            confidence_flag = "needs_review" if ec_status == "conflicting" else "ok"

            out_row["confidence_score"] = f"{confidence_score:.4f}"
            out_row["confidence_score_formula"] = (
                f"C1({c1:.4f})*{WEIGHTS['c1']} + C2({c2:.4f})*{WEIGHTS['c2']} + "
                f"C3({c3:.4f})*{WEIGHTS['c3']} + C4({c4:.4f})*{WEIGHTS['c4']} = {confidence_score:.4f}"
            )
            out_row["confidence_score_tier"] = confidence_score_tier(confidence_score)
            out_row["confidence_flag"] = confidence_flag

            writer.writerow(out_row)
            flag_counts[confidence_flag] = flag_counts.get(confidence_flag, 0) + 1
            tier_counts[out_row["confidence_score_tier"]] = tier_counts.get(out_row["confidence_score_tier"], 0) + 1
            n += 1

    print(f"[score-confidence-final] Scored {n} protein-coding genes, "
          f"skipped {n_skipped} non-coding → {output_path}")
    print("  confidence_flag:")
    for flag, count in sorted(flag_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {flag:15s} {count:6d} ({100.0 * count / n:.1f}%)")
    print("  confidence_score_tier:")
    for tier, count in sorted(tier_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {tier:15s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
