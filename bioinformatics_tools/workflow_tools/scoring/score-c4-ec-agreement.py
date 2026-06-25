#!/usr/bin/env python3
"""score-c4-ec-agreement.py — margie_sb phase11 (scoring), metric C4: EC
agreement.

Reads labeled-genes-ec-consensus.tsv (phase10, add-ec-consensus.py's
output, READ-ONLY) and converts ec_agreement_status into the C4 score.
Also joins in labeled-genes-confidence-tier.tsv (this same scoring phase,
score-confidence-tier.py's output, also READ-ONLY) purely for
PROVENANCE -- so a reviewer looking at this one file can see the FULL
step-by-step chain that led here (hierarchy tier -> ec agreement ->
combined score -> confidence tier), not just C4's own number in
isolation.

C4 mapping:
  full_consensus      -> 1.0   every independent tool that reported an EC
                                agrees -- strong corroboration.
  majority_consensus  -> 0.75  a shared core exists, one tool reports
                                something extra unconfirmed.
  single_source       -> 0.5   real evidence, nothing independent to
                                check it against -- neutral.
  no_evidence         -> 0.5   silence isn't conflict -- neutral, not a
                                penalty.
  conflicting         -> 0.0   independent sources actively disagree.

c4_score itself is deliberately NOT hierarchy-tier-aware -- folding
hierarchy trust into c4_score would double-count against C1, which
already measures tool-level support. hierarchy_tier_score/combined_score/
confidence_tier are carried here as separate, clearly-labeled PROVENANCE
columns (computed earlier by score-hierarchy-tier.py/
score-confidence-tier.py), not blended into c4_score.

confidence_flag: ec_agreement_status == "conflicting" sets
confidence_flag = "needs_review" regardless of any other signal --
carried here (not just inferable from the score) so this real
disagreement is visible from this file alone.

Output (labeled-genes-c4-ec-agreement.tsv): identity columns,
c4_ec_agreement_status, c4_score, c4_formula (literal text, e.g.
"ec_agreement_status=conflicting -> c4_score=0.0000"), confidence_flag,
then the joined provenance columns: hierarchy_tier_name,
hierarchy_tier_score, combined_score, confidence_tier, and
combined_score_formula (e.g. "hierarchy_tier_score(4) +
ec_agreement_score(-2) = combined_score(2) -> confidence_tier=
flagged_for_review").
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

EC_AGREEMENT_TO_C4: dict[str, float] = {
    "full_consensus": 1.0,
    "majority_consensus": 0.75,
    "single_source": 0.5,
    "no_evidence": 0.5,
    "conflicting": 0.0,
}

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ec-consensus-input", required=True, help="labeled-genes-ec-consensus.tsv")
    parser.add_argument("--confidence-tier-input", required=True,
                        help="score-confidence-tier.py's output TSV (labeled-genes-confidence-tier.tsv)"
                             " -- joined in for hierarchy/combined-score provenance only")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ec_path = Path(args.ec_consensus_input)
    tier_path = Path(args.confidence_tier_input)
    if not ec_path.is_file():
        print(f"[score-c4-ec-agreement] ERROR: input not found: {ec_path}", file=sys.stderr)
        raise SystemExit(1)
    if not tier_path.is_file():
        print(f"[score-c4-ec-agreement] ERROR: input not found: {tier_path}", file=sys.stderr)
        raise SystemExit(1)

    provenance_by_gene: dict[str, dict[str, str]] = {}
    with open(tier_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                provenance_by_gene[fid] = {
                    "hierarchy_tier_name": row.get("hierarchy_tier_name", ""),
                    "hierarchy_tier_score": row.get("hierarchy_tier_score", ""),
                    "ec_agreement_score": row.get("ec_agreement_score", ""),
                    "combined_score": row.get("combined_score", ""),
                    "confidence_tier": row.get("confidence_tier", ""),
                }

    out_columns = _IDENTITY_COLUMNS + [
        "c4_ec_agreement_status", "c4_score", "c4_formula", "confidence_flag",
        "hierarchy_tier_name", "hierarchy_tier_score", "combined_score",
        "confidence_tier", "combined_score_formula",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    flag_counts: dict[str, int] = {}
    n = 0
    with open(ec_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            status = row.get("ec_agreement_status", "no_evidence")
            c4 = EC_AGREEMENT_TO_C4.get(status, 0.5)
            flag = "needs_review" if status == "conflicting" else "ok"

            prov = provenance_by_gene.get(fid, {
                "hierarchy_tier_name": "", "hierarchy_tier_score": "",
                "ec_agreement_score": "", "combined_score": "", "confidence_tier": "",
            })

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c4_ec_agreement_status"] = status
            out_row["c4_score"] = f"{c4:.4f}"
            out_row["c4_formula"] = f"ec_agreement_status={status} -> c4_score={c4:.4f}"
            out_row["confidence_flag"] = flag
            out_row["hierarchy_tier_name"] = prov["hierarchy_tier_name"]
            out_row["hierarchy_tier_score"] = prov["hierarchy_tier_score"]
            out_row["combined_score"] = prov["combined_score"]
            out_row["confidence_tier"] = prov["confidence_tier"]
            if prov["hierarchy_tier_score"] and prov["ec_agreement_score"]:
                out_row["combined_score_formula"] = (
                    f"hierarchy_tier_score({prov['hierarchy_tier_score']}) + "
                    f"ec_agreement_score({prov['ec_agreement_score']}) = "
                    f"combined_score({prov['combined_score']}) -> confidence_tier={prov['confidence_tier']}"
                )
            else:
                out_row["combined_score_formula"] = ""

            writer.writerow(out_row)
            flag_counts[flag] = flag_counts.get(flag, 0) + 1
            n += 1

    print(f"[score-c4-ec-agreement] Wrote {n} genes → {output_path}")
    for flag, count in sorted(flag_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {flag:15s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
