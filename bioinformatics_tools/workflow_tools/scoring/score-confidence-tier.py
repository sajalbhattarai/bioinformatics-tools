#!/usr/bin/env python3
"""score-confidence-tier.py — scoring stage, step 2: combined confidence tier.

combined_score = hierarchy_tier_score + ec_agreement_score, where EC status maps
full_consensus +2, majority_consensus +1, single_source/no_evidence 0,
conflicting -2. Tiers: flagged_for_review whenever EC status is conflicting;
high at >= 5; moderate at 2-4; low below 2 or without a qualifying winner.
"""
import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

EC_AGREEMENT_SCORE = {
    "full_consensus": 2,
    "majority_consensus": 1,
    "single_source": 0,
    "no_evidence": 0,
    "conflicting": -2,
}

_HIGH_THRESHOLD = 5
_MODERATE_THRESHOLD = 2

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]


def score_confidence_tier(hierarchy_tier_score, ec_agreement_status):
    """Returns (combined_score, confidence_tier) from the hierarchy score and EC status."""
    ec_score = EC_AGREEMENT_SCORE.get(ec_agreement_status, 0)
    combined_score = hierarchy_tier_score + ec_score

    if ec_agreement_status == "conflicting":
        return combined_score, "flagged_for_review"
    if hierarchy_tier_score < 0:
        return combined_score, "low"
    if combined_score >= _HIGH_THRESHOLD:
        return combined_score, "high"
    if combined_score >= _MODERATE_THRESHOLD:
        return combined_score, "moderate"
    return combined_score, "low"


def main() -> None:
    """Joins hierarchy tiers to EC status with csv and writes the confidence tier table."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hierarchy-tier-input", required=True,
                        help="score-hierarchy-tier.py's output TSV (labeled-genes-hierarchy-tier.tsv)")
    parser.add_argument("--ec-consensus-input", required=True,
                        help="add-ec-consensus.py's output TSV (labeled-genes-ec-consensus.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    hierarchy_path = Path(args.hierarchy_tier_input)
    ec_path = Path(args.ec_consensus_input)
    if not hierarchy_path.is_file():
        print(f"[score-confidence-tier] ERROR: input not found: {hierarchy_path}", file=sys.stderr)
        raise SystemExit(1)
    if not ec_path.is_file():
        print(f"[score-confidence-tier] ERROR: input not found: {ec_path}", file=sys.stderr)
        raise SystemExit(1)

    ec_status_by_gene = {}
    with open(ec_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                ec_status_by_gene[fid] = row.get("ec_agreement_status", "no_evidence")

    out_columns = _IDENTITY_COLUMNS + [
        "hierarchy_tier_score", "hierarchy_tier_name",
        "ec_agreement_status", "ec_agreement_score",
        "combined_score", "confidence_tier",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tier_counts = {}
    n = 0
    with open(hierarchy_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            hierarchy_score = int(row.get("hierarchy_tier_score", "0") or "0")
            ec_status = ec_status_by_gene.get(fid, "no_evidence")
            combined_score, confidence_tier = score_confidence_tier(hierarchy_score, ec_status)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["hierarchy_tier_score"] = row.get("hierarchy_tier_score", "")
            out_row["hierarchy_tier_name"] = row.get("hierarchy_tier_name", "")
            out_row["ec_agreement_status"] = ec_status
            out_row["ec_agreement_score"] = str(EC_AGREEMENT_SCORE.get(ec_status, 0))
            out_row["combined_score"] = str(combined_score)
            out_row["confidence_tier"] = confidence_tier

            writer.writerow(out_row)
            tier_counts[confidence_tier] = tier_counts.get(confidence_tier, 0) + 1
            n += 1

    print(f"[score-confidence-tier] Wrote {n} genes → {output_path}")
    for tier, count in sorted(tier_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {tier:20s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
