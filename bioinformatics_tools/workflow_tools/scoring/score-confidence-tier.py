#!/usr/bin/env python3
"""score-confidence-tier.py — margie_sb phase11 (scoring), step 2:
combined confidence tier.

Reads score-hierarchy-tier.py's output (labeled-genes-hierarchy-tier.tsv)
and add-ec-consensus.py's output (labeled-genes-ec-consensus.tsv, from
phase10/labeling), joins on feature_id, and combines the two INDEPENDENT
confidence signals into one named confidence_tier -- how much trust the
winning tool's curation tier (hierarchy_tier_score) deserves, corroborated
or contradicted by what independent EC evidence says (ec_agreement_score).

Both READ-ONLY inputs -- this script never modifies either upstream
table, same "one core table, several focused derived views" pattern used
throughout this pipeline.

EC_AGREEMENT_SCORE mapping (see add-ec-consensus.py for the full
agreement-classification reasoning):
  full_consensus      -> +2  every independent tool that reported an EC
                              agrees -- strong corroboration.
  majority_consensus  -> +1  a shared core exists, but one tool reports
                              something extra unconfirmed.
  single_source       ->  0  real evidence, but nothing independent to
                              check it against -- neutral, not a penalty.
  no_evidence         ->  0  silence isn't evidence against the label --
                              most non-enzymatic genes legitimately have
                              no EC at all. Neutral, not a penalty.
  conflicting         -> -2  independent sources actively disagree with
                              EACH OTHER, regardless of how trusted the
                              winning tool's TEXT is. This is the one
                              status that overrides the combined score
                              entirely below.

CONFIDENCE_TIER (combined_score = hierarchy_tier_score + ec_agreement_score):
  flagged_for_review  -- ec_agreement_status == "conflicting", ALWAYS,
                          overriding combined_score. Independent evidence
                          disagreeing is worth a reviewer's attention no
                          matter how curated the winning tool is.
  high                -- combined_score >= 5 (e.g. tier4 hierarchy +
                          full_consensus, or tier4 + majority_consensus).
  moderate            -- combined_score in [2, 4].
  low                 -- combined_score < 2, OR hierarchy_tier_score was
                          already "no_qualifying_winner" (-1) to begin
                          with, regardless of any EC signal.

A single hierarchy-tier winner is not uniformly trustworthy -- among
genes won by the same tool, the underlying EC evidence can range from
full independent agreement to active conflict. This script turns that
split into an actionable per-gene flag instead of treating every winner
from one tool as equally reliable.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

EC_AGREEMENT_SCORE: dict[str, int] = {
    "full_consensus": 2,
    "majority_consensus": 1,
    "single_source": 0,
    "no_evidence": 0,
    "conflicting": -2,
}

_HIGH_THRESHOLD = 5
_MODERATE_THRESHOLD = 2

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def score_confidence_tier(hierarchy_tier_score: int, ec_agreement_status: str) -> tuple[int, str]:
    """Returns (combined_score, confidence_tier)."""
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

    ec_status_by_gene: dict[str, str] = {}
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

    tier_counts: dict[str, int] = {}
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
