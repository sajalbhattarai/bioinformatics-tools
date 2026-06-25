#!/usr/bin/env python3
"""score-hierarchy-tier.py — margie_sb phase11 (scoring), step 1: hierarchy
tier.

Reads labeled-genes.tsv (assign-canonical-label.py's output, READ-ONLY)
and maps each gene's label_source to a trust TIER -- the same grouping
already implicit in assign-canonical-label.py's _EVALUATORS priority
order, made explicit and numeric here so it can be combined with other
confidence signals later (EC agreement, etc.) without re-deriving it each
time.

This is intentionally the FIRST, standalone piece of the eventual
confidence-scoring step -- verify the tier bucketing alone reads as
sensible on real data before layering the EC-agreement signal on top in
a follow-up script. Lives in its own workflow_tools/scoring/ phase
folder (phase11), separate from labeling (phase10) -- labeling decides
WHICH tool's text wins; scoring assesses how much to trust a decision
already made, a different job that will keep growing as more signals
(EC agreement, eventually others) get folded in.

TIER RATIONALE (mirrors assign-canonical-label.py's own hierarchy
reasoning, just grouped into 5 buckets instead of 12 individual ranks):
  Tier 4 -- PGAP, TIGRFAM, HAMAP, NCBIFAM
            Prokaryote-built-from-day-one, equivalog/family-rule-level
            curation, --cut_tc/--cut_ga trusted cutoffs (or InterProScan's
            own internal cutoff) already enforced upstream.
  Tier 3 -- PIRSF, UNIPROT
            PIRSF: whole-protein (not domain-only) HMM classification.
            UniProt: Swiss-Prot's curation is gold-standard, but reached
            via a single best BLAST/DIAMOND alignment, gated at >=40%
            identity -- one step below a curated family HMM.
  Tier 2 -- PFAM, CDD
            Domain-level curation spanning all of life, not
            prokaryote-specific.
  Tier 1 -- KEGG, EGGNOG, COG
            Orthology/ortholog-group *inference*, broader/coarser calls
            than a direct family-membership HMM hit.
  Tier 0 -- RAST
            Always-available fallback; curation consistency varies a lot
            across RAST/SEED's subsystems.
  (no winner) -- label_source == "NONE": lower than every real tier --
            no qualifying evidence existed at all.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# tool_name -> (tier_score, tier_name). Same grouping logic as
# assign-canonical-label.py's _EVALUATORS order, just bucketed.
HIERARCHY_TIER: dict[str, tuple[int, str]] = {
    "PGAP": (4, "tier4_curated_prokaryote_family"),
    "TIGRFAM": (4, "tier4_curated_prokaryote_family"),
    "HAMAP": (4, "tier4_curated_prokaryote_family"),
    "NCBIFAM": (4, "tier4_curated_prokaryote_family"),
    "PIRSF": (3, "tier3_whole_protein_or_curated_best_hit"),
    "UNIPROT": (3, "tier3_whole_protein_or_curated_best_hit"),
    "PFAM": (2, "tier2_general_domain_curation"),
    "CDD": (2, "tier2_general_domain_curation"),
    "KEGG": (1, "tier1_orthology_inference"),
    "EGGNOG": (1, "tier1_orthology_inference"),
    "COG": (1, "tier1_orthology_inference"),
    "RAST": (0, "tier0_fallback"),
}
_NO_WINNER_SCORE, _NO_WINNER_NAME = -1, "no_qualifying_winner"

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def score_hierarchy_tier(label_source: str) -> tuple[int, str]:
    return HIERARCHY_TIER.get(label_source, (_NO_WINNER_SCORE, _NO_WINNER_NAME))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True,
                        help="assign-canonical-label.py's output TSV (labeled-genes.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    if not labeled_path.is_file():
        print(f"[score-hierarchy-tier] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = _IDENTITY_COLUMNS + ["hierarchy_tier_score", "hierarchy_tier_name"]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tier_counts: dict[str, int] = {}
    n = 0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            score, name = score_hierarchy_tier(row.get("label_source", ""))
            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["hierarchy_tier_score"] = str(score)
            out_row["hierarchy_tier_name"] = name
            writer.writerow(out_row)
            tier_counts[name] = tier_counts.get(name, 0) + 1
            n += 1

    print(f"[score-hierarchy-tier] Wrote {n} genes → {output_path}")
    for name, count in sorted(tier_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {name:40s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
