#!/usr/bin/env python3
"""score-confidence-final.py — margie_sb phase11 (scoring), final step:
combine C1/C2/C3/C4 into a TWO-STAGE 0-1 confidence score.

Reads all four individual metric tables (score-c1-tool-coverage.py,
score-c2-operon-probability.py, c3_score_organism.py,
score-c4-ec-agreement.py) -- ALL READ-ONLY -- joins on feature_id, and
writes ONE file carrying every component score, its human-readable
reasoning, both confidence stages, whether operon context helped, the
tier, and the review flag together. This is the single reference table
for "how confident are we in this gene's annotation, and why" -- every
number traces back to a reasoning string on the same row.

TWO-STAGE FORMULA (NON-FITTED, 0-1 SCALE):

  Stage 1 -- preliminary confidence from the gene's OWN evidence. C4 is now a
  graded EC-conflict CLEARANCE (c4_score = 1 - (m/5)*R, where R = conflicting
  EC tool-pairs / total EC tool-pairs and m = number of EC-reporting tools;
  1.0 = no conflict). It DISCOUNTS the C1 tool-coverage multiplicatively rather
  than being averaged co-equally with it, so silence/single-source (c4_score=1)
  never drags C1 down -- only a real EC conflict does:
    preliminary = C1 * c4_score           # C1 reduced only by real EC conflict

  Stage 2 -- fold in operon context (C3 gates C2 contribution). If C3 < 0.5
  (operon is NOT conserved across organisms), only C3 contributes (its penalty).
  If C3 >= 0.5 (operon IS conserved), both C2 and C3 contribute equally.
  This ensures low-OCC genes are penalized even if their operon probability
  is high, reflecting that species-specific operons don't prove universal function:
    if C3 < 0.5:
        context = (C3 - 0.5)            # Only use negative C3 signal (penalty)
    else:
        context = ((C2 - 0.5) + (C3 - 0.5)) / 2   # Both contribute (averaging)
    final   = clip(preliminary + context, 0.0, 1.0)

  does_context_improve_confidence? (±0.1 material threshold — trivial nudges
  are ignored so only operon contexts that actually move the needle count):
    final - preliminary >= +0.1 -> "increases"   (context corroborates)
    final - preliminary <= -0.1 -> "decreases"   (context undermines)
    otherwise                   -> "no effect"    (incl. singletons C2=C3=0.5)

Where each component is 0-1:
    C1 = Tool coverage: fraction of 7 independent tools with informative hits
    C2 = Operon probability: geometric-mean co-directionality (neutral 0.5)
    C3 = Operonic Context Confidence: three-level OCC (neutral 0.5)
    C4 = EC conflict clearance: 1 - (m/5)*R, the fraction of C1 that survives
         the EC-conflict check (neutral 1.0 = no conflict)

DESIGN RATIONALE:
  Two independent, non-fitted stages. Stage 1 is what the gene IS from its
  own homology + functional evidence; stage 2 asks whether its genomic
  (operon) context corroborates or undermines that. The C3 (operonic coherence)
  score GATES the C2 contribution: if the operon is poorly conserved across
  organisms (C3 < 0.5), the high C2 (operon probability) is ignored, reflecting
  that a species-specific operon doesn't prove universal function. Only when
  the operon is conserved (C3 >= 0.5) do both C2 and C3 boost confidence equally.

confidence_tier (on final, 0-1):
    > 0.9 highest | > 0.7 high | > 0.5 medium | > 0.3 fair | else low

needs_review = "yes" if ANY of:
  - EC conflict (c4_ec_agreement_status == "conflicting")
  - operon context DECREASES confidence by >= 0.1
  - low confidence (final < 0.5)
A context INCREASE is corroboration (good news) and is NOT flagged for review.

Row coloring (make-final-excel.py / ssh.py), priority red > blue > green:
  RED   : EC conflict OR operon context decreases confidence by >=0.1 (white text)
  BLUE  : final < 0.5
  GREEN : operon context increases confidence by >=0.1 (informational, not needs_review)
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# Two-stage non-fitted model (each component 0-1, final 0-1):
#   preliminary = C1 * c4_score        # c4_score = 1 - (m/5)*R EC-conflict clearance
#   context     = C2*C3 - C2*conflict  # C2 (operon prob) GATES; C3 (pure conservation)
#                                      # boosts, conflict (descriptor consensus) penalizes.
#                                      # Novelty (C3->0, no conflict) and non-operons -> 0.
#   final       = clip(preliminary + context, 0, 1)
# See CONFIDENCE_MODEL_DERIVATION.md + operon-context-neutral-derivation/ for why.
_NEUTRAL = 0.5
# Operon context must move final confidence by at least this much (either
# direction) to count as a material increase/decrease; smaller = "no effect".
_CONTEXT_MATERIAL_THRESHOLD = 0.1
# Operon-inference ambiguity (m/n = fraction of the operon's adjacent pairs
# blocked by a hypothetical/uncharacterized member). At or above this, and when
# operon context did NOT boost the call, we add a review COMMENT (never a score
# penalty — boost-only): the operon was too uncharacterized to corroborate.
_OPERON_AMBIGUITY_MIN = 0.5

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]

_OPERON_COLUMNS = ["operon_id", "operon_member_count", "operon_gene_position_in_operon"]

_LOCALIZATION_COLUMNS = [
    "SIGNALP6_prediction",
    "TMBED_topology",
    "PSORTB_localization",
    "PSORTB_score",
    "PSORTB_is_confident",
    "ENVELOPE_envelope_type",
    "ENVELOPE_inference_basis",
]

def safe_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _c1_reasoning(c1_row: dict) -> str:
    """Human-readable why-this-C1 string from the C1 tool-coverage row."""
    tools = c1_row.get("c1_informative_tools", "")
    n = c1_row.get("c1_informative_tool_count", "?")
    total = c1_row.get("c1_total_tools_considered", "7")
    base = (f"{n}/{total} independent databases gave an informative hit"
            f" (of RAST/COG/PFAM/KEGG/EGGNOG/UNIPROT/TIGRFAM_CLUSTER)")
    return f"{base}: {tools}" if tools else f"{base} — none informative"


def _simplify_tmbed(raw: str) -> str:
    """Collapse per-segment TMBED topology string to a compact set of unique types.

    Input: "0-31: signal_peptide; 32-1014: outside"
    Output: "signal_peptide + outside"
    """
    if not raw or raw == "inside":
        return raw
    types = []
    seen = set()
    for part in raw.split(";"):
        m = re.search(r":\s*(\S.*)", part.strip())
        if m:
            t = m.group(1).strip().rstrip(",")
            if t not in seen:
                seen.add(t)
                types.append(t)
    return " + ".join(types) if types else raw


def _specialized_db_hits(row: dict) -> str:
    """Compact, pipe-separated summary of the specialized-database calls
    (protease / transporter / CAZyme) for one gene, pulled from the consolidated
    merged table. Uses each database's clean classification code (the raw
    descriptions are verbose and inconsistent). Only databases with an actual
    hit are included; genes with none get an empty string.
    Format:  MEROPS:<family> | TCDB:<tc-number> | dbCAN:<CAZy-family> [EC ...]."""
    segs = []
    # MEROPS peptidase family (e.g. S85); fall back to the accession id.
    merops = row.get("MEROPS_family", "").strip() or row.get("MEROPS_id", "").strip()
    if merops:
        segs.append(f"MEROPS:{merops}")
    # TCDB transporter classification number (e.g. 2.A.1.2.20).
    tcdb = row.get("TCDB_id", "").strip()
    if tcdb:
        segs.append(f"TCDB:{tcdb}")
    # dbCAN CAZy family (e.g. GH73), with EC numbers when present.
    dbcan = row.get("DBCAN_id", "").strip()
    if dbcan:
        ec = row.get("DBCAN_ec_numbers", "").strip()
        # Only show EC if it carries a real number (skip placeholders like "-|-").
        segs.append(f"dbCAN:{dbcan}" + (f" [EC {ec}]" if any(c.isdigit() for c in ec) else ""))
    return " | ".join(segs)


# Confidence tier for the final 0-1 confidence score.
def confidence_score_tier(final: float) -> str:
    """Map final confidence (0-1) to tier: >0.9 highest ... <0.3 low."""
    if final > 0.9:
        return "highest"
    if final > 0.7:
        return "high"
    if final > 0.5:
        return "medium"
    if final > 0.3:
        return "fair"
    return "low"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--c1-input", required=True, help="labeled-genes-c1-tool-coverage.tsv")
    parser.add_argument("--c2-input", required=True, help="labeled-genes-c2-operon-probability.tsv")
    parser.add_argument("--c3-input", required=True, help="labeled-genes-c3-operonic-context-confidence.tsv")
    parser.add_argument("--c4-input", required=True, help="labeled-genes-c4-ec-agreement.tsv")
    parser.add_argument("--operon-input", required=True, help="labeled-genes-operon-info.tsv")
    parser.add_argument("--merged-input", required=True,
                        help="consolidated-merged-all-columns.tsv (localization columns)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--context-mode", choices=("c2-gated", "c3-only"),
                        default="c2-gated",
                        help="'c2-gated' (default, shipped): boost = C2*C3 (operon "
                             "probability gates the conservation boost). 'c3-only': "
                             "boost = C3 (conservation drives the boost directly; the "
                             "single-genome operon-probability gate is dropped).")
    args = parser.parse_args()

    paths = {
        "c1": Path(args.c1_input), "c2": Path(args.c2_input),
        "c3": Path(args.c3_input), "c4": Path(args.c4_input),
        "operon": Path(args.operon_input), "merged": Path(args.merged_input),
    }
    for name, path in paths.items():
        if not path.is_file():
            print(f"[score-confidence-final] ERROR: {name} input not found: {path}", file=sys.stderr)
            raise SystemExit(1)

    def load(path):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            return {row["feature_id"]: row for row in reader if row.get("feature_id")}

    c1_by_gene = load(paths["c1"])
    c2_by_gene = load(paths["c2"])
    c3_by_gene = load(paths["c3"])
    c4_by_gene = load(paths["c4"])
    operon_by_gene = load(paths["operon"])

    loc_by_gene = {}
    spec_by_gene = {}
    type_by_gene = {}
    with open(paths["merged"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if not fid or fid in loc_by_gene:
                continue
            loc_by_gene[fid] = {col: row.get(col, "") for col in _LOCALIZATION_COLUMNS}
            spec_by_gene[fid] = _specialized_db_hits(row)
            # RAST_feature_type: CDS / rna / prophage -- carried through so the FINAL
            # table, workbook and operon diagrams can distinguish coding features
            # from non-CDS (RNA) and colour/flag them accordingly.
            type_by_gene[fid] = (row.get("RAST_feature_type", "") or "").strip()

    out_columns = (
        _IDENTITY_COLUMNS
        + _OPERON_COLUMNS
        + _LOCALIZATION_COLUMNS
        + [
            # ── feature type (CDS / rna / prophage) for coloring + flagging ──
            "feature_type",
            # ── component scores with friendly names + reasoning ──────────
            "c1_score_database_coverage", "c1_score_reasoning",
            # raw UniOP operon probability, reported alongside C2 so the neutral
            # 0.5 fallback (singletons) is never mistaken for a real probability.
            "c2_uniop_probability_raw",
            "c2_score_operon_probability", "c2_score_reasoning",
            "c3_score_operon_context", "c3_score_operon_context_hybrid", "c3_reasoning",
            "c4_score_EC_agreement", "c4_reasoning", "c4_ec_agreement_status",
            # ── two-stage confidence ──────────────────────────────────────
            "preliminary_confidence_c1_c4",
            "final_confidence_operon_context",
            "final_confidence_operon_context_hybrid",
            "does_context_improve_confidence?",
            "confidence_tier",
            "confidence_tier_hybrid",
            "needs_review", "needs_review_reason",
            # ── specialized-database calls (protease/transporter/CAZyme),
            #    pipe-separated in one column, sits just before the copied
            #    alias columns below. ──
            "specialized_db_hits",
            # ── backward-compatible aliases for downstream fingerprint /
            #    evidence steps (add-gene-fingerprint, add-fingerprint-to-final,
            #    build-gene-report). confidence_score == final confidence. ──
            "c1_score", "c1_informative_tools", "c1_formula",
            "c2_score_from_operon_probability", "c2_formula",
            "c3_score", "c3_signal_breakdown", "c3_formula",
            "c4_score", "c4_formula",
            "confidence_score", "confidence_score_formula", "confidence_score_tier", "confidence_flag",
            # carried through for make-final-annotated's product_descriptor_hierarchy_tier_name column
            "hierarchy_tier_name",
        ]
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    flag_counts = {}
    tier_counts = {}
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
            op_row = operon_by_gene.get(fid, {})
            loc_row = loc_by_gene.get(fid, {})

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}

            # Operon presence
            out_row["operon_id"] = op_row.get("operon_id", "")
            out_row["operon_member_count"] = op_row.get("operon_member_count", "")
            out_row["operon_gene_position_in_operon"] = op_row.get("operon_gene_position_in_operon", "")

            # Localization summary
            for col in _LOCALIZATION_COLUMNS:
                val = loc_row.get(col, "")
                if col == "TMBED_topology":
                    val = _simplify_tmbed(val)
                out_row[col] = val

            # Specialized-database calls (MEROPS/TCDB/dbCAN), pipe-separated.
            out_row["specialized_db_hits"] = spec_by_gene.get(fid, "")

            # ── Component scores: friendly-named columns + reasoning, plus
            #    backward-compatible aliases read by downstream steps. ──
            c1_raw = c1_row.get("c1_score", "")
            c2_raw = c2_row.get("c2_score_from_operon_probability", "")
            c3_raw = c3_row.get("c3_score", "")
            conflict_raw = c3_row.get("c3_descriptor_conflict", "")
            significance_raw = c3_row.get("c3_operon_significance", "")
            ambiguity_raw = c3_row.get("c3_operon_ambiguity", "")
            c4_raw = row.get("c4_score", "")
            ec_status = row.get("c4_ec_agreement_status", "no_evidence")

            out_row["feature_type"] = type_by_gene.get(fid, "")
            out_row["c1_score_database_coverage"] = c1_raw
            out_row["c1_score_reasoning"] = _c1_reasoning(c1_row)
            out_row["c2_uniop_probability_raw"] = c2_row.get("c2_uniop_probability_raw", "")
            out_row["c2_score_operon_probability"] = c2_raw
            out_row["c2_score_reasoning"] = c2_row.get("c2_formula", "")
            out_row["c3_score_operon_context"] = c3_raw
            out_row["c3_reasoning"] = c3_row.get("c3_signal_breakdown", "") or c3_row.get("c3_formula", "")
            out_row["c4_score_EC_agreement"] = c4_raw
            out_row["c4_reasoning"] = row.get("c4_reasoning", "")
            out_row["c4_ec_agreement_status"] = ec_status
            out_row["hierarchy_tier_name"] = row.get("hierarchy_tier_name", "")

            # aliases
            out_row["c1_score"] = c1_raw
            out_row["c1_informative_tools"] = c1_row.get("c1_informative_tools", "")
            out_row["c1_formula"] = c1_row.get("c1_formula", "")
            out_row["c2_score_from_operon_probability"] = c2_raw
            out_row["c2_formula"] = c2_row.get("c2_formula", "")
            out_row["c3_score"] = c3_raw
            out_row["c3_signal_breakdown"] = c3_row.get("c3_signal_breakdown", "")
            out_row["c3_formula"] = c3_row.get("c3_formula", "")
            out_row["c4_score"] = c4_raw
            out_row["c4_formula"] = row.get("c4_formula", "")

            if not c2_raw:
                # Non-coding feature -- operon probability (and the two-stage
                # confidence) doesn't apply. confidence_score left EMPTY so
                # build-gene-report.py skips it as non-scored.
                for col in ("preliminary_confidence_c1_c4", "final_confidence_operon_context",
                            "final_confidence_operon_context_hybrid",
                            "c3_score_operon_context_hybrid"):
                    out_row[col] = ""
                out_row["does_context_improve_confidence?"] = "NOT_APPLICABLE_NON_CODING"
                out_row["confidence_tier"] = "NOT_APPLICABLE_NON_CODING"
                out_row["confidence_tier_hybrid"] = "NOT_APPLICABLE_NON_CODING"
                out_row["needs_review"] = "n/a"
                out_row["needs_review_reason"] = "non-coding — scoring not applicable"
                out_row["confidence_score"] = ""
                out_row["confidence_score_formula"] = ""
                out_row["confidence_score_tier"] = "NOT_APPLICABLE_NON_CODING"
                out_row["confidence_flag"] = "NOT_APPLICABLE_NON_CODING"
                writer.writerow(out_row)
                n_skipped += 1
                continue

            c1 = safe_float(c1_raw, 0.0)
            c2 = safe_float(c2_raw, 0.0)         # operon probability (0 = not/unknown operon)
            c3 = safe_float(c3_raw, 0.0)         # PURE conservation (0 = novel/no evidence)
            # hybrid C3 (per-gene max adjacency/co-member) -- scored in PARALLEL so the
            # output carries an adjacency final and a hybrid final side by side.
            c3_hyb = safe_float(c3_row.get("c3_score_hybrid", ""), 0.0)
            conflict = safe_float(conflict_raw, 0.0)   # descriptor-consensus contradiction
            significance = safe_float(significance_raw, 0.0)  # enrichment: chance-above-random co-occurrence
            ambiguity = safe_float(ambiguity_raw, 0.0)  # m/n operon pairs blocked by a hypothetical
            try:
                _omc = int(float(op_row.get("operon_member_count", 0) or 0))
            except (TypeError, ValueError):
                _omc = 0
            in_operon = _omc >= 2
            # c4 is the EC-conflict CLEARANCE (1 - (m/5)*R); neutral = 1.0
            # (no conflict = no penalty), NOT 0.5 like the operon-context terms.
            c4 = safe_float(c4_raw, 1.0)

            # ── Stage 1: preliminary from the gene's own evidence. C4 discounts
            #    C1 multiplicatively — only a real EC conflict (c4_score<1) lowers
            #    it; silence/single-source (c4_score=1) leaves base = C1. ──
            preliminary = c1 * c4
            # ── Stage 2: operon context, C2 as a MULTIPLICATIVE GATE. We score the
            #    functional call, so the EVIDENCE is cross-genome: C3 (conservation)
            #    corroborates it, `conflict` (descriptor consensus) contradicts it.
            #    C2 (operon probability) is only a GATE — it says whether there is a
            #    real operon to lend context, so it can DISCOUNT but never inflate.
            #    C3/conflict drive the magnitude; C2 scales it. Novelty (C3→0, no
            #    conflict) and non-operon genes are neutral. context ∈ [−1,+1]. ──
            if in_operon:
                if args.context_mode == "c3-only":
                    boost = max(0.0, c3)                     # conservation drives the boost directly
                else:
                    boost = max(0.0, c2) * max(0.0, c3)      # C3 (conservation) boosts; C2 gates
                # PENALTY DISABLED (boost-only). Every conflict formulation tried
                # (adjacency and co-membership) false-fires on lineage-specific
                # operons: it compares our descriptor to a mate's GLOBAL consensus,
                # which a novel module never matches, and the OCC includes the
                # scored genome (no leave-one-out). A trustworthy contradiction
                # signal needs leave-one-out + module-support gating + the operon
                # propensity prior -- an open problem, not this release. C3 stays
                # novelty-neutral, so nothing is penalized for merely being unseen.
                # c3_descriptor_conflict is still computed and emitted for study.
                penalty = 0.0
                context = boost
            else:
                boost = penalty = context = 0.0
            final = min(1.0, max(0.0, preliminary + context))

            # ── parallel HYBRID final: identical model, hybrid C3 in place of the
            #    adjacency C3.  Same context-mode gate. ──
            if in_operon:
                boost_hyb = (max(0.0, c3_hyb) if args.context_mode == "c3-only"
                             else max(0.0, c2) * max(0.0, c3_hyb))
            else:
                boost_hyb = 0.0
            final_hyb = min(1.0, max(0.0, preliminary + boost_hyb))

            delta = final - preliminary
            if delta >= _CONTEXT_MATERIAL_THRESHOLD:
                context_effect = "increases"
            elif delta <= -_CONTEXT_MATERIAL_THRESHOLD:
                context_effect = "decreases"
            else:
                context_effect = "no effect"

            tier = confidence_score_tier(final)

            # ── needs_review: EC conflict OR operon-driven DECREASE OR low.
            #    A context increase is corroboration (good news) → not flagged. ──
            ec_conflict = ec_status == "conflicting"
            context_drop = context_effect == "decreases"
            low_conf = final < _NEUTRAL
            reasons = []
            if ec_conflict:
                reasons.append("EC conflict — independent EC sources disagree")
            if context_drop:
                reasons.append(f"operon context lowers confidence by ≥{_CONTEXT_MATERIAL_THRESHOLD} "
                               f"({preliminary:.4f}→{final:.4f})")
            if low_conf:
                reasons.append(f"low confidence (final={final:.4f} < {_NEUTRAL})")
            # Weak operon probability: the gene is placed in an operon, but its
            # pairwise operon probability (C2, UniOP) with the neighbour is < 0.5,
            # so that operon assignment itself is doubtful -- flag for review.
            if in_operon and c2 < 0.5:
                reasons.append(f"weak operon probability (C2={c2:.2f} < 0.5) -- this gene's "
                               f"operon assignment with its neighbour is doubtful")
            # Operon member that operon context could NOT corroborate because the
            # operon is mostly uncharacterized: comment only (boost-only leaves the
            # score at the gene's own evidence; we never penalize for ignorance).
            operon_ambiguous = (in_operon and context_effect != "increases"
                                and ambiguity >= _OPERON_AMBIGUITY_MIN)
            if operon_ambiguous:
                reasons.append(
                    f"operon inference ambiguous — {ambiguity:.0%} of adjacent pairs "
                    f"involve an uncharacterized protein, so operon context could not "
                    f"corroborate this call (score left at own evidence, not penalized)")
            needs_review = "yes" if reasons else "no"

            out_row["preliminary_confidence_c1_c4"] = f"{preliminary:.4f}"
            out_row["final_confidence_operon_context"] = f"{final:.4f}"
            out_row["final_confidence_operon_context_hybrid"] = f"{final_hyb:.4f}"
            out_row["c3_score_operon_context_hybrid"] = f"{c3_hyb:.4f}"
            out_row["confidence_tier_hybrid"] = confidence_score_tier(final_hyb)
            out_row["does_context_improve_confidence?"] = context_effect
            out_row["confidence_tier"] = tier
            out_row["needs_review"] = needs_review
            out_row["needs_review_reason"] = "; ".join(reasons) if reasons else "no review triggers"

            # backward-compatible aliases (confidence_score == final confidence)
            out_row["confidence_score"] = f"{final:.4f}"
            
            # Build context formula string (C2-gated geometric mean)
            if in_operon:
                context_formula = (
                    f"context=C2·C3=({c2:.3f}·{c3:.3f})={boost:+.4f} "
                    f"[boost-only; C2 gates, C3 conservation. conflict={conflict:.3f}, "
                    f"sig={significance:.3f} computed but NOT applied -- penalty disabled]")
            else:
                context_formula = "context=0.0 [gene not in an operon]"
            
            out_row["confidence_score_formula"] = (
                f"preliminary=C1*c4_score=({c1:.4f})*({c4:.4f})={preliminary:.4f}; "
                f"{context_formula}; "
                f"final=clip(preliminary+context,0,1)={final:.4f} (range: 0-1)"
            )
            out_row["confidence_score_tier"] = tier
            out_row["confidence_flag"] = "needs_review" if needs_review == "yes" else "ok"

            writer.writerow(out_row)
            flag_counts[out_row["confidence_flag"]] = flag_counts.get(out_row["confidence_flag"], 0) + 1
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
            n += 1

    print(f"[score-confidence-final] Scored {n} protein-coding genes, "
          f"skipped {n_skipped} non-coding → {output_path}")
    print("  confidence_flag:")
    for flag, count in sorted(flag_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {flag:15s} {count:6d} ({100.0 * count / n:.1f}%)")
    print("  confidence_tier:")
    for tier_name, count in sorted(tier_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {tier_name:15s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
