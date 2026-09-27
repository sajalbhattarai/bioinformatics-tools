#!/usr/bin/env python3
"""score-confidence-final.py — scoring stage, final step: two-stage 0-1 confidence.

Joins the C1-C4 tables, operon info and localisation columns on feature_id and
writes every component with its reasoning, both stages, the tier and the review flag.
  preliminary = C1 * c4_score
  context     = C2 * C3 for operon members (C3 alone with --context-mode c3-only), else 0
  final       = clip(preliminary + context, 0, 1)
A parallel "hybrid" final uses the hybrid C3. Context only raises the score.
Tiers: > 0.9 highest | > 0.7 high | > 0.5 medium | > 0.3 fair | else low.
needs_review is set by EC conflict, final < 0.5, C2 < 0.5 in an operon, or an
ambiguous operon that gave no boost.
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_NEUTRAL = 0.5
# Minimum change from context that counts as "increases"/"decreases".
_CONTEXT_MATERIAL_THRESHOLD = 0.1
# Fraction of adjacent operon pairs involving an uncharacterized member at which
# a review note is added (never a score penalty).
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
    """Returns float(value), or default when it cannot be parsed."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _c1_reasoning(c1_row: dict) -> str:
    """Builds the C1 reasoning text from the C1 tool-coverage row."""
    tools = c1_row.get("c1_informative_tools", "")
    n = c1_row.get("c1_informative_tool_count", "?")
    total = c1_row.get("c1_total_tools_considered", "7")
    base = (f"{n}/{total} independent databases gave an informative hit"
            f" (of RAST/COG/PFAM/KEGG/EGGNOG/UNIPROT/TIGRFAM_CLUSTER)")
    return f"{base}: {tools}" if tools else f"{base} — none informative"


def _simplify_tmbed(raw: str) -> str:
    """Collapses a per-segment TMBED topology to its distinct types by regex.

    "0-31: signal_peptide; 32-1014: outside" -> "signal_peptide + outside"
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
    """Summarises MEROPS/TCDB/dbCAN classification codes for one gene.

    Format: "MEROPS:<family> | TCDB:<tc-number> | dbCAN:<CAZy-family> [EC ...]"; empty without hits.
    """
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
        # EC shown only when it contains a digit (skips placeholders like "-|-").
        segs.append(f"dbCAN:{dbcan}" + (f" [EC {ec}]" if any(c.isdigit() for c in ec) else ""))
    return " | ".join(segs)


def confidence_score_tier(final: float) -> str:
    """Maps the final 0-1 confidence to a tier: > 0.9 highest ... <= 0.3 low."""
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
    """Loads the component tables with csv, computes both confidence stages per gene and writes the final table."""
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
        """Reads a TSV with csv into a dict keyed by feature_id."""
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
            # RAST_feature_type (CDS / rna / prophage) lets later outputs colour non-coding rows.
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
            # Raw UniOP probability, so the 0.5 fallback is distinguishable.
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
            # ── specialized-database calls (protease/transporter/CAZyme) ──
            "specialized_db_hits",
            # ── alias columns read by the fingerprint and evidence steps;
            #    confidence_score == final confidence. ──
            "c1_score", "c1_informative_tools", "c1_formula",
            "c2_score_from_operon_probability", "c2_formula",
            "c3_score", "c3_signal_breakdown", "c3_formula",
            "c4_score", "c4_formula",
            "confidence_score", "confidence_score_formula", "confidence_score_tier", "confidence_flag",
            # read by make-final-annotated.py
            "hierarchy_tier_name",
        ]
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    flag_counts = {}
    tier_counts = {}
    n = 0
    n_skipped = 0
    # Iterates over C4, which has a row for every gene including non-coding ones.
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

            # ── Component scores, reasoning and alias columns ──
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
                # Non-coding feature: no confidence; empty confidence_score makes
                # build-gene-report.py skip it.
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
            c3 = safe_float(c3_raw, 0.0)         # conservation (0 = novel/no evidence)
            # Hybrid C3 (per-gene max of adjacency/co-member), scored in parallel.
            c3_hyb = safe_float(c3_row.get("c3_score_hybrid", ""), 0.0)
            conflict = safe_float(conflict_raw, 0.0)   # descriptor-consensus contradiction
            significance = safe_float(significance_raw, 0.0)  # enrichment: chance-above-random co-occurrence
            ambiguity = safe_float(ambiguity_raw, 0.0)  # m/n operon pairs blocked by a hypothetical
            try:
                _omc = int(float(op_row.get("operon_member_count", 0) or 0))
            except (TypeError, ValueError):
                _omc = 0
            in_operon = _omc >= 2
            # c4 is the EC-conflict clearance; its neutral value is 1.0, not 0.5.
            c4 = safe_float(c4_raw, 1.0)

            # ── Stage 1: gene's own evidence; C4 discounts C1 only on EC conflict. ──
            preliminary = c1 * c4
            # ── Stage 2: operon context; C3 (conservation) sets the boost and C2
            #    (operon probability) scales it. Non-operon genes get 0. ──
            if in_operon:
                if args.context_mode == "c3-only":
                    boost = max(0.0, c3)                     # conservation drives the boost directly
                else:
                    boost = max(0.0, c2) * max(0.0, c3)      # C3 (conservation) boosts; C2 gates
                # No penalty: descriptor-conflict signals misfire on lineage-specific
                # operons, so c3_descriptor_conflict is reported but not applied.
                penalty = 0.0
                context = boost
            else:
                boost = penalty = context = 0.0
            final = min(1.0, max(0.0, preliminary + context))

            # ── Hybrid final: same model with the hybrid C3. ──
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

            # ── needs_review triggers; a context increase is not a trigger. ──
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
            # Operon member with C2 < 0.5: the operon assignment itself is doubtful.
            if in_operon and c2 < 0.5:
                reasons.append(f"weak operon probability (C2={c2:.2f} < 0.5) -- this gene's "
                               f"operon assignment with its neighbour is doubtful")
            # Mostly uncharacterized operon that gave no boost: review note only.
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

            # alias: confidence_score == final confidence
            out_row["confidence_score"] = f"{final:.4f}"
            
            # Context formula text
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
