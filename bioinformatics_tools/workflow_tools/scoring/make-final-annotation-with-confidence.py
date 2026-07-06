#!/usr/bin/env python3
"""Build FINAL_ANNOTATION_WITH_CONFIDENCE.tsv for user download.

This file is a scoring-phase, user-facing table that presents gene identity,
descriptor provenance, localization/topology context, and C1-C4 confidence
explanations in one ordered sheet.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import c3_occ  # noqa: E402
import c3_lib as L  # noqa: E402
import c3_score_organism as C3S  # noqa: E402

_PAIRWISE_EMPTY = {"", "NA", "N/A", "None", "none", "nan", "NaN"}
_CONTEXT_THRESHOLD = 0.1

BASE_COLUMNS = [
    "organism_name",
    "domain",
    "envelope",
    "gene_id",
    "RAST_feature_id",
    "FEATURE_TYPE",
    "RAST_strand",
    "RAST_start",
    "RAST_end",
    "RAST_na_sequence",
    "RAST_aa_sequence",
    "IS_IN_OPERON?",
    "best_consensus_product_descriptor",
    "best_consensus_product_descriptor_source",
    "best_consensus_product_descriptor_source_hierarchy_order",
    "best_consensus_product_descriptor_source_audit_trail",
    "gene_description_fingerprint",
    "specialized_database_hits",
    "localization_and_topology_hits",
    "ENVELOPE",
    "ENVELOPE_reason",
    "C1_score_database_coverage",
    "C1_reasoning",
    "UniOP_OPERON_id",
    "UniOP_operon_probability",
    "C2_score_pairwise_genes_UniOP_probability",
    "C2_score_formula",
    "C2_score_reasoning",
    "C3_score_operon_context",
    "C3_score_operon_context_hybrid",
    "C3_score_formula",
    "C3_score_reasoning",
    "C4_score_EC_conflict",
    "C4_score_formula",
    "C4_score_reasoning",
    "EC_EVIDENCE_STATUS",
    "PRELIMINARY_confidence_C1_C4",
    "PRELIMINARY_confidence_C1_C4_formula",
    "PRELIMINARY_confidence_C1_C4_reasoning",
    "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT",
    "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_hybrid",
    "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_formula",
    "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_reasoning",
    "DOES_OPERON_CONTEXT_IMPROVE_CONFIDENCE?",
    "CONFIDENCE_TIER",
    "CONFIDENCE_TIER_hybrid",
    "NEEDS_REVIEW?",
    "NEEDS_REVIEW_REASON",
    "BEST_PRODUCT_DESCRIPTOR(copied_here_for_convenience)",
]


def _load_tsv(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = (row.get("feature_id") or "").strip()
            if fid and fid not in out:
                out[fid] = row
    return out


def _excel_col_name(i: int) -> str:
    n = i + 1
    out = []
    while n > 0:
        n, r = divmod(n - 1, 26)
        out.append(chr(ord("A") + r))
    return "".join(reversed(out))


def _prefixed_columns(cols: list[str]) -> list[str]:
    return [f"Column-{_excel_col_name(i)}: {name}" for i, name in enumerate(cols)]


def _safe_float(v: str, default: float) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _desc_fp(desc: str) -> str:
    if not desc:
        return ""
    return hashlib.sha256(desc.encode("utf-8")).hexdigest()[:16]


def _pairwise_map(operon_results: Path) -> dict[str, list[tuple[str, float]]]:
    out: dict[str, list[tuple[str, float]]] = {}
    with open(operon_results, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = (row.get("feature_id") or "").strip()
            if not fid:
                continue
            vals: list[tuple[str, float]] = []
            for gcol, pcol in (
                ("OPERON_upstream_gene_id", "OPERON_upstream_pairwise_probability"),
                ("OPERON_downstream_gene_id", "OPERON_downstream_pairwise_probability"),
            ):
                other = (row.get(gcol) or "").strip()
                rawp = (row.get(pcol) or "").strip()
                if not other or rawp in _PAIRWISE_EMPTY:
                    continue
                try:
                    vals.append((other, float(rawp)))
                except ValueError:
                    continue
            out[fid] = vals
    return out


def _localization_summary(cf: dict[str, str]) -> str:
    segs = []
    if cf.get("SIGNALP6_prediction"):
        segs.append(f"SIGNALP6={cf.get('SIGNALP6_prediction')}")
    if cf.get("TMBED_topology"):
        segs.append(f"TMBED={cf.get('TMBED_topology')}")
    if cf.get("PSORTB_localization"):
        segs.append(f"PSORTB={cf.get('PSORTB_localization')}")
    return " | ".join(segs)


def _hierarchy_order(hrow: dict[str, str]) -> str:
    s = (hrow.get("hierarchy_tier_score") or "").strip()
    n = (hrow.get("hierarchy_tier_name") or "").strip()
    if s and n:
        return f"{n} (score={s})"
    return n or s


def _c2_reasoning(
    fid: str,
    operon_id: str,
    member_count: str,
    pairwise: dict[str, list[tuple[str, float]]],
    descriptor_by_fid: dict[str, str],
    c2_score_text: str,
) -> str:
    if operon_id == "NOT_APPLICABLE_NON_CODING":
        return "non-coding feature: C2 not applicable (UniOP_operon_probability column left blank)"
    if operon_id == "NOT_IN_AN_OPERON":
        return (
            "singleton (not in a multi-gene operon): this gene has no UniOP pairwise "
            "probability, so no raw value is shown in the UniOP_operon_probability column. "
            "C2 falls back to the neutral value 0.5000 -- the minimum C2 can take -- because "
            "being a singleton is not evidence against the label, only a context in which "
            "operon corroboration is unavailable; 0.5 is neutral (neither corroborating nor "
            "contradicting) so it neither rewards nor penalises the gene."
        )
    pairs = pairwise.get(fid, [])
    if not pairs:
        return (
            f"{operon_id} (total number of genes in this operon = {member_count}): "
            "no usable UniOP pairwise probability for this gene's adjacent operon edge(s), so "
            "the UniOP_operon_probability column is left blank. C2 falls back to the neutral "
            "value 0.5000 -- the minimum C2 can take -- so a missing probability neither "
            "rewards nor penalises the gene."
        )
    own_desc = descriptor_by_fid.get(fid, "")
    details = []
    for other, p in pairs:
        other_desc = descriptor_by_fid.get(other, "")
        details.append(f"({own_desc}) <-> ({other_desc}) = {p:.4f}")
    return (
        f"{operon_id} (total number of genes in this operon = {member_count}): "
        + "; ".join(details)
        + f"; mean pairwise probability = {c2_score_text or '0.5000'}. Here C2 IS this raw "
        "UniOP operon probability (reported unchanged in the UniOP_operon_probability column); "
        "the 0.5 floor is never invoked because UniOP within-operon probabilities are always "
        ">= 0.5."
    )


def _operon_members(
    labeled_by_fid: dict[str, dict[str, str]],
    operon_by_fid: dict[str, dict[str, str]],
) -> dict[str, list[tuple[str, int]]]:
    by_operon: dict[str, list[tuple[str, int]]] = {}
    for fid, op in operon_by_fid.items():
        operon_id = (op.get("operon_id") or "").strip()
        if not operon_id.startswith("operon_"):
            continue
        lr = labeled_by_fid.get(fid, {})
        start = 0
        try:
            start = int((lr.get("gene_start") or "0").strip())
        except ValueError:
            start = 0
        by_operon.setdefault(operon_id, []).append((fid, start))
    for operon_id in by_operon:
        by_operon[operon_id].sort(key=lambda x: x[1])
    return by_operon


def _c3_pair_details_for_operon(
    operon_members: list[tuple[str, int]],
    descriptor_by_fid: dict[str, str],
    clean_by_fid: dict[str, str],
    uninf_by_fid: dict[str, bool],
    pair_prob_by_pair: dict[frozenset[str], float],
    ref: dict,
) -> str:
    parts = []
    for (a, _), (b, _) in zip(operon_members[:-1], operon_members[1:]):
        da = descriptor_by_fid.get(a, "")
        db = descriptor_by_fid.get(b, "")
        p = pair_prob_by_pair.get(frozenset((a, b)))
        if p is None:
            parts.append(f"({da})<->({db}): UniOP=NA (unlinked)")
            continue
        if uninf_by_fid.get(a, True) or uninf_by_fid.get(b, True):
            parts.append(f"({da})<->({db}): hypothetical member -> EXCLUDED from C3")
            continue
        ca = clean_by_fid.get(a, "")
        cb = clean_by_fid.get(b, "")
        rho_lb = c3_occ.rho_adj(ca, cb, ref)
        rho_mean = C3S.rho_adj_mean(ca, cb, ref)
        if rho_mean <= 0.0 and rho_lb <= 0.0:
            parts.append(f"({da})<->({db}): no OCC evidence (novel) -> EXCLUDED from C3")
            continue
        # C3 is PURE conservation now: term = rho_adj (no P_ab). Operon probability
        # C2 gates in the final score, not here; UniOP pairwise shown for reference.
        parts.append(
            f"({da})<->({db}): C3 term=rho_mean={rho_mean:.4f} (OCC_lb={rho_lb:.4f}; "
            f"UniOP pairwise={p:.4f} -> C2 gate in final)"
        )
    return " ; ".join(parts)


def _parse_c3_counts(text: str) -> tuple[int | None, int | None, int | None, int | None]:
    if not text:
        return None, None, None, None
    m = re.search(
        r"operon context:\s*(\d+)\s*/\s*(\d+)\s+OCC-supported adjacencies,\s*(\d+)\s+neutral(?:,\s*(\d+)\s+unlinked)?",
        text,
    )
    if not m:
        return None, None, None, None
    supported = int(m.group(1))
    evaluated = int(m.group(2))
    neutral = int(m.group(3))
    unlinked = int(m.group(4)) if m.group(4) is not None else 0
    return supported, evaluated, neutral, unlinked


def _clean_c3_reasoning_prefix(text: str) -> str:
    if not text:
        return ""
    t = text.strip()
    if t.startswith("operon context:"):
        return ""
    return t


def _c3_reasoning_expanded(
    base_breakdown: str,
    operon_id: str,
    operon_member_count: str,
    pair_details: str,
    reference_organism_count: int,
) -> str:
    if operon_id == "NOT_APPLICABLE_NON_CODING":
        return "non-coding feature: C3 not applicable; neutral handling used by scorer."
    if operon_id == "NOT_IN_AN_OPERON":
        return "singleton gene (not in a multi-gene operon): C3 is neutral by design."

    supported, evaluated, neutral, unlinked = _parse_c3_counts(base_breakdown)
    prefix = _clean_c3_reasoning_prefix(base_breakdown)
    if supported is not None and evaluated is not None:
        explanation = (
            f"{operon_id} (total number of genes in this operon = {operon_member_count}): "
            f"evaluated adjacent gene-pairs with usable UniOP probability = {evaluated}; "
            f"pairs with OCC support (>0 reliability in cross-organism OCC reference) = {supported}; "
            f"pairs EXCLUDED from C3 (hypothetical member or no OCC evidence -> novel, not penalized) = {neutral}; "
            f"adjacent pairs missing usable UniOP probability = {unlinked}. "
            f"OCC reference population size = {reference_organism_count} organisms."
        )
    else:
        explanation = (
            f"{operon_id} (total number of genes in this operon = {operon_member_count}): "
            f"OCC reference population size = {reference_organism_count} organisms."
        )
    if prefix:
        explanation += f" Base classifier note: {prefix}."
    if pair_details:
        explanation += f" Pair-level calculations: {pair_details}"
    return explanation


def _parse_c4_formula(c4_formula: str) -> tuple[str, str]:
    if not c4_formula:
        return "", ""
    m = re.search(r"\(\s*(\d+)\s*/\s*5\s*\)\s*\*\s*([0-9]*\.?[0-9]+)", c4_formula)
    if not m:
        return c4_formula, ""
    n_tools = m.group(1)
    conflict_fraction = m.group(2)
    pretty = (
        "C4_EC_conflict_clearance = 1 - (EC_reporting_tools_count / 5) * EC_conflict_fraction "
        f"= 1 - ({n_tools}/5) * {conflict_fraction}"
    )
    reasoning = (
        f"EC_reporting_tools_count = {n_tools}. "
        f"EC_conflict_fraction = {conflict_fraction} = conflicting tool-pairs / all EC-reporting tool-pairs. "
        "C4 is the fraction of confidence retained after EC-conflict penalty."
    )
    return pretty, reasoning


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labeled-input", required=True)
    ap.add_argument("--operon-info-input", required=True)
    ap.add_argument("--operon-results-input", required=True)
    ap.add_argument("--merged-input", required=True)
    ap.add_argument("--hierarchy-tier-input", required=True)
    ap.add_argument("--confidence-final-input", required=True)
    ap.add_argument("--c4-input", required=True)
    ap.add_argument("--occ-reference", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    paths = [
        Path(args.labeled_input),
        Path(args.operon_info_input),
        Path(args.operon_results_input),
        Path(args.merged_input),
        Path(args.hierarchy_tier_input),
        Path(args.confidence_final_input),
        Path(args.c4_input),
        Path(args.occ_reference),
    ]
    for p in paths:
        if not p.is_file():
            print(f"[make-final-annotation-with-confidence] ERROR: missing input {p}", file=sys.stderr)
            raise SystemExit(1)

    labeled_by_fid = _load_tsv(Path(args.labeled_input))
    operon_by_fid = _load_tsv(Path(args.operon_info_input))
    merged_by_fid = _load_tsv(Path(args.merged_input))
    hierarchy_by_fid = _load_tsv(Path(args.hierarchy_tier_input))
    conf_by_fid = _load_tsv(Path(args.confidence_final_input))
    c4_by_fid = _load_tsv(Path(args.c4_input))

    ref = c3_occ.load_reference(Path(args.occ_reference))
    if not ref.get("finalized"):
        c3_occ.finalize_reference(ref)
    ref_organism_count = len(ref.get("organisms_added", []))

    pairwise = _pairwise_map(Path(args.operon_results_input))
    pair_prob_by_pair = {
        frozenset((fid, other)): p
        for fid, vals in pairwise.items()
        for other, p in vals
    }

    descriptor_by_fid = {
        fid: (row.get("best_consensus_product_descriptor") or "").strip()
        for fid, row in labeled_by_fid.items()
    }
    clean_by_fid = {fid: L.clean_descriptor(desc) for fid, desc in descriptor_by_fid.items()}
    uninf_by_fid = {fid: L.is_uninformative(clean) for fid, clean in clean_by_fid.items()}
    operon_members = _operon_members(labeled_by_fid, operon_by_fid)
    c3_pairs_by_operon = {
        oid: _c3_pair_details_for_operon(
            members, descriptor_by_fid, clean_by_fid, uninf_by_fid, pair_prob_by_pair, ref
        )
        for oid, members in operon_members.items()
    }

    prefixed_cols = _prefixed_columns(BASE_COLUMNS)
    col_map = dict(zip(BASE_COLUMNS, prefixed_cols))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=prefixed_cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for fid, cf in conf_by_fid.items():
            lab = labeled_by_fid.get(fid, {})
            op = operon_by_fid.get(fid, {})
            mg = merged_by_fid.get(fid, {})
            ht = hierarchy_by_fid.get(fid, {})

            operon_id = (op.get("operon_id") or "").strip()
            in_operon = "yes" if operon_id.startswith("operon_") else "no"
            best_desc = (lab.get("best_consensus_product_descriptor") or "").strip()

            preliminary = (cf.get("preliminary_confidence_c1_c4") or "").strip()
            final = (cf.get("final_confidence_operon_context") or "").strip()
            c1 = (cf.get("c1_score_database_coverage") or "").strip()
            c2 = (cf.get("c2_score_operon_probability") or "").strip()
            c3 = (cf.get("c3_score_operon_context") or "").strip()
            c4 = (cf.get("c4_score_EC_agreement") or cf.get("c4_score") or "").strip()

            try:
                prelim_f = float(preliminary) if preliminary else None
                final_f = float(final) if final else None
            except ValueError:
                prelim_f = None
                final_f = None
            if prelim_f is not None and final_f is not None:
                delta = final_f - prelim_f
                cctx_reason = (
                    f"Context contribution to confidence = adjusted_confidence - preliminary_confidence = {delta:+.4f}. "
                    f"Interpretation threshold = ±{_CONTEXT_THRESHOLD:.1f}: "
                    ">= +0.1 means context increases confidence; <= -0.1 means context decreases confidence; "
                    "between -0.1 and +0.1 means no material effect."
                )
            else:
                cctx_reason = "non-coding or not scored"

            c2_reason = _c2_reasoning(
                fid=fid,
                operon_id=operon_id,
                member_count=(op.get("operon_member_count") or ""),
                pairwise=pairwise,
                descriptor_by_fid=descriptor_by_fid,
                c2_score_text=c2,
            )
            c3_pairs = c3_pairs_by_operon.get(operon_id, "")
            c3_reason = _c3_reasoning_expanded(
                base_breakdown=(cf.get("c3_signal_breakdown") or "").strip(),
                operon_id=operon_id,
                operon_member_count=(op.get("operon_member_count") or ""),
                pair_details=c3_pairs,
                reference_organism_count=ref_organism_count,
            )

            c4_row = c4_by_fid.get(fid, {})
            c4_formula_src = (cf.get("c4_formula") or c4_row.get("c4_formula") or "").strip()
            c4_formula_pretty, c4_extra_reasoning = _parse_c4_formula(c4_formula_src)
            c4_reasoning_src = (cf.get("c4_reasoning") or c4_row.get("c4_reasoning") or "").strip()
            if c4_extra_reasoning:
                c4_reasoning = f"{c4_extra_reasoning} {c4_reasoning_src}".strip()
            else:
                c4_reasoning = c4_reasoning_src

            row = {
                col_map["organism_name"]: lab.get("organism_name", ""),
                col_map["domain"]: lab.get("domain", ""),
                col_map["envelope"]: cf.get("ENVELOPE_envelope_type", ""),
                col_map["gene_id"]: lab.get("gene_id", ""),
                col_map["RAST_feature_id"]: fid,
                col_map["FEATURE_TYPE"]: cf.get("feature_type", "") or lab.get("RAST_feature_type", ""),
                col_map["RAST_strand"]: lab.get("RAST_strand", ""),
                col_map["RAST_start"]: lab.get("gene_start", ""),
                col_map["RAST_end"]: lab.get("gene_end", ""),
                col_map["RAST_na_sequence"]: lab.get("na_seq", ""),
                col_map["RAST_aa_sequence"]: lab.get("aa_seq", ""),
                col_map["IS_IN_OPERON?"]: in_operon,
                col_map["best_consensus_product_descriptor"]: best_desc,
                col_map["best_consensus_product_descriptor_source"]: lab.get("product_descriptor_source", ""),
                col_map["best_consensus_product_descriptor_source_hierarchy_order"]: _hierarchy_order(ht),
                col_map["best_consensus_product_descriptor_source_audit_trail"]: lab.get("product_descriptor_audit_trail", ""),
                col_map["gene_description_fingerprint"]: _desc_fp(best_desc),
                col_map["specialized_database_hits"]: cf.get("specialized_db_hits", ""),
                col_map["localization_and_topology_hits"]: _localization_summary(cf),
                col_map["ENVELOPE"]: cf.get("ENVELOPE_envelope_type", ""),
                col_map["ENVELOPE_reason"]: cf.get("ENVELOPE_inference_basis", ""),
                col_map["C1_score_database_coverage"]: c1,
                col_map["C1_reasoning"]: cf.get("c1_score_reasoning", ""),
                col_map["UniOP_OPERON_id"]: operon_id,
                col_map["UniOP_operon_probability"]: cf.get("c2_uniop_probability_raw", ""),
                col_map["C2_score_pairwise_genes_UniOP_probability"]: c2,
                col_map["C2_score_formula"]: cf.get("c2_formula", ""),
                col_map["C2_score_reasoning"]: c2_reason,
                col_map["C3_score_operon_context"]: c3,
                col_map["C3_score_operon_context_hybrid"]: cf.get("c3_score_operon_context_hybrid", ""),
                col_map["C3_score_formula"]: cf.get("c3_formula", ""),
                col_map["C3_score_reasoning"]: c3_reason,
                col_map["C4_score_EC_conflict"]: c4,
                col_map["C4_score_formula"]: c4_formula_pretty or c4_formula_src,
                col_map["C4_score_reasoning"]: c4_reasoning,
                col_map["EC_EVIDENCE_STATUS"]: cf.get("c4_ec_agreement_status", ""),
                col_map["PRELIMINARY_confidence_C1_C4"]: preliminary,
                col_map["PRELIMINARY_confidence_C1_C4_formula"]: f"C1({c1}) * C4({c4})",
                col_map["PRELIMINARY_confidence_C1_C4_reasoning"]: "Preliminary confidence is C1 discounted by EC-conflict clearance (C4).",
                col_map["ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT"]: final,
                col_map["ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_hybrid"]: cf.get("final_confidence_operon_context_hybrid", ""),
                col_map["ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_formula"]: cf.get("confidence_score_formula", ""),
                col_map["ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_reasoning"]: cctx_reason,
                col_map["DOES_OPERON_CONTEXT_IMPROVE_CONFIDENCE?"]: cf.get("does_context_improve_confidence?", ""),
                col_map["CONFIDENCE_TIER"]: cf.get("confidence_tier", ""),
                col_map["CONFIDENCE_TIER_hybrid"]: cf.get("confidence_tier_hybrid", ""),
                col_map["NEEDS_REVIEW?"]: cf.get("needs_review", ""),
                col_map["NEEDS_REVIEW_REASON"]: cf.get("needs_review_reason", ""),
                col_map["BEST_PRODUCT_DESCRIPTOR(copied_here_for_convenience)"]: best_desc,
            }
            writer.writerow(row)
            n += 1

    print(f"[make-final-annotation-with-confidence] Wrote {n} rows -> {out_path}")


if __name__ == "__main__":
    main()
