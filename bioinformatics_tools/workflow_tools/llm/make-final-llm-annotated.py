#!/usr/bin/env python3
"""make-final-llm-annotated.py — build FINAL_LLM_labeled-genes-annotated.tsv,
the publication-ready LLM-enriched annotation file.

Joins the full evidence file (labeled-genes-final-annotated.tsv, 41 cols) with
the LLM summary TSV on feature_id. Adds specialized-DB agreement flags, operon
coherence flags, and a composite flag_needs_review column.

Inputs:
  labeled-genes-final-annotated.tsv  (full evidence file, 41 cols)
  llm-summary.tsv                    (score-genes-llm.py output)

Output:
  FINAL_LLM_labeled-genes-annotated.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_OUT_COLS = [
    # ── Gene identity ──────────────────────────────────────────────────
    "feature_id",
    "organism_name",
    # ── Annotation result ─────────────────────────────────────────────
    "concordant_label",
    "label_derivation_logic",
    # ── Localization summary ──────────────────────────────────────────
    "SIGNALP6_prediction",
    "TMBED_topology",
    "PSORTB_localization",
    "ENVELOPE_envelope_type",
    "ENVELOPE_inference_basis",
    # ── Operon context ────────────────────────────────────────────────
    "operon_id",
    "operon_member_count",
    "operon_gene_position_in_operon",
    "operon_member_genes_with_labels_and_confidence_scores_in_position_order",
    # ── Mechanical confidence ─────────────────────────────────────────
    "confidence_score",
    "confidence_tier",
    "confidence_flag",
    "c1_score",
    "c2_score_from_operon_probability",
    "c3_score",
    "c4_score",
    # ── LLM confidence ────────────────────────────────────────────────
    "llm_assessment_score",
    "llm_c3_pathway_coherence_score",
    # ── Specialized DB agreement ──────────────────────────────────────
    "specialized_database_agreement_with_concordant_label",
    "flag_label_disagrees_with_specialized_databases",
    # ── Operon coherence ──────────────────────────────────────────────
    "operon_coherence_assessment",
    "flag_operon_possibly_incoherent",
    # ── Topology ──────────────────────────────────────────────────────
    "topology_consistency_with_localization_predictions",
    # ── Composite review flag ─────────────────────────────────────────
    "flag_needs_review",
    # ── Cross-genome fingerprint evidence ─────────────────────────────
    "fingerprint_hash",
    "fingerprint_consensus_label",
    "gene_fingerprint_exact_pattern_occurrence_count_in_database",
    "operon_label_ordered_pattern_occurrence_count_in_database",
    # ── LLM reasoning ─────────────────────────────────────────────────
    "llm_reasoning_text",
    # ── Provenance ────────────────────────────────────────────────────
    "full_evidence_source_file",
    "llm_evidence_source_file",
]


def _make_derivation_logic(row: dict) -> str:
    """Compact derivation: 'SOURCE [ID] (hierarchy_tier_name)'."""
    source    = (row.get("label_source") or "").strip()
    source_id = (row.get("label_source_id") or "").strip()
    tier      = (row.get("hierarchy_tier_name") or "").strip()
    if source and source_id:
        base = f"{source} [{source_id}]"
    elif source:
        base = source
    else:
        base = "unassigned"
    return f"{base} ({tier})" if tier else base


def _flag_needs_review(spec_agree: str, operon_coh: str,
                        topology: str, llm_verdict: str) -> str:
    if any([
        spec_agree == "disagrees",
        operon_coh == "no",
        topology == "inconsistent",
        llm_verdict == "DISAGREES",
    ]):
        reasons = []
        if spec_agree == "disagrees":
            reasons.append("specialized_db_conflict")
        if operon_coh == "no":
            reasons.append("operon_incoherent")
        if topology == "inconsistent":
            reasons.append("topology_mismatch")
        if llm_verdict == "DISAGREES":
            reasons.append("llm_score_disagrees")
        return "YES: " + "; ".join(reasons)
    return "no"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full-annotated-input", required=True,
                        help="labeled-genes-final-annotated.tsv (full 41-col evidence file)")
    parser.add_argument("--llm-summary-input", required=True,
                        help="llm-summary.tsv (score-genes-llm.py output)")
    parser.add_argument("--output", required=True,
                        help="FINAL_LLM_labeled-genes-annotated.tsv")
    args = parser.parse_args()

    final_path = Path(args.full_annotated_input)
    llm_path   = Path(args.llm_summary_input)
    out_path   = Path(args.output)

    for label, p in [("full-annotated", final_path), ("llm-summary", llm_path)]:
        if not p.is_file():
            print(f"[make-final-llm-annotated] ERROR: {label} not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    final_provenance = final_path.name
    llm_provenance   = llm_path.name

    # Load LLM summary keyed by feature_id
    llm_by_gene: dict[str, dict] = {}
    with open(llm_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            if fid and fid not in llm_by_gene:
                llm_by_gene[fid] = row

    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = n_llm_matched = 0
    with open(final_path, newline="") as in_fh, \
         open(out_path, "w", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter="\t")
        if reader.fieldnames is None:
            print("[make-final-llm-annotated] ERROR: full annotated input is empty", file=sys.stderr)
            raise SystemExit(1)

        writer = csv.DictWriter(out_fh, fieldnames=_OUT_COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            fid  = row.get("feature_id", "")
            llm  = llm_by_gene.get(fid, {})
            if llm:
                n_llm_matched += 1

            spec_agree  = llm.get("specialized_db_agreement", "")
            operon_coh  = llm.get("operon_coherent", "")
            topology    = llm.get("topology", "")
            llm_verdict = llm.get("llm_agreement_with_label", "")

            out_row = {
                # pass-through from FINAL_
                "feature_id":          fid,
                "organism_name":       row.get("organism_name", ""),
                "best_consensus_product_descriptor":    row.get("best_consensus_product_descriptor", ""),
                "label_derivation_logic": _make_derivation_logic(row),
                "SIGNALP6_prediction": row.get("SIGNALP6_prediction", ""),
                "TMBED_topology":      row.get("TMBED_topology", ""),
                "PSORTB_localization": row.get("PSORTB_localization", ""),
                "ENVELOPE_envelope_type":    row.get("ENVELOPE_envelope_type", ""),
                "ENVELOPE_inference_basis":  row.get("ENVELOPE_inference_basis", ""),
                "operon_id":                 row.get("operon_id", ""),
                "operon_member_count":       row.get("operon_member_count", ""),
                "operon_gene_position_in_operon": row.get("operon_gene_position_in_operon", ""),
                "operon_member_genes_with_labels_and_confidence_scores_in_position_order":
                    row.get("operon_member_genes_with_labels_and_confidence_scores_in_position_order", ""),
                "confidence_score":    row.get("confidence_score", ""),
                "confidence_tier":     row.get("confidence_score_tier", ""),
                "confidence_flag":     row.get("confidence_flag", ""),
                "c1_score":            row.get("c1_score", ""),
                "c2_score_from_operon_probability": row.get("c2_score_from_operon_probability", ""),
                "c3_score":            row.get("c3_score", ""),
                "c4_score":            row.get("c4_score", ""),
                "fingerprint_hash":    row.get("fingerprint_hash", ""),
                "fingerprint_consensus_label": row.get("fingerprint_consensus_label", ""),
                "gene_fingerprint_exact_pattern_occurrence_count_in_database":
                    row.get("gene_fingerprint_exact_pattern_occurrence_count_in_database", ""),
                "operon_label_ordered_pattern_occurrence_count_in_database":
                    row.get("operon_label_ordered_pattern_occurrence_count_in_database", ""),
                "full_evidence_source_file": final_provenance,
                # from LLM summary
                "llm_assessment_score":         llm.get("llm_assessment_score", ""),
                "llm_c3_pathway_coherence_score": llm.get("llm_c3", ""),
                "specialized_database_agreement_with_concordant_label": spec_agree,
                "flag_label_disagrees_with_specialized_databases":
                    "YES" if spec_agree == "disagrees" else ("no" if spec_agree else ""),
                "operon_coherence_assessment":  operon_coh,
                "flag_operon_possibly_incoherent":
                    "YES" if operon_coh == "no" else ("no" if operon_coh else ""),
                "topology_consistency_with_localization_predictions": topology,
                "flag_needs_review": _flag_needs_review(spec_agree, operon_coh,
                                                         topology, llm_verdict),
                "llm_reasoning_text":    llm.get("llm_reasoning_text", ""),
                "llm_evidence_source_file": llm_provenance,
            }
            writer.writerow(out_row)
            n_written += 1

    print(f"[make-final-llm-annotated] {n_written} genes written → {out_path}")
    print(f"  LLM summary matched: {n_llm_matched}/{n_written} genes")
    print(f"  {len(_OUT_COLS)} columns")


if __name__ == "__main__":
    main()
