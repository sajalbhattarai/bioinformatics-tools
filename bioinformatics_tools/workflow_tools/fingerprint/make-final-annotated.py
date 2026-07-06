#!/usr/bin/env python3
"""make-final-annotated.py — build the publication-ready FINAL-scored-labeled-genes-annotated.tsv.

Joins the full evidence file (labeled-genes-final-annotated.tsv, from
add-fingerprint-to-final.py) with labeled-genes.tsv (gene details + product
descriptor hierarchy) and phobius_top1.tsv to produce a clean, self-contained
publication file.

Inputs (all required):
  labeled-genes-final-annotated.tsv      primary evidence (from add-fingerprint-to-final.py);
                                           includes C1-C4 scores + formulas, EC evidence, and
                                           fingerprint frequency columns
  labeling/labeled-genes.tsv             gene details: na/aa seq, positions, product_descriptor_hierarchy,
                                           product_descriptor_audit_trail, product_descriptor_confirmatory_summary
  phobius/phobius_top1.tsv               per-gene phobius topology summary
  fingerprint/labeled-genes-fingerprint-full.tsv
                                           full fingerprint: pattern hash || label ||
                                           annotation pattern (all tool slots)

Output:
  scoring/FINAL-scored-labeled-genes-annotated.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_OUT_COLS = [
    # ── Gene identity and details ──────────────────────────────────────
    "organism_name",
    "domain",
    "envelope",
    "feature_id",
    "gene_id",
    "na_seq",
    "aa_seq",
    "na_length",
    "aa_length",
    "gene_start",
    "gene_end",
    "RAST_feature_type",
    "RAST_strand",
    # ── Operon context ────────────────────────────────────────────────
    "operon_id",
    "operon_member_count",
    "operon_gene_position_in_operon",
    # ── Gene fingerprint ──────────────────────────────────────────────
    "gene_fingerprint_with_hash",
    # ── Annotation result ─────────────────────────────────────────────
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
    "product_descriptor_hierarchy",
    "product_descriptor_audit_trail",
    "product_descriptor_confirmatory_summary_specialized_tools",
    "product_descriptor_hierarchy_tier_name",
    # ── C1: database coverage ─────────────────────────────────────────
    "c1_score_database_coverage",
    "c1_score_reasoning",
    # ── C2: operon probability ────────────────────────────────────────
    "c2_score_operon_probability",
    "c2_score_reasoning",
    # ── C3: operon context ────────────────────────────────────────────
    "c3_score_operon_context",
    "c3_reasoning",
    # ── C4: EC number agreement ───────────────────────────────────────
    "c4_score_EC_agreement",
    "c4_reasoning",
    "c4_ec_agreement_status",
    # ── Two-stage confidence ──────────────────────────────────────────
    "preliminary_confidence_c1_c4",
    "final_confidence_operon_context",
    "does_context_improve_confidence?",
    "confidence_tier",
    "needs_review",
    "needs_review_reason",
    # ── Quick-view result copies (repeated at the very end so the label +
    #    coordinates sit next to the confidence columns, no scrolling back) ─
    "gene_id_copy",
    "gene_start_copy",
    "gene_end_copy",
    "best_consensus_product_descriptor_copy",
]


def _load_tsv(path: Path, key_col: str) -> dict[str, dict]:
    result: dict[str, dict] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            k = row.get(key_col, "")
            if k and k not in result:
                result[k] = row
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full-evidence-input", required=True,
                        help="labeled-genes-final-annotated.tsv (from add-fingerprint-to-final.py)")
    parser.add_argument("--labeling-genes-input", required=True,
                        help="labeling/labeled-genes.tsv")
    parser.add_argument("--phobius-top1-input", required=True,
                        help="phobius/phobius_top1.tsv")
    parser.add_argument("--fingerprint-full-input", required=True,
                        help="fingerprint/labeled-genes-fingerprint-full.tsv")
    parser.add_argument("--output", required=True,
                        help="scoring/FINAL-scored-labeled-genes-annotated.tsv")
    args = parser.parse_args()

    in_path      = Path(args.full_evidence_input)
    lab_path     = Path(args.labeling_genes_input)
    phob_path    = Path(args.phobius_top1_input)
    fp_full_path = Path(args.fingerprint_full_input)
    out_path     = Path(args.output)

    for label, p in [
        ("full-evidence",    in_path),
        ("labeling-genes",   lab_path),
        ("phobius-top1",     phob_path),
        ("fingerprint-full", fp_full_path),
    ]:
        if not p.is_file():
            print(f"[make-final-annotated] ERROR: {label} not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    lab_by_gene  = _load_tsv(lab_path,     "feature_id")
    fp_by_gene   = _load_tsv(fp_full_path, "feature_id")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(in_path, newline="") as in_fh, \
         open(out_path, "w", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter="\t")
        if reader.fieldnames is None:
            print("[make-final-annotated] ERROR: input is empty", file=sys.stderr)
            raise SystemExit(1)

        writer = csv.DictWriter(out_fh, fieldnames=_OUT_COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            fid = row.get("feature_id", "")
            lab = lab_by_gene.get(fid, {})
            fp  = fp_by_gene.get(fid, {})

            # gene_fingerprint_with_hash: use fp_pattern directly — it already contains
            # "pattern hash: XXXX || label: YYYY || fingerprint: ..." so prefixing again is redundant.
            fp_pattern = fp.get("fingerprint", "")
            if fp_pattern:
                gene_fp = fp_pattern
            else:
                fp_hash  = row.get("fingerprint_hash", "")
                fp_label = row.get("fingerprint_consensus_label", "")
                gene_fp  = f"{fp_hash} || {fp_label}" if fp_hash else ""

            out_row = {
                "organism_name":    row.get("organism_name", ""),
                "domain":           lab.get("domain", ""),
                "envelope":         row.get("ENVELOPE_envelope_type", ""),
                "feature_id":       fid,
                "gene_id":          lab.get("gene_id", ""),
                "na_seq":           lab.get("na_seq", ""),
                "aa_seq":           lab.get("aa_seq", ""),
                "na_length":        lab.get("na_length", ""),
                "aa_length":        lab.get("aa_length", ""),
                "gene_start":       lab.get("gene_start", ""),
                "gene_end":         lab.get("gene_end", ""),
                "RAST_feature_type": lab.get("RAST_feature_type", ""),
                "RAST_strand":      lab.get("RAST_strand", ""),
                "operon_id":                     row.get("operon_id", ""),
                "operon_member_count":            row.get("operon_member_count", ""),
                "operon_gene_position_in_operon": row.get("operon_gene_position_in_operon", ""),
                "gene_fingerprint_with_hash":     gene_fp,
                "best_consensus_product_descriptor":   row.get("best_consensus_product_descriptor", ""),
                "product_descriptor_source":           lab.get("product_descriptor_source", ""),
                "product_descriptor_source_id":        lab.get("product_descriptor_source_id", ""),
                "product_descriptor_hierarchy":        lab.get("product_descriptor_hierarchy", ""),
                "product_descriptor_audit_trail":      lab.get("product_descriptor_audit_trail", ""),
                "product_descriptor_confirmatory_summary_specialized_tools":
                    lab.get("product_descriptor_confirmatory_summary", ""),
                "product_descriptor_hierarchy_tier_name": row.get("hierarchy_tier_name", ""),
                "c1_score_database_coverage":    row.get("c1_score_database_coverage", ""),
                "c1_score_reasoning":            row.get("c1_score_reasoning", ""),
                "c2_score_operon_probability":   row.get("c2_score_operon_probability", ""),
                "c2_score_reasoning":            row.get("c2_score_reasoning", ""),
                "c3_score_operon_context":       row.get("c3_score_operon_context", ""),
                "c3_reasoning":                  row.get("c3_reasoning", ""),
                "c4_score_EC_agreement":         row.get("c4_score_EC_agreement", ""),
                "c4_reasoning":                  row.get("c4_reasoning", ""),
                "c4_ec_agreement_status":        row.get("c4_ec_agreement_status", ""),
                "preliminary_confidence_c1_c4":  row.get("preliminary_confidence_c1_c4", ""),
                "final_confidence_operon_context": row.get("final_confidence_operon_context", ""),
                "does_context_improve_confidence?": row.get("does_context_improve_confidence?", ""),
                "confidence_tier":               row.get("confidence_tier", ""),
                "needs_review":                  row.get("needs_review", ""),
                "needs_review_reason":           row.get("needs_review_reason", ""),
                # quick-view copies (same values as the columns near the top)
                "gene_id_copy":                  lab.get("gene_id", ""),
                "gene_start_copy":               lab.get("gene_start", ""),
                "gene_end_copy":                 lab.get("gene_end", ""),
                "best_consensus_product_descriptor_copy":     row.get("best_consensus_product_descriptor", ""),
            }
            writer.writerow(out_row)
            n += 1

    print(f"[make-final-annotated] {n} genes written → {out_path}")
    print(f"  {len(_OUT_COLS)} columns")


if __name__ == "__main__":
    main()
