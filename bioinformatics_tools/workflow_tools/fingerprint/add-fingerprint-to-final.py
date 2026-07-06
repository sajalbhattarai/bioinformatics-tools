#!/usr/bin/env python3
"""add-fingerprint-to-final.py — enrich labeled-genes-confidence-final.tsv
with per-gene fingerprint hit columns from labeled-genes-fingerprint-hash-label.tsv
and cross-genome frequency counts from the persistent fingerprint databases.

Reads (all required):
  labeled-genes-confidence-final.tsv     (phase11/scoring)
  labeled-genes-fingerprint-hash-label.tsv  (phase12/fingerprint)
  labeled-genes-operon-fingerprint.tsv   (phase12/run_operon_fingerprint)
  fingerprint-database.tsv               (cross-genome gene fingerprint pool)
  operon-fingerprint-database-label-ordered.tsv    (cross-genome operon pool)
  operon-fingerprint-database-label-composition.tsv

Writes labeled-genes-final-annotated.tsv — the single user-facing output
that combines identity, localization, operon, scoring, fingerprint evidence,
cross-genome frequency counts, and ordered operon member details in one
flat table.

New columns appended (all human-readable names):
  fingerprint_hash
  fingerprint_consensus_label
  gene_fingerprint_exact_pattern_occurrence_count_in_database
  gene_fingerprint_exact_pattern_organism_count_in_database
  gene_fingerprint_consensus_label_occurrence_count_in_database
  operon_label_ordered_pattern_occurrence_count_in_database
  operon_label_composition_pattern_occurrence_count_in_database
  operon_member_genes_with_labels_and_confidence_scores_in_position_order
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_FINGERPRINT_RE = re.compile(
    r"pattern hash:\s*([0-9a-f]+).*?\|\|\s*label:\s*(.+)", re.DOTALL
)
_OPERON_HASH_RE = {
    "label_ordered":     re.compile(r"operon hash by label \(ordered\):\s*(\S+)"),
    "label_composition": re.compile(r"operon hash by label \(composition\):\s*(\S+)"),
}
_NOT_IN_OPERON = "NOT_IN_AN_OPERON"

_NEW_COLUMNS = [
    "fingerprint_hash",
    "fingerprint_consensus_label",
    "gene_fingerprint_exact_pattern_occurrence_count_in_database",
    "gene_fingerprint_exact_pattern_organism_count_in_database",
    "gene_fingerprint_consensus_label_occurrence_count_in_database",
    "operon_label_ordered_pattern_occurrence_count_in_database",
    "operon_label_composition_pattern_occurrence_count_in_database",
    "operon_member_genes_with_labels_and_confidence_scores_in_position_order",
    "ec_all_evidence",
    "ec_supporting_tools",
]


def parse_fingerprint(raw: str) -> tuple[str, str]:
    """Extract (hash, consensus_label) from 'pattern hash: XXXX || label: YYYY'."""
    m = _FINGERPRINT_RE.search(raw)
    if m:
        return m.group(1).strip(), m.group(2).strip()
    return "", raw.strip()


def load_fingerprint_db(path: Path) -> dict[str, dict]:
    """Returns {hash: {frequency, label_frequency, organism_count}}."""
    result: dict[str, dict] = {}
    if not path.is_file():
        return result
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            h = row.get("fingerprint_hash", "")
            if not h:
                continue
            orgs = row.get("organisms", "")
            result[h] = {
                "frequency": int(row.get("fingerprint_frequency", 0) or 0),
                "label_frequency": int(row.get("fingerprint_label_frequency", 0) or 0),
                "organism_count": len([o for o in orgs.split("|") if o]) if orgs else 0,
            }
    return result


def load_operon_fp_db(path: Path) -> dict[str, int]:
    """Returns {operon_hash: operon_count}."""
    result: dict[str, int] = {}
    if not path.is_file():
        return result
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            h = row.get("fingerprint_hash", "")
            if h:
                result[h] = int(row.get("fingerprint_frequency", 0) or 0)
    return result


def parse_operon_hashes(operon_fp_str: str) -> tuple[str, str]:
    """Extract (label_ordered_hash, label_composition_hash) from the operon fingerprint string."""
    label_ordered = ""
    label_composition = ""
    m = _OPERON_HASH_RE["label_ordered"].search(operon_fp_str)
    if m:
        label_ordered = m.group(1).strip()
    m = _OPERON_HASH_RE["label_composition"].search(operon_fp_str)
    if m:
        label_composition = m.group(1).strip()
    return label_ordered, label_composition


def build_member_detail(
    members: list[tuple[int, str, str, str]]
) -> str:
    """
    members = [(position, feature_id, concordant_label, confidence_score), ...]
    sorted by position already.
    Returns a human-readable string like:
      Position 1: peg.1 [formate dehydrogenase alpha] confidence=0.87
        → Position 2: peg.2 [formate dehydrogenase beta] confidence=0.89
    """
    parts = []
    for pos, fid, label, score in members:
        score_str = f"{float(score):.4f}" if score and score not in ("", "N/A", "NOT_APPLICABLE_NON_CODING") else score
        parts.append(f"Position {pos}: {fid} [{label}] confidence={score_str}")
    return " → ".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--confidence-final-input", required=True,
                        help="labeled-genes-confidence-final.tsv (phase11)")
    parser.add_argument("--fingerprint-hash-label-input", required=True,
                        help="labeled-genes-fingerprint-hash-label.tsv (phase12)")
    parser.add_argument("--operon-fingerprint-input", required=True,
                        help="labeled-genes-operon-fingerprint.tsv (phase12 run_operon_fingerprint)")
    parser.add_argument("--fingerprint-database", required=True,
                        help="Cross-genome gene fingerprint pool (fingerprint-database.tsv)")
    parser.add_argument("--operon-fp-label-ordered-database", required=True,
                        help="operon-fingerprint-database-label-ordered.tsv")
    parser.add_argument("--operon-fp-label-composition-database", required=True,
                        help="operon-fingerprint-database-label-composition.tsv")
    parser.add_argument("--ec-consensus-input", required=True,
                        help="labeled-genes-ec-consensus.tsv (phase10/labeling) — "
                             "provides ec_all_evidence and ec_supporting_tools per gene")
    parser.add_argument("--output", required=True,
                        help="labeled-genes-final-annotated.tsv")
    args = parser.parse_args()

    for label, p in [
        ("confidence-final", Path(args.confidence_final_input)),
        ("fingerprint-hash-label", Path(args.fingerprint_hash_label_input)),
        ("operon-fingerprint", Path(args.operon_fingerprint_input)),
        ("ec-consensus", Path(args.ec_consensus_input)),
    ]:
        if not Path(p).is_file():
            print(f"[add-fingerprint-to-final] ERROR: {label} not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    # Load per-gene fingerprint (hash_label) → hash + consensus label
    fp_by_gene: dict[str, tuple[str, str]] = {}
    with open(args.fingerprint_hash_label_input, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            if fid and fid not in fp_by_gene:
                fp_by_gene[fid] = parse_fingerprint(row.get("fingerprint", ""))

    # Load operon fingerprint → per-gene operon hashes (label-ordered and label-composition)
    operon_hashes_by_gene: dict[str, tuple[str, str]] = {}
    with open(args.operon_fingerprint_input, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            fp_str = row.get("operon_fingerprint", "")
            if fid and fp_str and fid not in operon_hashes_by_gene:
                operon_hashes_by_gene[fid] = parse_operon_hashes(fp_str)

    # Load EC consensus evidence per gene (ec_all_evidence, ec_supporting_tools)
    ec_by_gene: dict[str, dict[str, str]] = {}
    with open(args.ec_consensus_input, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            if fid and fid not in ec_by_gene:
                ec_by_gene[fid] = {
                    "ec_all_evidence":    row.get("ec_all_evidence", ""),
                    "ec_supporting_tools": row.get("ec_supporting_tools", ""),
                }

    # Load cross-genome databases (may not exist yet on first run — tolerated)
    gene_fp_db = load_fingerprint_db(Path(args.fingerprint_database))
    operon_fp_label_ordered_db = load_operon_fp_db(Path(args.operon_fp_label_ordered_database))
    operon_fp_label_composition_db = load_operon_fp_db(Path(args.operon_fp_label_composition_database))

    # First pass over confidence_final: build operon-member index
    # operon_id → sorted list of (position, feature_id, concordant_label, confidence_score)
    operon_members: dict[str, list[tuple[int, str, str, str]]] = {}
    cf_rows: list[dict] = []

    with open(args.confidence_final_input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            print("[add-fingerprint-to-final] ERROR: confidence-final is empty", file=sys.stderr)
            raise SystemExit(1)
        out_cols = list(reader.fieldnames) + _NEW_COLUMNS
        for row in reader:
            cf_rows.append(row)
            oid = row.get("operon_id", _NOT_IN_OPERON) or _NOT_IN_OPERON
            if oid != _NOT_IN_OPERON:
                try:
                    pos = int(row.get("operon_gene_position_in_operon", "0") or "0")
                except ValueError:
                    pos = 0
                fid = row.get("feature_id", "")
                label = row.get("best_consensus_product_descriptor", "") or row.get("canonical_label", "")
                score = row.get("confidence_score", "")
                operon_members.setdefault(oid, []).append((pos, fid, label, score))

    # Sort each operon's members by position
    for oid in operon_members:
        operon_members[oid].sort(key=lambda t: t[0])

    # Second pass: write enriched output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with open(output_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in cf_rows:
            fid = row.get("feature_id", "")
            oid = row.get("operon_id", _NOT_IN_OPERON) or _NOT_IN_OPERON

            # Gene-level fingerprint
            fp_hash, fp_label = fp_by_gene.get(fid, ("", ""))
            gene_db = gene_fp_db.get(fp_hash, {}) if fp_hash else {}

            # Operon-level fingerprint hashes
            label_ordered_hash, label_composition_hash = operon_hashes_by_gene.get(fid, ("", ""))
            operon_ordered_count = operon_fp_label_ordered_db.get(label_ordered_hash, "")
            operon_composition_count = operon_fp_label_composition_db.get(label_composition_hash, "")

            # Operon member detail
            if oid != _NOT_IN_OPERON and oid in operon_members:
                member_detail = build_member_detail(operon_members[oid])
            else:
                member_detail = _NOT_IN_OPERON

            row["fingerprint_hash"] = fp_hash
            row["fingerprint_consensus_label"] = fp_label
            row["gene_fingerprint_exact_pattern_occurrence_count_in_database"] = (
                gene_db.get("frequency", "")
            )
            row["gene_fingerprint_exact_pattern_organism_count_in_database"] = (
                gene_db.get("organism_count", "")
            )
            row["gene_fingerprint_consensus_label_occurrence_count_in_database"] = (
                gene_db.get("label_frequency", "")
            )
            row["operon_label_ordered_pattern_occurrence_count_in_database"] = (
                operon_ordered_count
            )
            row["operon_label_composition_pattern_occurrence_count_in_database"] = (
                operon_composition_count
            )
            row["operon_member_genes_with_labels_and_confidence_scores_in_position_order"] = (
                member_detail
            )
            ec_row = ec_by_gene.get(fid, {})
            row["ec_all_evidence"]    = ec_row.get("ec_all_evidence", "")
            row["ec_supporting_tools"] = ec_row.get("ec_supporting_tools", "")
            writer.writerow(row)
            n_written += 1

    n_fp_matched = sum(1 for v in fp_by_gene.values() if v[0])
    n_operonic = sum(1 for row in cf_rows
                     if (row.get("operon_id") or _NOT_IN_OPERON) != _NOT_IN_OPERON)
    print(f"[add-fingerprint-to-final] {n_written} genes written → {output_path}")
    print(f"  fingerprint hashes matched: {n_fp_matched}/{len(fp_by_gene)}")
    print(f"  genes in operons (member detail generated): {n_operonic}/{n_written}")
    print(f"  gene fingerprint DB entries loaded: {len(gene_fp_db)}")
    print(f"  operon label-ordered DB entries loaded: {len(operon_fp_label_ordered_db)}")
    print(f"  operon label-composition DB entries loaded: {len(operon_fp_label_composition_db)}")


if __name__ == "__main__":
    main()
