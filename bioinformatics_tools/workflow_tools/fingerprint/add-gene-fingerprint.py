#!/usr/bin/env python3
"""add-gene-fingerprint.py — fingerprint stage, per-gene fingerprints.

Builds each gene's fingerprint from labeled-genes.tsv: a fixed-position,
" | "-joined list of RAST_description plus id and description slots for each
tool in _ALL_IDS_DESC_TOOLS (empty slots kept), its 16-hex SHA-256 hash, and
its consensus label. Writes five views combining hash, label and pattern;
the full-with-scores view adds C1-C4 and the confidence score from
labeled-genes-confidence-final.tsv when the gene has a scoring row.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_KEPT_COLUMNS = [
    "organism_name", "feature_id", "domain", "gene_id", "gene_start", "gene_end",
    "RAST_feature_type", "RAST_strand",
]

# Fixed slot order (same as assign-canonical-label.py's trust hierarchy); each
# tool gives an id and a description slot, kept even when empty.
_ALL_IDS_DESC_TOOLS = [
    "PGAP", "TIGRFAM", "HAMAP", "NCBIFAM", "PIRSF", "UNIPROT", "PFAM",
    "CDD", "KEGG", "EGGNOG", "COG", "MEROPS", "TCDB", "DBCAN",
]


def _hash16(s: str) -> str:
    """Returns the first 16 hex characters of the string's SHA-256 (hashlib)."""
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def _all_slots_empty(fingerprint_values: str) -> bool:
    """Returns True when every "{field}: {value}" slot has an empty value."""
    return all(slot.split(": ", 1)[-1] == "" for slot in fingerprint_values.split(" | "))


def build_fingerprint_values(row: dict[str, str]) -> str:
    """Builds the " | "-joined fingerprint string of "{field}: {value}" slots in fixed tool order."""
    slots = [f"RAST_description: {row.get('RAST_description', '').strip()}"]
    for tool in _ALL_IDS_DESC_TOOLS:
        slots.append(f"{tool}_id: {row.get(f'{tool}_all_ids', '').strip()}")
        slots.append(f"{tool}_description: {row.get(f'{tool}_all_descriptions', '').strip()}")
    return " | ".join(slots)


def build_scores_token(score_row: dict[str, str] | None) -> str:
    """Formats C1-C4 and the confidence score and tier as one "|"-joined token, empty without a scoring row."""
    if score_row is None:
        return ""
    return (
        f"C1:{score_row.get('c1_score', '')}"
        f"|C2:{score_row.get('c2_score_from_operon_probability', '')}"
        f"|C3:{score_row.get('c3_score', '')}"
        f"|C4:{score_row.get('c4_score', '')}"
        f"|CONFIDENCE:{score_row.get('confidence_score', '')}:{score_row.get('confidence_score_tier', '')}"
    )


def main() -> None:
    """Streams labeled-genes.tsv with csv and writes the five fingerprint views."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True, help="labeled-genes.tsv")
    parser.add_argument("--confidence-final-input", required=True,
                        help="labeled-genes-confidence-final.tsv (phase11/scoring) -- "
                             "only consumed by --output-full-with-scores")
    parser.add_argument("--output-hash-pattern", required=True)
    parser.add_argument("--output-hash-label", required=True)
    parser.add_argument("--output-label-pattern", required=True)
    parser.add_argument("--output-full", required=True)
    parser.add_argument("--output-full-with-scores", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    confidence_path = Path(args.confidence_final_input)
    for p in (labeled_path, confidence_path):
        if not p.is_file():
            print(f"[add-gene-fingerprint] ERROR: input not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    score_by_fid: dict[str, dict[str, str]] = {}
    with open(confidence_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            if fid:
                score_by_fid[fid] = row

    outputs = {
        "hash_pattern": Path(args.output_hash_pattern),
        "hash_label": Path(args.output_hash_label),
        "label_pattern": Path(args.output_label_pattern),
        "full": Path(args.output_full),
        "full_with_scores": Path(args.output_full_with_scores),
    }
    for p in outputs.values():
        p.parent.mkdir(parents=True, exist_ok=True)

    out_columns = _KEPT_COLUMNS + ["fingerprint"]

    n = 0
    no_hit_count = 0
    with open(labeled_path, newline="") as fh, \
         open(outputs["hash_pattern"], "w", newline="") as f_hp, \
         open(outputs["hash_label"], "w", newline="") as f_hl, \
         open(outputs["label_pattern"], "w", newline="") as f_lp, \
         open(outputs["full"], "w", newline="") as f_full, \
         open(outputs["full_with_scores"], "w", newline="") as f_fws:

        reader = csv.DictReader(fh, delimiter="\t")
        writers = {
            key: csv.DictWriter(f, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
            for key, f in (("hash_pattern", f_hp), ("hash_label", f_hl),
                           ("label_pattern", f_lp), ("full", f_full),
                           ("full_with_scores", f_fws))
        }
        for w in writers.values():
            w.writeheader()

        for row in reader:
            values = build_fingerprint_values(row)
            h = _hash16(values)
            label = row.get("best_consensus_product_descriptor", "")
            identity = {col: row.get(col, "") for col in _KEPT_COLUMNS}
            scores = build_scores_token(score_by_fid.get(row.get("feature_id", "")))

            writers["hash_pattern"].writerow({**identity, "fingerprint": f"pattern hash: {h} || fingerprint: {values}"})
            writers["hash_label"].writerow({**identity, "fingerprint": f"pattern hash: {h} || label: {label}"})
            writers["label_pattern"].writerow({**identity, "fingerprint": f"label: {label} || fingerprint: {values}"})
            writers["full"].writerow({**identity, "fingerprint": f"pattern hash: {h} || label: {label} || fingerprint: {values}"})
            full_with_scores = (f"pattern hash: {h} || label: {label} || scores: {scores} || fingerprint: {values}"
                                 if scores else f"pattern hash: {h} || label: {label} || fingerprint: {values}")
            writers["full_with_scores"].writerow({**identity, "fingerprint": full_with_scores})

            if _all_slots_empty(values):
                no_hit_count += 1
            n += 1

    print(f"[add-gene-fingerprint] Wrote {n} genes x 5 files:")
    for key, p in outputs.items():
        print(f"    {key}: {p}")
    print(f"    genes with zero tool hits (empty fingerprint values): {no_hit_count} ({100.0*no_hit_count/n:.1f}%)")


if __name__ == "__main__":
    main()
