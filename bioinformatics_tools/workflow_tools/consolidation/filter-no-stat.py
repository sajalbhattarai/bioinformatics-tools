#!/usr/bin/env python3
"""filter-no-stat.py — stage 3 of consolidation.

Reduces merge-all-columns.py's table to identity, id/description columns and
core localisation/topology calls. Annotation id/description columns pass
through with their accession keys; topology columns (TMBED, Phobius, DeepSig)
lose their coordinate keys and keep distinct labels in first-seen order.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)


def dedup_labels_from_keyed_string(value: str) -> str:
    """Strips "key: " prefixes from a "; "-joined string and keeps distinct labels in order.

    '0-24: inside; 25-50: transmembrane_helix; 51-61: inside' -> 'inside; transmembrane_helix'
    """
    if not value:
        return ""
    seen: set[str] = set()
    ordered: list[str] = []
    for part in value.split("; "):
        label = part.split(": ", 1)[1] if ": " in part else part
        label = label.strip()
        if label and label not in seen:
            seen.add(label)
            ordered.append(label)
    return "; ".join(ordered)


def coalesce_phobius_topology(row: dict[str, str]) -> str:
    """Combines Phobius segment_type and segment_label into one deduplicated topology string.

    Each segment uses its label when present, else its type.
    """
    type_val = row.get("PHOBIUS_segment_type", "")
    label_val = row.get("PHOBIUS_segment_label", "")
    if not type_val and not label_val:
        return ""
    # Both columns share the same ordered segment keys, so entries pair by index.
    type_parts = type_val.split("; ") if type_val else []
    label_parts = label_val.split("; ") if label_val else []
    combined: list[str] = []
    for i, tpart in enumerate(type_parts):
        coord = tpart.split(": ", 1)[0] if ": " in tpart else ""
        tval = tpart.split(": ", 1)[1] if ": " in tpart else tpart
        lval = ""
        if i < len(label_parts):
            lpart = label_parts[i]
            lval = lpart.split(": ", 1)[1] if ": " in lpart else lpart
        chosen = lval.strip() if lval.strip() else tval.strip()
        combined.append(f"{coord}: {chosen}" if coord else chosen)
    return dedup_labels_from_keyed_string("; ".join(combined))


# ─── Column plan ────────────────────────────────────────────────────────────────
# Entries are (output column, source column, kind); kind is "copy",
# "dedup" (strip keys, dedup labels) or "phobius" (two-column coalesce).

IDENTITY_COLUMNS: list[tuple[str, str, str]] = [
    ("feature_id", "feature_id", "copy"),
    ("organism_name", "organism_name", "copy"),
    ("domain", "domain", "copy"),
    ("gene_id", "gene_id", "copy"),
    ("gene_start", "gene_start", "copy"),
    ("gene_end", "gene_end", "copy"),
    ("na_length", "na_length", "copy"),
    ("aa_length", "aa_length", "copy"),
    ("RAST_feature_type", "RAST_feature_type", "copy"),
    ("RAST_strand", "RAST_strand", "copy"),
    ("RAST_description", "RAST_description", "copy"),
    ("na_seq", "na_seq", "copy"),
    ("aa_seq", "aa_seq", "copy"),
]

ANNOTATION_ID_DESCRIPTION_TOOLS: list[str] = [
    "COG", "KEGG", "EGGNOG", "PFAM", "TIGRFAM", "PGAP",
    "MEROPS", "TCDB", "DBCAN", "UNIPROT",
    "INTERPRO_HAMAP", "INTERPRO_NCBIFAM", "INTERPRO_PFAM", "INTERPRO_PANTHER",
    "INTERPRO_PIRSF", "INTERPRO_PIRSR", "INTERPRO_GENE3D", "INTERPRO_CDD",
    "INTERPRO_COILS", "INTERPRO_PRINTS", "INTERPRO_SMART", "INTERPRO_SFLD",
    "INTERPRO_SUPERFAMILY", "INTERPRO_MOBIDB", "INTERPRO_PROSITE_PATTERNS",
    "INTERPRO_PROSITE_PROFILES", "INTERPRO_FUNFAM", "INTERPRO_ANTIFAM",
]

GENEPROP_COLUMNS: list[tuple[str, str, str]] = [
    ("GENEPROP_tigrfam_hit", "GENEPROP_tigrfam_hit", "copy"),
    ("GENEPROP_id", "GENEPROP_id", "copy"),
    ("GENEPROP_description", "GENEPROP_description", "copy"),
]

TOPOLOGY_COLUMNS: list[tuple[str, str, str]] = [
    ("TMBED_topology", "TMBED_topology", "dedup"),
    ("PHOBIUS_topology", "", "phobius"),
    ("DEEPSIG_feature_type", "DEEPSIG_feature_type", "dedup"),
]

SINGLE_ROW_PREDICTION_COLUMNS: list[tuple[str, str, str]] = [
    ("PSORTB_localization", "PSORTB_localization", "copy"),
    ("SIGNALP4_is_signal_peptide", "SIGNALP4_is_signal_peptide", "copy"),
    ("SIGNALP6_prediction", "SIGNALP6_prediction", "copy"),
    ("OPERON_id", "OPERON_id", "copy"),
]

ENVELOPE_COLUMNS: list[tuple[str, str, str]] = [
    ("ENVELOPE_envelope_type", "ENVELOPE_envelope_type", "copy"),
    ("ENVELOPE_inference_basis", "ENVELOPE_inference_basis", "copy"),
]


def build_column_plan(header: set[str]) -> list[tuple[str, str, str]]:
    """Returns the output column plan restricted to columns present in the merged header."""
    plan: list[tuple[str, str, str]] = list(IDENTITY_COLUMNS)
    for prefix in ANNOTATION_ID_DESCRIPTION_TOOLS:
        id_col, desc_col = f"{prefix}_id", f"{prefix}_description"
        if id_col in header:
            plan.append((id_col, id_col, "copy"))
        if desc_col in header:
            plan.append((desc_col, desc_col, "copy"))
    if "GENEPROP_id" in header:
        plan.extend(GENEPROP_COLUMNS)
    plan.extend(c for c in TOPOLOGY_COLUMNS if c[0] == "PHOBIUS_topology" or c[1] in header)
    plan.extend(c for c in SINGLE_ROW_PREDICTION_COLUMNS if c[1] in header)
    plan.extend(c for c in ENVELOPE_COLUMNS if c[1] in header)
    return plan


def transform_value(row: dict[str, str], source_col: str, kind: str) -> str:
    """Returns one output cell for a row according to the plan entry's kind."""
    if kind == "copy":
        return row.get(source_col, "")
    if kind == "dedup":
        return dedup_labels_from_keyed_string(row.get(source_col, ""))
    if kind == "phobius":
        return coalesce_phobius_topology(row)
    return ""


def main() -> None:
    """Streams the merged TSV through the column plan with csv and writes the filtered TSV."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[filter-no-stat] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = set(reader.fieldnames or [])
        plan = build_column_plan(header)
        out_columns = [out_col for out_col, _, _ in plan]

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        with open(output_path, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t")
            writer.writeheader()
            for row in reader:
                out_row = {
                    out_col: transform_value(row, src_col, kind)
                    for out_col, src_col, kind in plan
                }
                writer.writerow(out_row)
                n += 1

    print(f"[filter-no-stat] Wrote {n} rows × {len(out_columns)} cols → {output_path}")


if __name__ == "__main__":
    main()
