#!/usr/bin/env python3
"""filter-no-stat.py — Stage 3 of margie_sb's consolidation pipeline.

Reads merge-all-columns.py's output (not the original per-tool files) and
produces a clean, stat-free, location-free, provenance-free view: just
identity, ids/descriptions, and core localization/topology calls.

Two different simplifications happen here, by tool shape:

  - Annotation-database tools (COG, KEGG, eggNOG, PFAM, TIGRFAM, PGAP,
    MEROPS, TCDB, dbCAN, UniProt, every interpro_<db>, GeneProp): their
    _id/_description columns are copied through as-is. merge-all-columns.py
    already keyed multi-hit values by accession ("PF00712.25: ...;
    PF02768.21: ..."), and that pairing stays -- you still need it to know
    which description belongs to which id when there's more than one hit.

  - Topology/localization tools (TMBED, Phobius, DeepSig): merge-all-
    columns.py keyed their multi-segment values by coordinate range
    ("0-24: inside; 25-50: transmembrane_helix; ..."). Here those
    coordinates get stripped entirely -- this view only answers "which
    kinds of region does this protein have," not "score" or "where."
    Repeated labels collapse to one entry, in first-occurrence order.

PSORTb/SignalP4/SignalP6/Operon are always exactly one row per gene in
merge-all-columns.py already (no coordinate-keying ever applied to them),
so their one core prediction column is copied through directly.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)


def dedup_labels_from_keyed_string(value: str) -> str:
    """'0-24: inside; 25-50: transmembrane_helix; 51-61: inside' -> 'inside; transmembrane_helix'

    Strips the "key: " prefix from each "; "-separated entry and keeps
    only the first occurrence of each distinct label.
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
    """Phobius splits topology across two columns (segment_type:
    DOMAIN/TRANSMEM, segment_label: CYTOPLASMIC/NON CYTOPLASMIC -- only
    populated for DOMAIN rows). Use the label when present, else the type,
    matching TMBED's single-column shape, then dedup like any other
    coordinate-keyed topology column."""
    type_val = row.get("PHOBIUS_segment_type", "")
    label_val = row.get("PHOBIUS_segment_label", "")
    if not type_val and not label_val:
        return ""
    # Both columns are coordinate-keyed the same way ("idx: value; ..."),
    # built from the same ordered segment list -- zip entry-by-entry.
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
#
# Each entry: (output_column_name, source_column_name_in_merged_tsv, kind)
#   kind == "copy"        -- passthrough, no transform
#   kind == "dedup"       -- strip coordinate/accession keys, dedup labels
#   kind == "phobius"     -- special two-column coalesce + dedup

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
    if kind == "copy":
        return row.get(source_col, "")
    if kind == "dedup":
        return dedup_labels_from_keyed_string(row.get(source_col, ""))
    if kind == "phobius":
        return coalesce_phobius_topology(row)
    return ""


def main() -> None:
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
