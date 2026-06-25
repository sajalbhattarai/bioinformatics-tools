#!/usr/bin/env python3
"""extract-all-tools-evalue.py — derives consolidated-all-tools-evalue.tsv
directly from merge-all-columns.py's output: id + description + every
available confidence metric (e-value, bit-score/identity, or whatever
that tool's own model actually produces), for every tool that has a
meaningful per-hit confidence signal.

Broader than "general tools" -- by request, also includes the specialist
databases (MEROPS, TCDB, dbCAN) and three InterPro member databases
chosen specifically for prokaryotic relevance:
  - HAMAP: curated specifically for bacterial/archaeal proteomes.
  - NCBIfam: NCBI's curated prokaryotic family HMMs (the modern
    superset of TIGRFAM used in NCBI's own genome annotation pipeline).
  - CDD: NCBI's Conserved Domain Database, broadly trusted across all
    organisms including prokaryotes.
  The other 15 InterPro member databases are deliberately excluded here:
  Coils/MobiDBLite are algorithm-based predictions with no e-value at
  all (not a "confidence" hit in this sense); AntiFam exists to flag
  spurious predictions, not annotate real ones; PANTHER/SMART/FunFam
  skew eukaryote-curated; Pfam/Gene3D/PIRSF/SUPERFAMILY are general-
  purpose, not prokaryote-specific, and Pfam is already covered by the
  standalone PFAM tool above. All 18 remain available in
  consolidated-merged-all-columns.tsv if more are needed later.
  The *shared* INTERPRO_id/INTERPRO_description (the cross-database
  entry an analysis-specific signature maps to) has no e-value of its
  own -- included anyway, conceptually treated as a high-confidence
  annotation in its own right (a curated InterPro entry, not a raw
  per-analysis score).

Confidence column(s) used per tool -- not every tool has the same shape:
  COG:      evalue, identity        (no bit-score field exists for COG)
  KEGG:     evalue, score, threshold (KofamScan's real significance test
                                      is score-vs-adaptive-threshold, not
                                      e-value alone)
  EGGNOG:   evalue, score
  PFAM:     full_seq_evalue, full_seq_score
  TIGRFAM:  full_seq_evalue, full_seq_score
  PGAP:     full_seq_evalue, full_seq_score
  UNIPROT:  evalue, bitscore
  MEROPS:   evalue, bitscore
  TCDB:     evalue, bitscore
  DBCAN:    no e-value/score at all -- its model is a 3-tool voting
            consensus; number_of_tools_hit is the closest confidence
            proxy it has (more sub-tools agreeing = more confidence)
  INTERPRO_HAMAP / INTERPRO_NCBIFAM / INTERPRO_CDD: evalue only

For multi-hit tools (PFAM/TIGRFAM/PGAP/InterPro sub-databases), every
confidence column is already "accession: value; accession2: value2" in
the merged table, same alignment as the id/description columns -- copied
through verbatim. "No hit" stays blank, never 0 -- a literal 0 would
misleadingly read as the best possible e-value rather than "no data".

Every column is always written, even if this genome had zero hits for
that tool -- empty rather than omitted.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from _shared import IDENTITY_COLUMNS

csv.field_size_limit(10_000_000)

# tool_label -> (id_col, description_col, [confidence columns in output order])
TOOL_CONFIDENCE_CONFIG: dict[str, tuple[str, str, list[str]]] = {
    "COG":     ("COG_id", "COG_description", ["COG_evalue", "COG_identity"]),
    "KEGG":    ("KEGG_id", "KEGG_description", ["KEGG_evalue", "KEGG_score", "KEGG_threshold"]),
    "EGGNOG":  ("EGGNOG_id", "EGGNOG_description", ["EGGNOG_evalue", "EGGNOG_score"]),
    "PFAM":    ("PFAM_id", "PFAM_description", ["PFAM_full_seq_evalue", "PFAM_full_seq_score"]),
    "TIGRFAM": ("TIGRFAM_id", "TIGRFAM_description", ["TIGRFAM_full_seq_evalue", "TIGRFAM_full_seq_score"]),
    "PGAP":    ("PGAP_id", "PGAP_description", ["PGAP_full_seq_evalue", "PGAP_full_seq_score"]),
    "UNIPROT": ("UNIPROT_id", "UNIPROT_description", ["UNIPROT_evalue", "UNIPROT_bitscore"]),
    "MEROPS":  ("MEROPS_id", "MEROPS_description", ["MEROPS_evalue", "MEROPS_bitscore"]),
    "TCDB":    ("TCDB_id", "TCDB_description", ["TCDB_evalue", "TCDB_bitscore"]),
    "DBCAN":   ("DBCAN_id", "DBCAN_description", ["DBCAN_number_of_tools_hit"]),
    "INTERPRO_HAMAP":   ("INTERPRO_HAMAP_id", "INTERPRO_HAMAP_description", ["INTERPRO_HAMAP_evalue"]),
    "INTERPRO_NCBIFAM": ("INTERPRO_NCBIFAM_id", "INTERPRO_NCBIFAM_description", ["INTERPRO_NCBIFAM_evalue"]),
    "INTERPRO_CDD":     ("INTERPRO_CDD_id", "INTERPRO_CDD_description", ["INTERPRO_CDD_evalue"]),
}

# Shared cross-database InterPro entry: no confidence column of its own.
INTERPRO_SHARED_COLUMNS: list[str] = ["INTERPRO_id", "INTERPRO_description"]


def build_column_plan() -> list[str]:
    plan: list[str] = list(IDENTITY_COLUMNS)
    for id_col, desc_col, confidence_cols in TOOL_CONFIDENCE_CONFIG.values():
        plan += [id_col, desc_col] + confidence_cols
    plan += INTERPRO_SHARED_COLUMNS
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[extract-all-tools-evalue] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    columns = build_column_plan()

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in columns if c not in (reader.fieldnames or [])]
        if missing:
            print(f"[extract-all-tools-evalue] WARNING: not in merged input, writing empty: {missing}",
                  file=sys.stderr)

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        with open(output_path, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                writer.writerow({col: row.get(col, "") for col in columns})
                n += 1

    print(f"[extract-all-tools-evalue] Wrote {n} rows × {len(columns)} cols → {output_path}")


if __name__ == "__main__":
    main()
