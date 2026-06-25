#!/usr/bin/env python3
"""extract-gene-symbols.py — derives consolidated-only-gene-symbols.tsv
directly from merge-all-columns.py's output (not any filtered/truncated
view -- this casts as wide a net as possible across every tool that
carries any kind of short symbolic/name evidence, not just strict formal
gene symbols).

Each tool's own native name-like column is copied through verbatim under
a uniform {TOOL}_gene_symbol name, e.g.:
  COG_gene_name          -> COG_gene_symbol      ("DnaN")
  EGGNOG_preferred_name  -> EGGNOG_gene_symbol    ("dnaN", or "-" when absent)
  UNIPROT_gene_name      -> UNIPROT_gene_symbol   ("dnaN")
  PGAP_name              -> PGAP_gene_symbol      ("dnan", sometimes a locus-like "MG010")
  PFAM_name              -> PFAM_gene_symbol      (Pfam family short name, e.g. "DNA_pol3_beta"
                                                    -- a domain name, not strictly a gene symbol,
                                                    but kept since it's still real symbolic evidence)
  TIGRFAM_name           -> TIGRFAM_gene_symbol
  RAST_BVBRC_Name        -> RAST_gene_symbol      (often empty on minimal genomes -- no other
                                                    RAST field is a true short symbol)
  UNIPROT_entry_name     -> UNIPROT_entry_name     (kept under its own name, not relabeled --
                                                    UniProt's structured ID code, e.g.
                                                    "DPO3B_MYCGE", a different kind of evidence
                                                    than UNIPROT_gene_name, both worth keeping)

Plus one derived column, RAST_description_gene_symbol: RAST/SEED-style
descriptions often carry a short symbol-like token in parentheses,
distinct from any "(EC ...)" parenthetical -- e.g. "Replicative DNA
helicase (DnaB) (EC 3.6.4.12)" -> "DnaB", or ribosomal proteins' universal
nomenclature alongside the bacterial name, "LSU ribosomal protein L13p
(L13Ae)" -> "L13Ae". The same parenthetical-extraction approach produces
mostly noise on every OTHER tool's description (PGAP/COG/PFAM/MEROPS
parentheses are mostly database IDs, organism names, or domain
abbreviations like "DUF240"/"PRK"/"HIT", not gene symbols) -- so this
extraction is RAST_description-only, not applied generically. Some false
positives still occur even here (e.g. "Asn" from "tRNA(Asn)", an amino
acid code, not a gene symbol) -- acceptable given the goal is broad
coverage, not perfect precision.

Every listed tool gets its column regardless of whether this particular
genome happened to have any hits -- written empty rather than omitted.
"-" (eggNOG's "no value" placeholder) is normalised to empty.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from _shared import IDENTITY_COLUMNS

csv.field_size_limit(10_000_000)

EC_PATTERN = re.compile(r"\d+\.\d+\.\d+\.(?:\d+|-)")
PAREN_PATTERN = re.compile(r"\(([^()]+)\)")
SYMBOL_LIKE_PATTERN = re.compile(r"^[A-Z][a-zA-Z0-9]{1,9}$")


def extract_rast_description_symbol(description: str) -> str:
    if not description:
        return ""
    for match in PAREN_PATTERN.finditer(description):
        inner = match.group(1).strip()
        if EC_PATTERN.search(inner):
            continue
        if SYMBOL_LIKE_PATTERN.match(inner):
            return inner
    return ""


# output_column_name -> source column in the merged table.
SYMBOL_SOURCE_COLUMN: dict[str, str] = {
    "RAST_gene_symbol": "RAST_BVBRC_Name",
    "COG_gene_symbol": "COG_gene_name",
    "EGGNOG_gene_symbol": "EGGNOG_preferred_name",
    "PFAM_gene_symbol": "PFAM_name",
    "TIGRFAM_gene_symbol": "TIGRFAM_name",
    "PGAP_gene_symbol": "PGAP_name",
    "UNIPROT_gene_symbol": "UNIPROT_gene_name",
    "UNIPROT_entry_name": "UNIPROT_entry_name",
}

NO_VALUE_PLACEHOLDERS: frozenset[str] = frozenset({"-", "."})


def clean_symbol(value: str) -> str:
    value = value.strip()
    return "" if value in NO_VALUE_PLACEHOLDERS else value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[extract-gene-symbols] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = list(IDENTITY_COLUMNS) + list(SYMBOL_SOURCE_COLUMN.keys()) + ["RAST_description_gene_symbol"]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = set(reader.fieldnames or [])
        missing_sources = [c for c in SYMBOL_SOURCE_COLUMN.values() if c not in header]
        if missing_sources:
            print(f"[extract-gene-symbols] WARNING: source columns not in merged input, "
                  f"those columns will be empty: {missing_sources}", file=sys.stderr)

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        per_col_hits = {col: 0 for col in list(SYMBOL_SOURCE_COLUMN) + ["RAST_description_gene_symbol"]}
        with open(output_path, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                out_row = {col: row.get(col, "") for col in IDENTITY_COLUMNS}
                for out_col, src_col in SYMBOL_SOURCE_COLUMN.items():
                    value = clean_symbol(row.get(src_col, ""))
                    out_row[out_col] = value
                    if value:
                        per_col_hits[out_col] += 1
                rast_symbol = extract_rast_description_symbol(row.get("RAST_description", ""))
                out_row["RAST_description_gene_symbol"] = rast_symbol
                if rast_symbol:
                    per_col_hits["RAST_description_gene_symbol"] += 1
                writer.writerow(out_row)
                n += 1

    print(f"[extract-gene-symbols] Wrote {n} rows × {len(out_columns)} cols → {output_path}")
    for col, count in per_col_hits.items():
        print(f"    {col}: {count} genes with a value")


if __name__ == "__main__":
    main()
