#!/usr/bin/env python3
"""extract-ec-numbers.py — derives consolidated-only-ECs.tsv directly from
merge-all-columns.py's output (not any filtered/truncated view -- per-tool
EC evidence needs every column intact to extract from).

For each tool, EC numbers are pulled from whichever column actually
carries them: a few tools (RAST, dbCAN, eggNOG) have a dedicated EC
column already; the rest (KEGG, TCDB, COG, TIGRFAM, PFAM, PGAP, MEROPS,
UniProt, InterPro, GeneProp) embed EC numbers inline in their description
text (e.g. KEGG: "DNA polymerase III subunit beta [EC:2.7.7.7]").
Every listed tool gets its own {TOOL}_EC column regardless of whether
this particular genome happened to have any hits -- written empty rather
than omitted, so the schema is stable across genomes.

Normalises every match to the bare "X.X.X.X" / "X.X.X.-" form (drops the
"EC"/"EC:" prefix text and any surrounding parens/brackets), deduplicated,
semicolon-joined if a tool has more than one distinct EC number for a gene.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from _shared import IDENTITY_COLUMNS

csv.field_size_limit(10_000_000)

EC_PATTERN = re.compile(r"(\d+\.\d+\.\d+\.(?:\d+|-))")

# tool_name -> source column to scan for EC numbers in the merged table.
EC_SOURCE_COLUMN: dict[str, str] = {
    "RAST": "RAST_EC_numbers",
    "DBCAN": "DBCAN_ec_numbers",
    "EGGNOG": "EGGNOG_EC_numbers",
    "KEGG": "KEGG_description",
    "TCDB": "TCDB_description",
    "COG": "COG_description",
    "TIGRFAM": "TIGRFAM_description",
    "PFAM": "PFAM_description",
    "PGAP": "PGAP_description",
    "MEROPS": "MEROPS_description",
    "UNIPROT": "UNIPROT_description",
    "INTERPRO": "INTERPRO_description",
    "GENEPROP": "GENEPROP_description",
}


def extract_ec_numbers(text: str) -> str:
    if not text:
        return ""
    seen: set[str] = set()
    ordered: list[str] = []
    for match in EC_PATTERN.finditer(text):
        ec = match.group(1)
        if ec not in seen:
            seen.add(ec)
            ordered.append(ec)
    return ";".join(ordered)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="merge-all-columns.py's output TSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[extract-ec-numbers] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = list(IDENTITY_COLUMNS) + [f"{tool}_EC" for tool in EC_SOURCE_COLUMN]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = set(reader.fieldnames or [])
        missing_sources = [c for c in EC_SOURCE_COLUMN.values() if c not in header]
        if missing_sources:
            print(f"[extract-ec-numbers] WARNING: source columns not in merged input, "
                  f"those tools' EC columns will be empty: {missing_sources}", file=sys.stderr)

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        n = 0
        per_tool_hits = {tool: 0 for tool in EC_SOURCE_COLUMN}
        with open(output_path, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                out_row = {col: row.get(col, "") for col in IDENTITY_COLUMNS}
                for tool, src_col in EC_SOURCE_COLUMN.items():
                    ec_value = extract_ec_numbers(row.get(src_col, ""))
                    out_row[f"{tool}_EC"] = ec_value
                    if ec_value:
                        per_tool_hits[tool] += 1
                writer.writerow(out_row)
                n += 1

    print(f"[extract-ec-numbers] Wrote {n} rows × {len(out_columns)} cols → {output_path}")
    for tool, count in per_tool_hits.items():
        print(f"    {tool}_EC: {count} genes with at least one EC number")


if __name__ == "__main__":
    main()
