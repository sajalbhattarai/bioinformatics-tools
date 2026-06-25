#!/usr/bin/env python3
"""score-c1-tool-coverage.py — margie_sb phase11 (scoring), metric C1:
tool coverage.

Reads labeled-genes.tsv (phase10, READ-ONLY) and computes, per gene, how
many of the 12 decision tools (the ones in assign-canonical-label.py's
_EVALUATORS chain) gave an INFORMATIVE hit -- not just any hit. A bare
"Domain of unknown function" counts as a tool finding nothing useful,
same is_uninformative() gate assign-canonical-label.py already applies
when picking the best hit per tool.

C1 = informative_tool_count / 12

Deliberately excludes MEROPS/TCDB/DBCAN from the denominator -- those are
narrow specialist DBs that structurally cannot hit most genes (a DNA
polymerase will never match MEROPS), so including them would cap every
non-specialist gene below 1.0 regardless of how completely it's
annotated elsewhere. Cross-checking canonical_label against MEROPS/TCDB/
DBCAN confirmatory hits is deferred to a later LLM-based review pass, not
folded into this coverage metric.

Output (labeled-genes-c1-tool-coverage.tsv): identity columns, c1_score,
c1_informative_tool_count, c1_total_tools_considered, c1_informative_tools
(which specific tools counted, for traceability -- not just the bare
number), and c1_formula (the literal arithmetic as text, e.g.
"10/12 = 0.8333", so the score is auditable from this column alone).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

DECISION_TOOLS = ["PGAP", "TIGRFAM", "HAMAP", "NCBIFAM", "PIRSF", "UNIPROT",
                  "PFAM", "CDD", "KEGG", "EGGNOG", "COG", "RAST"]

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]

# Same uninformative-text gate as assign-canonical-label.py.
_UNINFORMATIVE = frozenset({
    "", "-", ".", "na", "n/a", "none", "null",
    "unknown", "uncharacterized", "uncharacterised", "putative", "predicted",
    "hypothetical protein", "conserved hypothetical protein", "conserved protein",
    "conserved domain protein", "predicted protein", "function unknown",
    "domain of unknown function", "general function prediction only",
    "poorly characterized", "open reading frame",
})
_UNINFORMATIVE_PREFIXES = (
    "domain of unknown function", "protein of unknown function",
    "family of unknown function", "region of unknown function",
    "repeat of unknown function", "module of unknown function",
    "duf", "upf", "uncharacteri", "putative uncharacteri",
    "conserved hypothetical", "hypothetical", "unknown protein",
    "unknown function", "orf", "pfam uncharacteri",
)
_UNKNOWN_FUNC_RE = re.compile(
    r'^\s*(?:(?:bacterial|viral|archaeal|eukaryotic|fungal|plant|marine|'
    r'transmembrane|integral membrane|membrane)\s+)?'
    r'(?:domain|protein|family|repeat|region|module)\s+of\s+unknown\s+function',
    re.IGNORECASE,
)
_DB_ID_HYPOTHETICAL_RE = re.compile(
    r'^(?:fig\d+|tigr\d+)[:\s].*(?:hypothetical|conserved hypothetical)', re.IGNORECASE,
)


def is_uninformative(val: str) -> bool:
    v = val.strip().lower()
    if not v or v in _UNINFORMATIVE:
        return True
    for prefix in _UNINFORMATIVE_PREFIXES:
        if v.startswith(prefix):
            return True
    if _UNKNOWN_FUNC_RE.match(v) or _DB_ID_HYPOTHETICAL_RE.match(v):
        return True
    return False


def compute_c1(labeled_row: dict[str, str]) -> tuple[float, list[str]]:
    informative_tools: list[str] = []
    for tool in DECISION_TOOLS:
        desc = labeled_row.get("RAST_description", "") if tool == "RAST" \
            else labeled_row.get(f"{tool}_best_hit_description", "")
        if desc and not is_uninformative(desc):
            informative_tools.append(tool)
    return len(informative_tools) / len(DECISION_TOOLS), informative_tools


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True, help="labeled-genes.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    if not labeled_path.is_file():
        print(f"[score-c1-tool-coverage] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)

    out_columns = _IDENTITY_COLUMNS + [
        "c1_score", "c1_informative_tool_count", "c1_total_tools_considered",
        "c1_informative_tools", "c1_formula",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    score_sum = 0.0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            c1, informative_tools = compute_c1(row)
            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c1_score"] = f"{c1:.4f}"
            out_row["c1_informative_tool_count"] = str(len(informative_tools))
            out_row["c1_total_tools_considered"] = str(len(DECISION_TOOLS))
            out_row["c1_informative_tools"] = ";".join(informative_tools)
            out_row["c1_formula"] = f"{len(informative_tools)}/{len(DECISION_TOOLS)} = {c1:.4f}"
            writer.writerow(out_row)
            score_sum += c1
            n += 1

    print(f"[score-c1-tool-coverage] Wrote {n} genes → {output_path}")
    print(f"    mean C1 = {score_sum / n:.4f}" if n else "    no genes scored")


if __name__ == "__main__":
    main()
