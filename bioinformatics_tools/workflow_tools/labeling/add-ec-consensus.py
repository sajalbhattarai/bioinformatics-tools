#!/usr/bin/env python3
"""add-ec-consensus.py — margie_sb phase10 (labeling), Tier-1 consensus
layer: EC numbers.

Reads labeled-genes.tsv (assign-canonical-label.py's output, READ-ONLY --
this script never modifies or extends that file) and
consolidated-merged-all-columns.tsv (for EC evidence across multiple
independent tools), joins on feature_id, and writes its OWN standalone
table: a small set of identity/decision columns for cross-reference, plus
the EC consensus columns. Kept as a separate output file, not merged back
into labeled-genes.tsv -- the same "one core table, several focused
derived views" pattern as consolidation/extract-ec-numbers.py relative to
merge-all-columns.py. This keeps labeled-genes.tsv itself untouched and
auditable on its own, and keeps this consensus question (do independent
sources agree with each other?) clearly separate from the hierarchical
decision itself (which tool's text wins canonical_label).

EC EVIDENCE: same 13 sources as consolidation/extract-ec-numbers.py (a
few tools have a dedicated EC column, the rest have EC numbers embedded
inline in their description text).

AGREEMENT CLASSIFICATION compares each tool's FULL SET of reported ECs
against the others (not one EC at a time) -- needed for two real reasons:
  1. Bifunctional enzymes genuinely have >1 correct EC number (e.g. FolD,
     EC 1.5.1.5 AND EC 3.5.4.9, both real). When multiple tools each
     report the SAME pair, that's perfect agreement -- voting on
     individual EC numbers as if they were competing alternatives ties
     each EC at the same vote count and wrongly calls it "conflicting".
     Comparing full sets also fixes a single tool legitimately reporting
     2 activities with no second tool to compare against at all -- that's
     single_source, not a tie needing breaking.
  2. EC nomenclature itself has a wildcard: "2.3.4.-" means "subclass
     known, exact reaction not yet pinned down", and is the SAME enzyme
     as a more specific "2.3.4.1" from the same subclass, not a different
     one. ec_compatible() treats these as one compatibility class (the
     fully-resolved value is the reported representative) rather than
     two competing, unrelated ECs.
  Genuinely different EC numbers (e.g. EGGNOG: 2.3.1.15 vs KEGG: 2.3.1.275
  for the same gene) remain correctly flagged conflicting -- this only
  removes FALSE conflicts, not real ones. One known residual case this
  can't resolve: EC class 7 ("translocases") was only introduced in 2018
  to reclassify ABC-transporter ATPases that used to sit under 3.6.3.x --
  two tools of different vintage reporting e.g. "3.6.3.27" vs "7.3.2.1"
  for the same transporter may be the SAME call under old vs. new
  nomenclature, not a real disagreement, but detecting that needs an
  external EC-to-EC cross-reference table this pipeline doesn't have --
  left as a known limitation, still reported as "conflicting".

  full_consensus      -- every tool that reported anything agrees on
                         exactly the same set of ECs (zero disagreement).
  majority_consensus  -- tools share a non-empty common core, but at
                         least one tool reports an extra EC the others
                         don't confirm.
  single_source       -- exactly one tool reported anything at all -- real
                         evidence, but nothing independent to check it
                         against.
  conflicting         -- multiple tools reported genuinely different ECs
                         with no shared core at all.
  no_evidence         -- no tool reported any EC for this gene.

label_ec_consistent: does the WINNING tool's (label_source) own EC set
overlap at all with the consensus EC set? True/False/NA -- NA when
there's no EC evidence at all, or the winning tool has no EC of its own
to check (e.g. PGAP/TIGRFAM/HAMAP winning on a non-enzymatic structural
family).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

EC_PATTERN = re.compile(r"(\d+\.\d+\.\d+\.(?:\d+|-))")

# tool_name -> source column to scan for EC numbers in the merged table.
# Same mapping as consolidation/extract-ec-numbers.py.
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

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def extract_ec_numbers(text: str) -> set[str]:
    if not text:
        return set()
    return set(EC_PATTERN.findall(text))


def collect_ec_evidence(merged_row: dict[str, str]) -> dict[str, set[str]]:
    evidence: dict[str, set[str]] = {}
    for tool, src_col in EC_SOURCE_COLUMN.items():
        values = extract_ec_numbers(merged_row.get(src_col, ""))
        if values:
            evidence[tool] = values
    return evidence


def ec_compatible(a: str, b: str) -> bool:
    """True if two EC numbers are the same enzyme at a coarser/finer
    specificity -- "2.3.4.-" (subclass known, exact reaction not pinned
    down) is compatible with "2.3.4.1" (the same subclass, fully
    resolved). EC_PATTERN only ever places "-" in the 4th position (the
    first three levels are always digits), so compatibility only needs:
    same first three levels, and either side's 4th level is the wildcard
    "-" or they match exactly. Two fully-resolved-but-different 4th
    levels (e.g. "6.3.5.6" vs "6.3.5.7") are NOT compatible -- those are
    two distinct, specific reactions, not a coarse/fine pair."""
    a_parts, b_parts = a.split("."), b.split(".")
    if a_parts[:3] != b_parts[:3]:
        return False
    return a_parts[3] == "-" or b_parts[3] == "-" or a_parts[3] == b_parts[3]


def most_specific(ec_values: set[str]) -> str:
    """Within a compatibility class, prefer a fully-resolved member
    ("2.3.4.1") over a coarser "X.X.X.-" placeholder as the reported
    representative."""
    resolved = sorted(v for v in ec_values if not v.endswith(".-"))
    return resolved[0] if resolved else sorted(ec_values)[0]


def cluster_ec_values(all_values: set[str]) -> list[set[str]]:
    """Groups EC strings into ec_compatible() classes. Pairwise/greedy,
    not a full union-find -- fine here since a single gene has at most a
    handful of distinct EC values across all tools."""
    clusters: list[set[str]] = []
    for val in all_values:
        for cluster in clusters:
            if any(ec_compatible(val, member) for member in cluster):
                cluster.add(val)
                break
        else:
            clusters.append({val})
    return clusters


def classify_ec_agreement(evidence: dict[str, set[str]]) -> tuple[str, str, int, int, str]:
    """Returns (status, consensus_ecs, supporting_tool_count,
    total_distinct_classes, supporting_tools). consensus_ecs may be a
    ";"-joined list of more than one EC when the gene is genuinely
    multi-functional and every tool agrees on the same set."""
    if not evidence:
        return "no_evidence", "", 0, 0, ""

    all_values: set[str] = set()
    for values in evidence.values():
        all_values |= values
    clusters = cluster_ec_values(all_values)
    cluster_repr = {id(c): most_specific(c) for c in clusters}

    tool_to_classes: dict[str, set[str]] = {}
    for tool, values in evidence.items():
        touched = {cluster_repr[id(c)] for c in clusters if values & c}
        tool_to_classes[tool] = touched

    total_distinct = len(clusters)
    tools = list(tool_to_classes.keys())

    if len(tools) == 1:
        only_tool = tools[0]
        consensus = ";".join(sorted(tool_to_classes[only_tool]))
        return "single_source", consensus, 1, total_distinct, only_tool

    class_sets = list(tool_to_classes.values())
    intersection = set.intersection(*class_sets)
    union = set.union(*class_sets)

    if intersection and intersection == union:
        consensus = ";".join(sorted(intersection))
        return "full_consensus", consensus, len(tools), total_distinct, ";".join(sorted(tools))

    if intersection:
        consensus = ";".join(sorted(intersection))
        supporting = [t for t in tools if intersection <= tool_to_classes[t]]
        return "majority_consensus", consensus, len(supporting), total_distinct, ";".join(sorted(supporting))

    # No shared class at all -- report whichever class has the most tool
    # support as the representative "consensus" pick, but the status
    # itself stays conflicting since nothing is truly agreed.
    class_support: dict[str, set[str]] = {}
    for tool, classes in tool_to_classes.items():
        for c in classes:
            class_support.setdefault(c, set()).add(tool)
    top_class, top_tools = max(class_support.items(), key=lambda kv: len(kv[1]))
    return "conflicting", top_class, len(top_tools), total_distinct, ";".join(sorted(top_tools))


def check_label_consistency(label_source: str, evidence: dict[str, set[str]], consensus_ecs: str) -> str:
    """True/False/NA -- does the winning tool's own EC set overlap AT
    ALL with the consensus EC set? (Overlap, not exact match -- a
    winning tool reporting just one of two agreed bifunctional ECs is
    still consistent, not a mismatch.)"""
    if not consensus_ecs or label_source not in EC_SOURCE_COLUMN:
        return "NA"
    winner_values = evidence.get(label_source)
    if not winner_values:
        return "NA"
    consensus_set = set(consensus_ecs.split(";"))
    return "True" if winner_values & consensus_set else "False"


def build_evidence_string(evidence: dict[str, set[str]]) -> str:
    """Full per-tool breakdown of everything reported, e.g.
    'EGGNOG: 2.3.1.15; KEGG: 2.3.1.275' -- shows exactly what every tool
    said, including the ones that disagreed with the consensus, not just
    the winning count."""
    parts = [f"{tool}: {','.join(sorted(values))}" for tool, values in sorted(evidence.items())]
    return "; ".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True,
                        help="assign-canonical-label.py's output TSV (labeled-genes.tsv)")
    parser.add_argument("--merged-input", required=True,
                        help="merge-all-columns.py's output TSV (consolidated-merged-all-columns.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    merged_path = Path(args.merged_input)
    if not labeled_path.is_file():
        print(f"[add-ec-consensus] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)
    if not merged_path.is_file():
        print(f"[add-ec-consensus] ERROR: input not found: {merged_path}", file=sys.stderr)
        raise SystemExit(1)

    ec_evidence_by_gene: dict[str, dict[str, set[str]]] = {}
    with open(merged_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                ec_evidence_by_gene[fid] = collect_ec_evidence(row)

    out_columns = _IDENTITY_COLUMNS + [
        "ec_consensus_number", "ec_consensus_count", "ec_total_distinct",
        "ec_agreement_status", "ec_supporting_tools", "ec_all_evidence", "label_ec_consistent",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    status_counts: dict[str, int] = {}
    n = 0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            label_source = row.get("label_source", "")
            evidence = ec_evidence_by_gene.get(fid, {})
            status, value, count, distinct, supporting = classify_ec_agreement(evidence)
            consistent = check_label_consistency(label_source, evidence, value)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["ec_consensus_number"] = value
            out_row["ec_consensus_count"] = str(count) if value else ""
            out_row["ec_total_distinct"] = str(distinct)
            out_row["ec_agreement_status"] = status
            out_row["ec_supporting_tools"] = supporting
            out_row["ec_all_evidence"] = build_evidence_string(evidence)
            out_row["label_ec_consistent"] = consistent

            writer.writerow(out_row)
            status_counts[status] = status_counts.get(status, 0) + 1
            n += 1

    print(f"[add-ec-consensus] Wrote {n} genes → {output_path}")
    for status, count in sorted(status_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {status:20s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
