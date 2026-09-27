#!/usr/bin/env python3
"""add-ec-consensus.py — labeling stage: EC-number consensus.

Collects EC numbers from 13 tools in the merged table (dedicated EC columns or
regex on descriptions) and compares each tool's full EC set. "X.X.X.-" is
compatible with any fully resolved EC in the same subclass. Status is
full_consensus, majority_consensus, single_source, conflicting or no_evidence;
product_descriptor_ec_consistent says whether the winning label source's ECs
overlap the consensus (True/False/NA). Old 3.6.3.x vs new class 7 numbers
still count as conflicting.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

EC_PATTERN = re.compile(r"(\d+\.\d+\.\d+\.(?:\d+|-))")

# tool_name -> merged-table column scanned for EC numbers.
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

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]


def extract_ec_numbers(text: str) -> set[str]:
    """Returns the set of EC numbers found in text by regex."""
    if not text:
        return set()
    return set(EC_PATTERN.findall(text))


def collect_ec_evidence(merged_row: dict[str, str]) -> dict[str, set[str]]:
    """Returns {tool: EC set} for every tool with at least one EC in the merged row."""
    evidence: dict[str, set[str]] = {}
    for tool, src_col in EC_SOURCE_COLUMN.items():
        values = extract_ec_numbers(merged_row.get(src_col, ""))
        if values:
            evidence[tool] = values
    return evidence


def ec_compatible(a: str, b: str) -> bool:
    """Returns True when two ECs share the first three levels and the fourth matches or either is "-"."""
    a_parts, b_parts = a.split("."), b.split(".")
    if a_parts[:3] != b_parts[:3]:
        return False
    return a_parts[3] == "-" or b_parts[3] == "-" or a_parts[3] == b_parts[3]


def most_specific(ec_values: set[str]) -> str:
    """Returns a class representative, preferring a fully resolved EC over "X.X.X.-"."""
    resolved = sorted(v for v in ec_values if not v.endswith(".-"))
    return resolved[0] if resolved else sorted(ec_values)[0]


def cluster_ec_values(all_values: set[str]) -> list[set[str]]:
    """Groups EC strings into ec_compatible() classes by greedy pairwise assignment."""
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
    """Classifies per-tool EC class sets by set intersection and union.

    Returns (status, consensus_ecs, supporting_tool_count, total_distinct_classes,
    supporting_tools); consensus_ecs is ";"-joined when several ECs agree.
    """
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

    # No shared class: reports the best-supported class, status conflicting.
    class_support: dict[str, set[str]] = {}
    for tool, classes in tool_to_classes.items():
        for c in classes:
            class_support.setdefault(c, set()).add(tool)
    top_class, top_tools = max(class_support.items(), key=lambda kv: len(kv[1]))
    return "conflicting", top_class, len(top_tools), total_distinct, ";".join(sorted(top_tools))


def check_label_consistency(product_descriptor_source: str, evidence: dict[str, set[str]], consensus_ecs: str) -> str:
    """Returns "True"/"False"/"NA" for whether the winning tool's ECs overlap the consensus set."""
    if not consensus_ecs or product_descriptor_source not in EC_SOURCE_COLUMN:
        return "NA"
    winner_values = evidence.get(product_descriptor_source)
    if not winner_values:
        return "NA"
    consensus_set = set(consensus_ecs.split(";"))
    return "True" if winner_values & consensus_set else "False"


def build_evidence_string(evidence: dict[str, set[str]]) -> str:
    """Formats every tool's ECs as 'EGGNOG: 2.3.1.15; KEGG: 2.3.1.275'."""
    parts = [f"{tool}: {','.join(sorted(values))}" for tool, values in sorted(evidence.items())]
    return "; ".join(parts)


def main() -> None:
    """Joins labeled genes to merged-table EC evidence with csv and writes the EC consensus table."""
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
        "ec_agreement_status", "ec_supporting_tools", "ec_all_evidence", "product_descriptor_ec_consistent",
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
            product_descriptor_source = (
                row.get("product_descriptor_source", "") or row.get("label_source", "")
            )
            evidence = ec_evidence_by_gene.get(fid, {})
            status, value, count, distinct, supporting = classify_ec_agreement(evidence)
            consistent = check_label_consistency(product_descriptor_source, evidence, value)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["ec_consensus_number"] = value
            out_row["ec_consensus_count"] = str(count) if value else ""
            out_row["ec_total_distinct"] = str(distinct)
            out_row["ec_agreement_status"] = status
            out_row["ec_supporting_tools"] = supporting
            out_row["ec_all_evidence"] = build_evidence_string(evidence)
            out_row["product_descriptor_ec_consistent"] = consistent

            writer.writerow(out_row)
            status_counts[status] = status_counts.get(status, 0) + 1
            n += 1

    print(f"[add-ec-consensus] Wrote {n} genes → {output_path}")
    for status, count in sorted(status_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {status:20s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
