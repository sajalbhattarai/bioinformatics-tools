#!/usr/bin/env python3
"""score-c4-ec-agreement.py — scoring stage, metric C4: EC conflict.

Each tool's EC numbers form one set; among the m tools reporting an EC, R is
the fraction of tool pairs where both sets keep an EC incompatible with the
other ("-" levels are wildcards). c4_score = 1 - (m/5)*R over the 5 EC-capable
databases, and the final scorer uses base = C1 * c4_score. The review flag comes
from add-ec-consensus.py's categorical "conflicting" status. Hierarchy-tier
columns are joined in for provenance.
"""
import argparse
import csv
import itertools
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# EC-capable databases (EGGNOG, RAST, KEGG, dbCAN, TCDB); the penalty weight is m/5.
N_EC_CAPABLE_TOOLS = 5

_EC_BLANK = {"", "-", "--", "n/a", "na", "none", "null", "*"}

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]


def parse_ec_evidence(ec_all_evidence: str) -> dict:
    """Parses 'EGGNOG: 3.6.3.14; KEGG: 7.4.2.8,1.1.1.1' into {tool: frozenset(ec)}."""
    out = {}
    if not ec_all_evidence:
        return out
    for part in ec_all_evidence.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        tool, ecs = part.split(":", 1)
        tool = tool.strip()
        s = {e.strip() for e in ecs.split(",") if e.strip().lower() not in _EC_BLANK}
        if s:
            out[tool] = frozenset(s)
    return out


def ec_tuple(ec: str) -> tuple:
    """Splits an EC into levels, '1.2.3.-' -> ('1','2','3',None), with blank levels as None."""
    return tuple(None if p.strip().lower() in _EC_BLANK else p.strip()
                 for p in ec.split("."))


def ec_compatible(ec1: str, ec2: str) -> bool:
    """Returns True when two ECs agree on every level both specify."""
    for a, b in zip(ec_tuple(ec1), ec_tuple(ec2)):
        if a is None or b is None:
            continue
        if a != b:
            return False
    return True


def pair_conflicts(set_a, set_b) -> bool:
    """Returns True when each EC set has an EC incompatible with every EC of the other."""
    a_left = [a for a in set_a if not any(ec_compatible(a, b) for b in set_b)]
    if not a_left:
        return False
    return any(not any(ec_compatible(a, b) for a in set_a) for b in set_b)


def c4_conflict_fraction(evmap: dict):
    """Returns (R, n_conflict_pairs, n_total_pairs, n_ec_tools) over itertools tool pairs; R is 0 below two tools."""
    tools = [t for t, s in evmap.items() if s]
    m = len(tools)
    if m < 2:
        return 0.0, 0, 0, m
    total = m * (m - 1) // 2
    conf = sum(1 for a, b in itertools.combinations(tools, 2)
               if pair_conflicts(evmap[a], evmap[b]))
    return conf / total, conf, total, m


def c4_clearance(conflict_fraction: float, n_ec_tools: int) -> float:
    """Returns c4_score = 1 - (m/5)*R, floored at 0."""
    weight = min(n_ec_tools, N_EC_CAPABLE_TOOLS) / float(N_EC_CAPABLE_TOOLS)
    return max(0.0, 1.0 - weight * conflict_fraction)


def _build_c4_reasoning(status: str, evidence: str, conflict_fraction: float,
                        n_conflict_pairs: int, n_total_pairs: int,
                        n_ec_tools: int, c4_score: float) -> str:
    """Builds the human-readable explanation of a gene's C4 score."""
    ev = evidence or "(no evidence)"
    if n_ec_tools == 0:
        return ("no EC number assigned by any tool — no conflict measurable — "
                "c4_score=1.0000 (no penalty on C1)")
    if n_ec_tools == 1:
        return (f"only one tool reported an EC ({ev}) — nothing independent to "
                f"conflict with — c4_score=1.0000 (no penalty on C1)")
    if n_conflict_pairs == 0:
        return (f"{n_ec_tools} tools reported ECs and none disagree "
                f"(0/{n_total_pairs} tool-pairs conflict): {ev} — "
                f"c4_score=1.0000 (no penalty on C1)")
    weight = min(n_ec_tools, N_EC_CAPABLE_TOOLS) / float(N_EC_CAPABLE_TOOLS)
    return (f"{n_conflict_pairs}/{n_total_pairs} EC tool-pairs conflict among "
            f"{n_ec_tools} tools (R={conflict_fraction:.4f}): {ev} — "
            f"penalty=(m/5)·R=({n_ec_tools}/5)·{conflict_fraction:.4f}"
            f"={weight * conflict_fraction:.4f} → c4_score={c4_score:.4f} "
            f"(discounts C1 by {100.0 * (1.0 - c4_score):.1f}%); "
            f"categorical status={status}")


def main() -> None:
    """Streams the EC consensus table with csv, scores C4 per gene and writes it with tier provenance."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ec-consensus-input", required=True, help="labeled-genes-ec-consensus.tsv")
    parser.add_argument("--confidence-tier-input", required=True,
                        help="score-confidence-tier.py's output TSV (labeled-genes-confidence-tier.tsv)"
                             " -- joined in for hierarchy/combined-score provenance only")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ec_path = Path(args.ec_consensus_input)
    tier_path = Path(args.confidence_tier_input)
    if not ec_path.is_file():
        print(f"[score-c4-ec-agreement] ERROR: input not found: {ec_path}", file=sys.stderr)
        raise SystemExit(1)
    if not tier_path.is_file():
        print(f"[score-c4-ec-agreement] ERROR: input not found: {tier_path}", file=sys.stderr)
        raise SystemExit(1)

    provenance_by_gene = {}
    with open(tier_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                provenance_by_gene[fid] = {
                    "hierarchy_tier_name": row.get("hierarchy_tier_name", ""),
                    "hierarchy_tier_score": row.get("hierarchy_tier_score", ""),
                    "ec_agreement_score": row.get("ec_agreement_score", ""),
                    "combined_score": row.get("combined_score", ""),
                    "confidence_tier": row.get("confidence_tier", ""),
                }

    out_columns = _IDENTITY_COLUMNS + [
        "c4_ec_agreement_status", "c4_ec_conflict_fraction",
        "c4_n_conflict_pairs", "c4_n_total_pairs", "c4_n_ec_tools",
        "c4_score", "c4_reasoning", "c4_formula", "confidence_flag",
        "hierarchy_tier_name", "hierarchy_tier_score", "combined_score",
        "confidence_tier", "combined_score_formula",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    flag_counts = {}
    n = 0
    n_conflicted = 0
    with open(ec_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            status = row.get("ec_agreement_status", "no_evidence")
            ec_evidence = row.get("ec_all_evidence", "")

            evmap = parse_ec_evidence(ec_evidence)
            conflict_fraction, n_conf, n_total, m = c4_conflict_fraction(evmap)
            c4 = c4_clearance(conflict_fraction, m)

            # The flag follows the categorical status; R drives the score.
            flag = "needs_review" if status == "conflicting" else "ok"

            prov = provenance_by_gene.get(fid, {
                "hierarchy_tier_name": "", "hierarchy_tier_score": "",
                "ec_agreement_score": "", "combined_score": "", "confidence_tier": "",
            })

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c4_ec_agreement_status"] = status
            out_row["c4_ec_conflict_fraction"] = f"{conflict_fraction:.4f}"
            out_row["c4_n_conflict_pairs"] = str(n_conf)
            out_row["c4_n_total_pairs"] = str(n_total)
            out_row["c4_n_ec_tools"] = str(m)
            out_row["c4_score"] = f"{c4:.4f}"
            out_row["c4_reasoning"] = _build_c4_reasoning(
                status, ec_evidence, conflict_fraction, n_conf, n_total, m, c4)
            out_row["c4_formula"] = (
                f"c4_score = 1 - (m/5)*R = 1 - ({m}/5)*{conflict_fraction:.4f} = {c4:.4f}"
            )
            out_row["confidence_flag"] = flag
            out_row["hierarchy_tier_name"] = prov["hierarchy_tier_name"]
            out_row["hierarchy_tier_score"] = prov["hierarchy_tier_score"]
            out_row["combined_score"] = prov["combined_score"]
            out_row["confidence_tier"] = prov["confidence_tier"]
            if prov["hierarchy_tier_score"] and prov["ec_agreement_score"]:
                out_row["combined_score_formula"] = (
                    f"hierarchy_tier_score({prov['hierarchy_tier_score']}) + "
                    f"ec_agreement_score({prov['ec_agreement_score']}) = "
                    f"combined_score({prov['combined_score']}) -> confidence_tier={prov['confidence_tier']}"
                )
            else:
                out_row["combined_score_formula"] = ""

            writer.writerow(out_row)
            flag_counts[flag] = flag_counts.get(flag, 0) + 1
            if conflict_fraction > 0.0:
                n_conflicted += 1
            n += 1

    print(f"[score-c4-ec-agreement] Wrote {n} genes → {output_path}")
    print(f"    EC-conflicted (R>0)  {n_conflicted:6d} ({100.0 * n_conflicted / n:.1f}%)"
          if n else "    (no genes)")
    for flag, count in sorted(flag_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {flag:15s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
