#!/usr/bin/env python3
"""score-c4-ec-agreement.py — margie_sb phase11 (scoring), metric C4: EC
conflict.

Reads labeled-genes-ec-consensus.tsv (phase10, add-ec-consensus.py's
output, READ-ONLY) and turns the per-tool EC evidence into a DIRECT,
GRADED MEASURE OF EC CONFLICT that discounts C1 (tool coverage). Also joins
in labeled-genes-confidence-tier.tsv (this same scoring phase,
score-confidence-tier.py's output, also READ-ONLY) purely for PROVENANCE,
so a reviewer looking at this one file can see the FULL step-by-step chain
that led here (hierarchy tier -> ec agreement -> combined score ->
confidence tier), not just C4's own number in isolation.

WHAT C4 MEASURES
  C4's only genuinely independent signal is CONFLICT — EC-capable tools that
  disagree on the EC number. It is a graded conflict penalty on C1, not a
  co-equal average term: silence and single-source agreement carry no
  independent information and never penalise C1.

DEFINITION  (conflict fraction R, then a participation-weighted clearance)
  Each database tool contributes ONE EC *set* (a tool reporting several EC
  numbers is still ONE node — its ECs are its set). Among the m tools that
  reported >=1 EC, form all C(m,2) tool-pairs:

      R = (# conflicting tool-pairs) / (m*(m-1)/2)          # 0 if m < 2

  A tool-pair (A,B) CONFLICTS iff, after matching the ECs they share, BOTH
  tools keep a private EC the other lacks — i.e. NEITHER EC set is a subset
  of the other under EC-hierarchy compatibility ('-' and missing trailing
  levels are wildcards, so 1.2.3.4 vs 1.2.3.- is compatible while 1.2.3.4
  vs 1.2.3.5 conflicts). m < 2 -> R = 0: silence and single-source are NOT
  conflict (nothing independent to disagree).

  The penalty applied to C1 is weighted by how much of the EC-capable panel
  actually reported an EC (participation), so there is NO tuned constant:

      c4_score = 1 - (m / 5) * R          # 5 = EC-capable databases
                                          #     (EGGNOG, RAST, KEGG,
                                          #      dbCAN, TCDB)

  and the final scorer (score-confidence-final.py) applies it
  multiplicatively:  base = C1 * c4_score. c4_score is thus the fraction of
  C1's coverage-confidence that SURVIVES the EC-conflict check: 1.0 = no
  conflict (base = C1), lower = more of the panel disagrees. The weight
  m/5 answers "why this penalty size?" with the fraction of EC-capable
  tools that weighed in — a conflict corroborated by more tools weighs
  more. A fully-contradicted gene seen by all 5 EC tools reaches
  c4_score = 0; in the current data m tops out at 4, so the largest weight
  applied is 4/5.

confidence_flag: ec_agreement_status == "conflicting" (the UPSTREAM
categorical call from add-ec-consensus.py) sets confidence_flag =
"needs_review" — the review TRIGGER, so it stays coherent with
score-confidence-tier.py and make-final-excel.py's row colouring. The
graded conflict fraction R drives the SCORE; the categorical status drives
the FLAG.

Output (labeled-genes-c4-ec-agreement.tsv): identity columns,
c4_ec_agreement_status, c4_ec_conflict_fraction (R), c4_n_conflict_pairs,
c4_n_total_pairs, c4_n_ec_tools (m), c4_score, c4_reasoning, c4_formula
(literal text, e.g. "c4_score = 1 - (m/5)*R = 1 - (3/5)*0.5000 = 0.7000"),
confidence_flag, then the joined provenance columns: hierarchy_tier_name,
hierarchy_tier_score, combined_score, confidence_tier, and
combined_score_formula.
"""
import argparse
import csv
import itertools
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# Number of databases capable of assigning an EC number (EGGNOG, RAST,
# KEGG, dbCAN, TCDB). The conflict penalty is weighted by how many of
# these actually reported an EC for a gene (m/5), so there is no tuned
# constant a reviewer could challenge.
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
    """'EGGNOG: 3.6.3.14; KEGG: 7.4.2.8,1.1.1.1' -> {tool: frozenset(ec)}.

    A tool with several ECs becomes ONE node whose value is the set of its ECs.
    """
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
    """'1.2.3.-' -> ('1','2','3',None); trailing '-'/blank levels become wildcards."""
    return tuple(None if p.strip().lower() in _EC_BLANK else p.strip()
                 for p in ec.split("."))


def ec_compatible(ec1: str, ec2: str) -> bool:
    """True iff ec1, ec2 do NOT contradict: equal on every mutually-specified
    level ('-' and missing trailing levels are wildcards -> EC hierarchy
    respected)."""
    for a, b in zip(ec_tuple(ec1), ec_tuple(ec2)):
        if a is None or b is None:
            continue
        if a != b:
            return False
    return True


def pair_conflicts(set_a, set_b) -> bool:
    """A tool-pair conflicts iff BOTH sides keep a private (incompatible) EC —
    i.e. neither EC set is a subset of the other under EC-compatibility."""
    a_left = [a for a in set_a if not any(ec_compatible(a, b) for b in set_b)]
    if not a_left:
        return False
    return any(not any(ec_compatible(a, b) for a in set_a) for b in set_b)


def c4_conflict_fraction(evmap: dict):
    """Return (R, n_conflict_pairs, n_total_pairs, n_ec_tools).

    R = conflicting tool-pairs / total EC-reporting tool-pairs; 0 when < 2
    tools reported an EC (no conflict is measurable).
    """
    tools = [t for t, s in evmap.items() if s]
    m = len(tools)
    if m < 2:
        return 0.0, 0, 0, m
    total = m * (m - 1) // 2
    conf = sum(1 for a, b in itertools.combinations(tools, 2)
               if pair_conflicts(evmap[a], evmap[b]))
    return conf / total, conf, total, m


def c4_clearance(conflict_fraction: float, n_ec_tools: int) -> float:
    """c4_score = 1 - (m/5)*R, clamped to [0,1].

    Weight m/5 = participation of the EC-capable panel; the final scorer
    applies base = C1 * c4_score."""
    weight = min(n_ec_tools, N_EC_CAPABLE_TOOLS) / float(N_EC_CAPABLE_TOOLS)
    return max(0.0, 1.0 - weight * conflict_fraction)


def _build_c4_reasoning(status: str, evidence: str, conflict_fraction: float,
                        n_conflict_pairs: int, n_total_pairs: int,
                        n_ec_tools: int, c4_score: float) -> str:
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

            # Review FLAG stays on the upstream categorical call (keeps
            # score-confidence-tier.py / make-final-excel.py coherent); the
            # graded conflict fraction drives the SCORE.
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
