#!/usr/bin/env python3
"""verify_report.py -- independently re-derive each figure's numbers from the
finished scoring outputs and confirm the companion TSVs match.

Reviewers will not tolerate mistakes, so every figure is checked: (1) its PNG
and companion TSV(s) exist and are non-empty; (2) headline aggregates written
in the TSV are re-computed from the raw scoring tables and must agree. Writes a
plain-text PASS/FAIL report into the figures folder. Never raises on a data
mismatch by itself failing the pipeline -- it reports; the caller decides.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reportfig_lib as L  # noqa: E402

TOL = 1e-4


class Report:
    def __init__(self):
        self.checks: list[tuple[str, bool, str]] = []

    def ok(self, name, cond, detail=""):
        self.checks.append((name, bool(cond), detail))

    def approx(self, name, a, b, detail=""):
        try:
            cond = abs(float(a) - float(b)) <= TOL
        except (TypeError, ValueError):
            cond = False
        self.checks.append((name, cond, detail or f"{a} vs {b}"))

    @property
    def passed(self):
        return all(c for _, c, _ in self.checks)

    def write(self, path: Path, header: str):
        lines = [header, "=" * len(header), ""]
        n_pass = sum(c for _, c, _ in self.checks)
        lines.append(f"{n_pass}/{len(self.checks)} checks passed"
                     f"  ->  {'PASS' if self.passed else 'FAIL'}")
        lines.append("")
        for name, cond, detail in self.checks:
            mark = "PASS" if cond else "FAIL"
            lines.append(f"  [{mark}] {name}" + (f"  ({detail})" if detail else ""))
        path.write_text("\n".join(lines) + "\n")
        print("\n".join(lines))


def _files_ok(rep: Report, outdir: Path, pngs: list[str], tsvs: list[str]):
    for f in pngs:
        p = outdir / f
        rep.ok(f"png exists & non-empty: {f}", p.is_file() and p.stat().st_size > 0)
    for f in tsvs:
        p = outdir / f
        nonempty = p.is_file() and p.stat().st_size > 0
        rows = 0
        if nonempty:
            try:
                rows = len(pd.read_csv(p, sep="\t"))
            except Exception:
                nonempty = False
        rep.ok(f"tsv exists & has rows: {f}", nonempty and rows > 0, f"{rows} rows")


def verify_organism(run_root: Path, organism: str, outdir: Path) -> Report:
    rep = Report()
    genes = L.load_organism_genes(run_root, organism)

    _files_ok(rep, outdir,
              pngs=["fig01_confidence_tiers_and_scores.png", "fig02_confidence_stages.png",
                    "fig03_operon_context_by_size.png", "fig04_component_distributions.png",
                    "fig05_operon_neighbourhood.png", "fig06_top_reproduced_operons.png",
                    "fig07_most_unique_operons.png"],
              tsvs=["fig01_confidence_tier_counts.tsv", "fig02_confidence_stage_means.tsv",
                    "fig03_c3_by_operon_size.tsv", "fig04_component_summary.tsv",
                    "fig05_operon_neighbourhood.tsv"])

    # fig01: tier counts must match value_counts
    try:
        t = pd.read_csv(outdir / "fig01_confidence_tier_counts.tsv", sep="\t")
        src = genes[genes["confidence_tier"] != L.NONCODING_TIER]["confidence_tier"].value_counts()
        for _, r in t.iterrows():
            rep.approx(f"fig01 tier count [{r['confidence_tier']}]",
                       r["n_genes"], int(src.get(r["confidence_tier"], -1)))
    except Exception as e:
        rep.ok("fig01 tier counts re-derived", False, str(e))

    # fig02: stage means must match
    try:
        m = pd.read_csv(outdir / "fig02_confidence_stage_means.tsv", sep="\t").set_index("stage")
        mask = (pd.to_numeric(genes["preliminary_confidence_c1_c4"], errors="coerce").notna()
                & pd.to_numeric(genes["final_confidence_operon_context"], errors="coerce").notna())
        exp = {
            "preliminary_c1_c4": pd.to_numeric(genes["preliminary_confidence_c1_c4"], errors="coerce")[mask].mean(),
            "after_operon_context": pd.to_numeric(genes["final_confidence_operon_context"], errors="coerce")[mask].mean(),
            "final_score": pd.to_numeric(genes["confidence_score"], errors="coerce")[mask].mean(),
        }
        for stage, v in exp.items():
            rep.approx(f"fig02 mean [{stage}]", m.loc[stage, "mean_confidence"], v)
    except Exception as e:
        rep.ok("fig02 stage means re-derived", False, str(e))

    # fig04: component medians must match
    try:
        s = pd.read_csv(outdir / "fig04_component_summary.tsv", sep="\t").set_index("component")
        cmap = {"C1": "c1_score", "C2": "c2_score_from_operon_probability",
                "C3": "c3_score", "C4": "c4_score"}
        for comp, col in cmap.items():
            exp = float(pd.to_numeric(genes[col], errors="coerce").dropna().median())
            rep.approx(f"fig04 median [{comp}]", s.loc[comp, "median"], exp)
    except Exception as e:
        rep.ok("fig04 medians re-derived", False, str(e))

    return rep


def verify_global(run_root: Path, outdir: Path) -> Report:
    rep = Report()
    genes = L.load_all_genes(run_root)

    _files_ok(rep, outdir,
              pngs=["fig01_operon_context_by_size.png", "fig02_most_conserved_operons.png",
                    "fig03_operon_probability.png", "fig04_confidence_tiers_across_genomes.png",
                    "fig05_confidence_relationship.png", "fig06_component_pca.png"],
              tsvs=["fig01_operon_context_by_size.tsv", "fig02_most_conserved_operons.tsv",
                    "fig04_confidence_tiers_per_genome.tsv", "fig06_pca_variance.tsv"])

    # fig01: median C3 per size bin
    try:
        t = pd.read_csv(outdir / "fig01_operon_context_by_size.tsv", sep="\t")
        d = genes[genes["in_operon"]].copy()
        d["size_bin"] = d["operon_member_count"].map(L.size_bin)
        for _, r in t.iterrows():
            exp = float(pd.to_numeric(
                d.loc[d["size_bin"] == r["operon_size_bin"], "c3_score"],
                errors="coerce").dropna().median())
            rep.approx(f"global fig01 median C3 [{r['operon_size_bin']}]",
                       r["median_c3_score"], exp)
    except Exception as e:
        rep.ok("global fig01 medians re-derived", False, str(e))

    # fig06: variance ratios sum to ~100%
    try:
        v = pd.read_csv(outdir / "fig06_pca_variance.tsv", sep="\t")
        rep.approx("global fig06 variance sums to 100%",
                   v["variance_explained_pct"].sum(), 100.0)
    except Exception as e:
        rep.ok("global fig06 variance re-derived", False, str(e))

    return rep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--organism", default=None, help="omit for the global report")
    ap.add_argument("--figures-dir", default=None)
    args = ap.parse_args()
    run_root = Path(args.run_root)

    if args.organism:
        outdir = Path(args.figures_dir) if args.figures_dir else \
            run_root / args.organism / "scoring" / "figures"
        rep = verify_organism(run_root, args.organism, outdir)
        header = f"Report-figure verification — {args.organism}"
    else:
        outdir = Path(args.figures_dir) if args.figures_dir else \
            run_root / "scoring" / "figures" / "global"
        rep = verify_global(run_root, outdir)
        header = "Report-figure verification — global (pangenome)"

    rep.write(outdir / "_verification_report.txt", header)
    raise SystemExit(0 if rep.passed else 3)


if __name__ == "__main__":
    main()
