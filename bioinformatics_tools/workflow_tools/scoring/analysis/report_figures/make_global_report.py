#!/usr/bin/env python3
"""make_global_report.py -- pangenome (all-organism) report figures + TSVs.

Runs AFTER all genomes in a run are scored. Reads only the finished per-organism
scoring outputs across the run plus the depot pangenome operon reference
(read-only); writes only into  <run>/scoring/figures/global/ . Cannot affect
scoring. Presentation only -- no conclusions in any title/label. Each figure
emits a >=400 dpi PNG and a companion TSV with the exact plotted numbers.
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reportfig_lib as L  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_RUN_ROOT = None    # set in main(); used by the per-figure sources footer
_OPERON_DB = None


# ---------------------------------------------------------------------------
# fig01 -- operon-context score by operon size, across the pangenome
# ---------------------------------------------------------------------------
def fig01(genes: pd.DataFrame, outdir: Path, n_org: int) -> None:
    d = genes[genes["in_operon"]].copy()
    d["size_bin"] = d["operon_member_count"].map(L.size_bin)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    labels, ns, meds, means = L.box_by_bin(ax, d, "size_bin", "c3_score", L.SIZE_BIN_ORDER,
                                           L.GREEN, "operon-context score (C3)",
                                           "operon size (number of member genes)")
    ax.set_ylim(0, 1.02)
    L.finish(fig, f"Operon-context score by operon size — {n_org} genomes")
    L.draw_sources_footer(fig, _RUN_ROOT, L.global_source_lines(_RUN_ROOT, n_org))
    L.savefig(fig, outdir / "fig01_operon_context_by_size.png")
    L.write_tsv(pd.DataFrame({"operon_size_bin": labels, "n_genes": ns,
                              "median_c3_score": meds, "mean_c3_score": means}),
                outdir / "fig01_operon_context_by_size.tsv")


# ---------------------------------------------------------------------------
# fig02 -- most-conserved operons across the pangenome
# ---------------------------------------------------------------------------
def fig02(recurrence: dict, outdir: Path, n_org: int) -> None:
    rows = []
    for mio, rec in recurrence.items():
        members = [m.strip() for m in mio.split(" -> ") if m.strip()]
        if len(members) < 2 or not L.operon_members_informative(mio):
            continue  # show operons of named-function genes only
        rows.append({"members_in_order": mio, "n_members": len(members),
                     "pangenome_organisms": rec["organism_count"],
                     "pangenome_occurrences": rec["label_frequency"]})
    df = pd.DataFrame(rows)
    if df.empty:
        return
    top = df.sort_values(["pangenome_organisms", "pangenome_occurrences", "n_members"],
                         ascending=False).head(10).reset_index(drop=True)

    _CAP = 8
    n = len(top)
    _TRACK_IN, _LINE_IN = 0.72, 0.15
    members_per, units_per = [], []
    for _, r in top.iterrows():
        parts = [m.strip() for m in r["members_in_order"].split(" -> ") if m.strip()]
        members = [{"label": p, "operon_id": "operon_syn", "strand": "+"} for p in parts[:_CAP]]
        members_per.append(members)
        units_per.append(L.member_table_units(members))
    block_h = [_TRACK_IN + u * _LINE_IN + 0.3 for u in units_per]
    fig = plt.figure(figsize=(15.0, sum(block_h) + 1.2))
    outer = fig.add_gridspec(n, 1, height_ratios=block_h, hspace=0.4)
    for i, r in top.iterrows():
        members = members_per[i]
        inner = outer[i].subgridspec(2, 1, height_ratios=[_TRACK_IN, units_per[i] * _LINE_IN],
                                     hspace=0.1)
        axT = fig.add_subplot(inner[0])
        axTab = fig.add_subplot(inner[1])
        entries = L.draw_gene_track(axT, members, show_gaps=False, min_span=_CAP + 1)
        shown = f"  (first {_CAP} of {int(r['n_members'])} shown)" if r["n_members"] > _CAP else ""
        axT.set_title(f"{int(r['n_members'])}-gene operon  ·  in "
                      f"{int(r['pangenome_organisms'])} pangenome genomes{shown}",
                      fontsize=11.5, pad=6, loc="left", fontweight="bold")
        L.render_gene_table(axTab, entries, fontsize=7.8, show_location=False,
                            show_scores=False)
    fh = max(fig.get_figheight(), 3.0)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.subplots_adjust(top=1 - (1.15 + (0.3 if L._PROVENANCE else 0)) / fh)
    st = fig.suptitle("Operons of named-function genes shared across the most genomes — "
                      f"{n_org} genomes", y=1 - 0.45 / fh, fontweight="bold")
    pv = L.draw_provenance_line(fig, 1 - 0.85 / fh)
    if pv is not None:
        fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [st, pv]
    L.draw_sources_footer(fig, _RUN_ROOT,
                          L.global_source_lines(_RUN_ROOT, n_org, scored=False,
                                                operon_db=_OPERON_DB))
    L.savefig(fig, outdir / "fig02_most_conserved_operons.png")
    L.write_tsv(top, outdir / "fig02_most_conserved_operons.tsv")


# ---------------------------------------------------------------------------
# fig03 -- operon probability: distribution + relation to operon-context score
# ---------------------------------------------------------------------------
def fig03(genes: pd.DataFrame, outdir: Path, n_org: int) -> None:
    d = genes.copy()
    c2 = pd.to_numeric(d["c2_score_from_operon_probability"], errors="coerce")
    c3 = pd.to_numeric(d["c3_score"], errors="coerce")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.0, 5.8))
    axA.hist(c2.dropna(), bins=np.linspace(0, 1, 26), color=L.ORANGE,
             edgecolor="white", linewidth=0.5)
    axA.set_xlabel("operon probability (C2)")
    axA.set_ylabel("number of genes")
    axA.set_xlim(0, 1)
    L.set_title(axA, "Distribution of operon probability")
    L.panel_letter(axA, "a")

    # mean operon-context score across operon-probability bins
    bins = np.linspace(0, 1, 11)
    d = d.assign(_c2=c2, _c3=c3).dropna(subset=["_c2", "_c3"])
    d["c2_bin"] = pd.cut(d["_c2"], bins, include_lowest=True)
    grp = d.groupby("c2_bin", observed=True)["_c3"]
    centers = [iv.mid for iv in grp.mean().index]
    axB.plot(centers, grp.mean().values, "-o", color=L.GREEN, lw=2.2, markersize=7)
    axB.fill_between(centers, grp.quantile(0.25).values, grp.quantile(0.75).values,
                     color=L.GREEN, alpha=0.18)
    axB.set_xlabel("operon probability (C2)")
    axB.set_ylabel("operon-context score (C3)")
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)
    L.set_title(axB, "Operon-context score across operon-probability bins")
    L.panel_letter(axB, "b")

    L.finish(fig, f"Operon probability — {n_org} genomes")
    L.draw_sources_footer(fig, _RUN_ROOT, L.global_source_lines(_RUN_ROOT, n_org))
    L.savefig(fig, outdir / "fig03_operon_probability.png")
    L.write_tsv(pd.DataFrame({"c2_bin_mid": centers,
                              "mean_c3": grp.mean().values,
                              "q25_c3": grp.quantile(0.25).values,
                              "q75_c3": grp.quantile(0.75).values,
                              "n_genes": grp.size().values}),
                outdir / "fig03_c3_by_operon_probability.tsv")


# ---------------------------------------------------------------------------
# fig04 -- confidence tiers per organism + pooled
# ---------------------------------------------------------------------------
def fig04(genes: pd.DataFrame, outdir: Path, n_org: int) -> None:
    d = genes[genes["confidence_tier"] != L.NONCODING_TIER].copy()
    tiers = [t for t in L.CONF_TIER_ORDER if t in set(d["confidence_tier"])]
    orgs = sorted(d["organism"].unique(), key=lambda o: -(d["organism"] == o).sum())
    frac = []
    for org in orgs:
        sub = d[d["organism"] == org]
        frac.append([100.0 * (sub["confidence_tier"] == t).mean() for t in tiers])
    frac = np.array(frac)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16.0, 7.2),
                                   gridspec_kw={"width_ratios": [2.4, 1]})
    left = np.zeros(len(orgs))
    for j, t in enumerate(tiers):
        axA.barh(range(len(orgs)), frac[:, j], left=left, height=0.78,
                 color=L.CONF_TIER_COLOR[t], edgecolor="white", label=t)
        left += frac[:, j]
    axA.set_yticks(range(len(orgs)))
    axA.set_yticklabels([L.short_organism(o) for o in orgs], fontsize=8.6)
    axA.set_xlabel("share of coding genes (%)")
    axA.set_xlim(0, 100)
    axA.invert_yaxis()
    axA.grid(False)
    L.set_title(axA, "Confidence tiers per genome")
    L.panel_letter(axA, "a")
    handles, labels = axA.get_legend_handles_labels()
    _fh = max(fig.get_figheight(), 3.0)
    fig.legend(handles, labels, ncol=len(tiers), loc="upper center",
               bbox_to_anchor=(0.5, 1 - 0.92 / _fh), fontsize=10, frameon=False)

    pooled = [int((d["confidence_tier"] == t).sum()) for t in tiers]
    axB.bar(tiers, pooled, color=[L.CONF_TIER_COLOR[t] for t in tiers],
            edgecolor="white", width=0.72)
    for i, v in enumerate(pooled):
        axB.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9.5)
    axB.set_ylabel("number of genes")
    axB.set_ylim(0, max(pooled) * 1.15)
    axB.grid(False)
    L.set_title(axB, "All genomes pooled")
    L.panel_letter(axB, "b")

    # extra band: suptitle + a shared tier legend both sit above the panels
    L.finish(fig, f"Annotation confidence across genomes — {n_org} genomes", band=1.95)
    L.draw_sources_footer(fig, _RUN_ROOT, L.global_source_lines(_RUN_ROOT, n_org))
    L.savefig(fig, outdir / "fig04_confidence_tiers_across_genomes.png")
    tsv = pd.DataFrame(frac, columns=[f"pct_{t}" for t in tiers])
    tsv.insert(0, "organism", orgs)
    L.write_tsv(tsv, outdir / "fig04_confidence_tiers_per_genome.tsv")


# ---------------------------------------------------------------------------
# fig05 -- preliminary -> operon-adjusted -> final (relationship surface)
# ---------------------------------------------------------------------------
def fig05(genes: pd.DataFrame, outdir: Path, n_org: int) -> None:
    prelim = pd.to_numeric(genes["preliminary_confidence_c1_c4"], errors="coerce")
    c3 = pd.to_numeric(genes["c3_score"], errors="coerce")
    final = pd.to_numeric(genes["confidence_score"], errors="coerce")
    m = prelim.notna() & final.notna()

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.4))

    hb = axA.hexbin(prelim[m], final[m], gridsize=34, cmap="Blues", mincnt=1)
    axA.plot([0, 1], [0, 1], color=L.VERMILLION, lw=1.4, ls="--")
    axA.set_xlabel("preliminary confidence (C1–C4)")
    axA.set_ylabel("final confidence score")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    fig.colorbar(hb, ax=axA, label="number of genes", shrink=0.85)
    axA.grid(False)
    L.set_title(axA, "Preliminary vs final")
    L.panel_letter(axA, "a")

    # (b) mean final by preliminary bin
    bins = np.linspace(0, 1, 11)
    dd = pd.DataFrame({"p": prelim[m], "f": final[m]})
    dd["pb"] = pd.cut(dd["p"], bins, include_lowest=True)
    g = dd.groupby("pb", observed=True)["f"]
    centers = [iv.mid for iv in g.mean().index]
    axB.plot(centers, g.mean().values, "-o", color=L.BLUE, lw=2.3, markersize=7)
    axB.fill_between(centers, g.quantile(0.25).values, g.quantile(0.75).values,
                     color=L.BLUE, alpha=0.18)
    axB.plot([0, 1], [0, 1], color="#999999", lw=1.2, ls="--")
    axB.set_xlabel("preliminary confidence (C1–C4)")
    axB.set_ylabel("final confidence score")
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)
    L.set_title(axB, "Mean final score by preliminary bin")
    L.panel_letter(axB, "b")

    # (c) response surface: mean final over (preliminary, operon-context) grid
    dd2 = pd.DataFrame({"p": prelim[m], "c3": c3[m], "f": final[m]}).dropna()
    edges = np.linspace(0, 1, 11)
    dd2["pb"] = pd.cut(dd2["p"], edges, labels=False, include_lowest=True)
    dd2["cb"] = pd.cut(dd2["c3"], edges, labels=False, include_lowest=True)
    grid = dd2.groupby(["cb", "pb"], observed=True)["f"].mean().unstack("pb")
    grid = grid.reindex(index=range(10), columns=range(10))
    im = axC.imshow(grid.values, origin="lower", aspect="auto", cmap="viridis",
                    extent=(0, 1, 0, 1), vmin=0, vmax=1)
    fig.colorbar(im, ax=axC, label="mean final score", shrink=0.85)
    axC.set_xlabel("preliminary confidence (C1–C4)")
    axC.set_ylabel("operon-context score (C3)")
    axC.grid(False)
    L.set_title(axC, "Final score: C1–C4 × operon context")
    L.panel_letter(axC, "c")

    L.finish(fig, f"How the final confidence relates to its inputs — {n_org} genomes", top=0.88)
    L.draw_sources_footer(fig, _RUN_ROOT, L.global_source_lines(_RUN_ROOT, n_org))
    L.savefig(fig, outdir / "fig05_confidence_relationship.png")
    L.write_tsv(pd.DataFrame({"preliminary_bin_mid": centers,
                              "mean_final": g.mean().values,
                              "q25_final": g.quantile(0.25).values,
                              "q75_final": g.quantile(0.75).values,
                              "n_genes": g.size().values}),
                outdir / "fig05_final_by_preliminary.tsv")


# ---------------------------------------------------------------------------
# fig06 -- PCA of the confidence components (what varies with confidence)
# ---------------------------------------------------------------------------
def fig06(genes: pd.DataFrame, outdir: Path, n_org: int) -> None:
    cols = [("C1", "c1_score"), ("C2", "c2_score_from_operon_probability"),
            ("C3", "c3_score"), ("C4", "c4_score")]
    X = np.column_stack([pd.to_numeric(genes[c], errors="coerce").values for _, c in cols])
    finalv = pd.to_numeric(genes["confidence_score"], errors="coerce").values
    ok = ~np.isnan(X).any(axis=1) & ~np.isnan(finalv)
    X, finalv = X[ok], finalv[ok]
    if len(X) < 10:
        return
    scores, loadings, evr = L.pca_svd(X)
    names = [n for n, _ in cols]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.5, 6.4))

    axA.bar(range(1, len(evr) + 1), evr * 100,
            color=[L.CATEGORICAL[i % len(L.CATEGORICAL)] for i in range(len(evr))],
            edgecolor="white", linewidth=1.2, width=0.6)
    for i, e in enumerate(evr):
        axA.text(i + 1, e * 100, f"{e*100:.0f}%", ha="center", va="bottom", fontsize=9.5)
    axA.set_xticks(range(1, len(evr) + 1))
    axA.set_xlabel("principal component")
    axA.set_ylabel("share of variance explained (%)")
    axA.set_ylim(0, max(evr) * 100 * 1.18)
    axA.grid(False)
    L.set_title(axA, "How much each component axis explains")
    L.panel_letter(axA, "a")

    idx = np.argsort(finalv)
    sc = axB.scatter(scores[idx, 0], scores[idx, 1], c=finalv[idx], s=7,
                     cmap="viridis", alpha=0.6, vmin=0, vmax=1)
    scale = np.abs(scores[:, :2]).max() * 0.9
    for j, nm in enumerate(names):
        axB.arrow(0, 0, loadings[j, 0] * scale, loadings[j, 1] * scale,
                  color=L.VERMILLION, width=0.003, head_width=scale * 0.05, zorder=5)
        axB.text(loadings[j, 0] * scale * 1.12, loadings[j, 1] * scale * 1.12, nm,
                 color=L.VERMILLION, fontsize=12, fontweight="bold", ha="center", zorder=6)
    fig.colorbar(sc, ax=axB, label="final confidence score", shrink=0.85)
    axB.set_xlabel(f"component 1 ({evr[0]*100:.0f}% of variance)")
    axB.set_ylabel(f"component 2 ({evr[1]*100:.0f}% of variance)")
    axB.grid(False)
    L.set_title(axB, "Genes in component space (arrows = C1–C4)")
    L.panel_letter(axB, "b")

    L.finish(fig, f"What varies together in the confidence components — {n_org} genomes", top=0.9)
    L.draw_sources_footer(fig, _RUN_ROOT, L.global_source_lines(_RUN_ROOT, n_org))
    L.savefig(fig, outdir / "fig06_component_pca.png")
    ld = pd.DataFrame(loadings, columns=[f"PC{i+1}" for i in range(loadings.shape[1])])
    ld.insert(0, "component", names)
    L.write_tsv(ld, outdir / "fig06_pca_loadings.tsv")
    L.write_tsv(pd.DataFrame({"principal_component": range(1, len(evr) + 1),
                              "variance_explained_pct": evr * 100}),
                outdir / "fig06_pca_variance.tsv")


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--operon-db", default=str(L.DEFAULT_OPERON_DB))
    ap.add_argument("--output-dir", default=None,
                    help="default: <run>/scoring/figures/global")
    args = ap.parse_args()

    run_root = Path(args.run_root)
    outdir = Path(args.output_dir) if args.output_dir else \
        run_root / "scoring" / "figures" / "global"
    outdir.mkdir(parents=True, exist_ok=True)

    global _RUN_ROOT, _OPERON_DB
    _RUN_ROOT, _OPERON_DB = run_root, Path(args.operon_db)
    L.apply_style()
    organisms = L.discover_organisms(run_root)
    genes = L.load_all_genes(run_root, organisms)
    # Scope pangenome recurrence + the "N genomes" caption to the ACTUAL OCC pool
    # (occ_reference.pkl's organisms) -- the persistent baseline the C3/OCC scores
    # were computed against -- NOT just this run's scored-so-far organisms (which
    # under-reports mid-run). This keeps the 21-organism baseline visible and grows
    # as new organisms enter the OCC. Falls back to the run if the OCC is unreadable.
    occ_pool = L.load_occ_organisms() or set(organisms)
    recurrence = L.load_operon_recurrence(Path(args.operon_db), restrict_to=occ_pool)
    # Global report spans the WHOLE pool (no leave-one-out -- it is not about any
    # single organism), so the caption is the full OCC: 21 genomes. Count + gene/
    # operon tallies come from the OCC genome-stats sidecar (complete even mid-run).
    pool_list = sorted(occ_pool)
    n_org = len(pool_list)
    pstats = L.aggregate_pool_stats(L.load_pool_stats(), pool_list)
    L.set_provenance(L.provenance_text(n_org, pstats))

    jobs = [
        ("fig01", lambda: fig01(genes, outdir, n_org)),
        ("fig02", lambda: fig02(recurrence, outdir, n_org)),
        ("fig03", lambda: fig03(genes, outdir, n_org)),
        ("fig04", lambda: fig04(genes, outdir, n_org)),
        ("fig05", lambda: fig05(genes, outdir, n_org)),
        ("fig06", lambda: fig06(genes, outdir, n_org)),
    ]
    failures = []
    for name, fn in jobs:
        try:
            fn()
        except Exception:
            failures.append(name)
            print(f"[make_global_report] {name} FAILED:\n{traceback.format_exc()}",
                  file=sys.stderr)
    L.write_sources_manifest(outdir, run_root, organisms, Path(args.operon_db))
    print(f"[make_global_report] {n_org} genomes: {len(jobs) - len(failures)}/{len(jobs)} "
          f"figures ok → {outdir}")
    if failures:
        print(f"[make_global_report] FAILURES: {failures}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
