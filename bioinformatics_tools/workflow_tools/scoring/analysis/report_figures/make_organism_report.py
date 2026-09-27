#!/usr/bin/env python3
"""make_organism_report.py -- per-organism report figures and TSVs.

Runs after a genome is scored; reads its scoring outputs and the depot operon
reference and writes <run>/<organism>/scoring/figures/. Each figure shows the
organism against the pangenome as a matplotlib PNG (>= 400 dpi) with a TSV of
the plotted numbers.
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


_RUN_ROOT = None    # set in main(); used by the per-figure sources footer
_OPERON_DB = None


def _table_units(members) -> float:
    """Returns the number of text lines render_gene_table needs for these members."""
    return L.member_table_units(members)


# ---------------------------------------------------------------------------
# fig01 -- confidence tiers and the confidence-score distribution
# ---------------------------------------------------------------------------
def fig01(genes: pd.DataFrame, outdir: Path, org_label: str) -> None:
    """Draws genes per confidence tier and the final-score histogram (numpy)."""
    coding = genes[genes["confidence_tier"] != L.NONCODING_TIER]
    counts = coding["confidence_tier"].value_counts()
    tiers = [t for t in L.CONF_TIER_ORDER if t in counts.index]
    vals = [int(counts.get(t, 0)) for t in tiers]
    cols = [L.CONF_TIER_COLOR[t] for t in tiers]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.6, 5.6))

    bars = axA.bar(tiers, vals, color=cols, edgecolor="white", linewidth=1.2, width=0.72)
    for b, v in zip(bars, vals):
        axA.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                 va="bottom", fontsize=10, color=L.INK)
    axA.set_ylabel("number of genes")
    axA.set_xlabel("confidence tier")
    axA.set_ylim(0, max(vals) * 1.15 if vals else 1)
    axA.grid(False, axis="x")
    L.set_title(axA, "Genes per confidence tier")
    L.panel_letter(axA, "a")

    score = pd.to_numeric(coding["confidence_score"], errors="coerce").dropna()
    axB.hist(score, bins=np.linspace(0, 1, 26), color=L.BLUE, edgecolor="white", linewidth=0.5)
    med = float(score.median())
    axB.axvline(med, color=L.VERMILLION, lw=2, label=f"median = {med:.2f}")
    axB.set_xlabel("final confidence score")
    axB.set_ylabel("number of genes")
    axB.set_xlim(0, 1)
    axB.legend(loc="best", fontsize=10, frameon=True, facecolor="white",
               edgecolor="#bbbbbb", framealpha=0.95)
    L.set_title(axB, "Distribution of the final confidence score")
    L.panel_letter(axB, "b")

    L.finish(fig, "Confidence of gene annotations", organism=org_label)
    L.draw_sources_footer(fig, _RUN_ROOT, L.organism_source_lines(_RUN_ROOT, org_label))
    L.savefig(fig, outdir / "fig01_confidence_tiers_and_scores.png")
    L.write_tsv(pd.DataFrame({"confidence_tier": tiers, "n_genes": vals}),
                outdir / "fig01_confidence_tier_counts.tsv")
    hist, edges = np.histogram(score, bins=np.linspace(0, 1, 26))
    L.write_tsv(pd.DataFrame({"score_bin_low": edges[:-1], "score_bin_high": edges[1:],
                              "n_genes": hist}),
                outdir / "fig01_confidence_score_histogram.tsv")


# ---------------------------------------------------------------------------
# fig02 -- preliminary -> operon-adjusted -> final confidence
# ---------------------------------------------------------------------------
def fig02(genes: pd.DataFrame, outdir: Path, org_label: str) -> None:
    """Draws preliminary vs final confidence, the operon-correction histogram and before/after means.

    The final score equals the operon-adjusted score, so the correction is final - preliminary.
    """
    d = genes.copy()
    prelim = pd.to_numeric(d["preliminary_confidence_c1_c4"], errors="coerce")
    adj = pd.to_numeric(d["final_confidence_operon_context"], errors="coerce")
    final = pd.to_numeric(d["confidence_score"], errors="coerce")
    mask = prelim.notna() & final.notna()
    corr = final - prelim
    up = mask & (corr > 1e-9)
    dn = mask & (corr < -1e-9)
    eq = mask & (~up) & (~dn)
    r = float(np.corrcoef(prelim[mask], final[mask])[0, 1]) if mask.sum() > 1 else float("nan")

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.2, 5.4))

    # (a) preliminary vs final -- olive=raised, dark red=lowered, hollow black=unchanged
    axA.plot([0, 1], [0, 1], color="black", lw=1.3, ls="--", zorder=1)
    axA.scatter(prelim[eq], final[eq], s=13, facecolor="none", edgecolor="black",
                linewidths=0.5, label="unchanged", zorder=2)
    axA.scatter(prelim[up], final[up], s=16, facecolor=L.OLIVE, edgecolor="black",
                linewidths=0.4, alpha=0.9, label="raised by operon context", zorder=3)
    axA.scatter(prelim[dn], final[dn], s=16, facecolor=L.DARKRED, edgecolor="black",
                linewidths=0.4, alpha=0.9, label="lowered by operon context", zorder=4)
    axA.text(0.04, 0.965, f"r = {r:.2f}", transform=axA.transAxes, va="top", ha="left",
             fontsize=12.5, fontweight="bold", color=L.INK)
    axA.set_xlabel("preliminary confidence (C1–C4)")
    axA.set_ylabel("final confidence score")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.legend(loc="lower right", fontsize=8, frameon=True, facecolor="white",
               edgecolor="#bbbbbb", framealpha=0.95)
    L.set_title(axA, "Preliminary vs final")
    L.panel_letter(axA, "a")

    # (b) distribution of the operon correction
    cc = corr[mask].values
    edges = np.linspace(-0.5, 0.5, 41)
    counts, _ = np.histogram(cc, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    cols = [L.OLIVE if c > 1e-9 else (L.DARKRED if c < -1e-9 else "#888888") for c in centers]
    axB.bar(centers, counts, width=(edges[1] - edges[0]) * 0.9, color=cols,
            edgecolor="black", linewidth=0.3)
    axB.axvline(0, color=L.INK, lw=1.3, ls="--")
    axB.set_yscale("log")
    axB.set_xlabel("operon correction  (final − preliminary)")
    axB.set_ylabel("number of genes (log)")
    axB.set_xlim(-0.5, 0.5)
    L.set_title(axB, "How much operon context changed the score")
    L.panel_letter(axB, "b")

    # (c) mean before vs after
    means = [float(prelim[mask].mean()), float(final[mask].mean())]
    axC.plot([0, 1], means, "-o", color=L.BLUE, lw=2.4, markersize=12)
    for i, v in enumerate(means):
        axC.text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=11, color=L.INK)
    axC.set_xticks([0, 1])
    axC.set_xticklabels(["preliminary\n(C1–C4)", "final\n(after operon context)"], fontsize=10)
    axC.set_ylabel("mean confidence")
    axC.set_xlim(-0.35, 1.35)
    axC.set_ylim(min(means) - 0.05, max(means) + 0.05)
    axC.set_xlabel(f"{int(up.sum()):,} genes raised · {int(dn.sum()):,} lowered")
    L.set_title(axC, "Mean confidence before vs after")
    L.panel_letter(axC, "c")

    L.finish(fig, "How operon context changes confidence", organism=org_label, top=0.9)
    L.draw_sources_footer(fig, _RUN_ROOT, L.organism_source_lines(_RUN_ROOT, org_label))
    L.savefig(fig, outdir / "fig02_confidence_stages.png")
    out = pd.DataFrame({
        "feature_id": d["feature_id"], "preliminary_c1_c4": prelim,
        "after_operon_context": adj, "final_score": final,
        "operon_correction": corr,
        "operon_context_effect": np.where(up, "raised", np.where(dn, "lowered", "unchanged")),
    })
    L.write_tsv(out, outdir / "fig02_confidence_stages_per_gene.tsv")
    L.write_tsv(pd.DataFrame({"stage": ["preliminary_c1_c4", "after_operon_context", "final_score"],
                              "mean_confidence": [float(prelim[mask].mean()),
                                                  float(adj[mask].mean()),
                                                  float(final[mask].mean())],
                              "n_genes": [int(mask.sum())] * 3}),
                outdir / "fig02_confidence_stage_means.tsv")


# ---------------------------------------------------------------------------
# fig03 -- operon-context score and pangenome breadth by operon size
# ---------------------------------------------------------------------------
def _box_by_bin(ax, df, value_col, color, ylabel):
    """Draws box plots of value_col by operon size bin via reportfig_lib.box_by_bin."""
    return L.box_by_bin(ax, df, "size_bin", value_col, L.SIZE_BIN_ORDER, color,
                        ylabel, "operon size (number of member genes)")


def fig03(genes: pd.DataFrame, operons: pd.DataFrame, recurrence: dict,
          outdir: Path, org_label: str) -> None:
    """Draws C3 by operon size and the share of operons recurring in >= 2 pangenome organisms."""
    d = genes[genes["in_operon"]].copy()
    d["size_bin"] = d["operon_member_count"].map(L.size_bin)

    ops = operons.copy()
    ops["size_bin"] = ops["size"].map(L.size_bin)
    ops["pangenome_organisms"] = ops["members_in_order"].map(
        lambda m: recurrence.get(m, {}).get("organism_count", np.nan))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.0, 5.8))
    labA, nA, medA, meanA = _box_by_bin(axA, d, "c3_score", L.GREEN,
                                        "operon-context score (C3)")
    axA.set_ylim(0, 1.02)
    axA.grid(False)
    L.set_title(axA, "Operon-context score by operon size")
    L.panel_letter(axA, "a")

    # (b) share of this organism's operons that also occur in >=2 pangenome
    # organisms, by operon size (matched operons only).
    matched = ops[ops["pangenome_organisms"].notna()]
    labB, pctB, nB = [], [], []
    for b in L.SIZE_BIN_ORDER:
        sub = matched[matched["size_bin"] == b]
        if len(sub):
            shared = int((sub["pangenome_organisms"] >= 2).sum())
            labB.append(b); nB.append(len(sub))
            pctB.append(100.0 * shared / len(sub))
    bars = axB.bar(range(len(labB)), pctB,
                   color=[L.CATEGORICAL[i % len(L.CATEGORICAL)] for i in range(len(labB))],
                   edgecolor="white",
                   linewidth=1.2, width=0.72)
    for bar, p in zip(bars, pctB):
        axB.text(bar.get_x() + bar.get_width() / 2, p, f"{p:.0f}%", ha="center",
                 va="bottom", fontsize=9.5, color=L.INK)
    axB.set_xticks(range(len(labB)))
    axB.set_xticklabels([f"{l}\n(n={n:,})" for l, n in zip(labB, nB)], fontsize=9)
    axB.set_ylabel("share of operons also seen in\n≥2 pangenome organisms (%)")
    axB.set_xlabel("operon size (number of member genes)")
    axB.set_ylim(0, max(pctB + [1]) * 1.2)
    axB.grid(False)
    L.set_title(axB, "How often each operon recurs across the pangenome")
    L.panel_letter(axB, "b")

    L.finish(fig, "Operon context vs operon size", organism=org_label)
    L.draw_sources_footer(fig, _RUN_ROOT,
                          L.organism_source_lines(_RUN_ROOT, org_label, operon_db=_OPERON_DB))
    L.savefig(fig, outdir / "fig03_operon_context_by_size.png")
    L.write_tsv(pd.DataFrame({"operon_size_bin": labA, "n_genes": nA,
                              "median_c3_score": medA, "mean_c3_score": meanA}),
                outdir / "fig03_c3_by_operon_size.tsv")
    L.write_tsv(pd.DataFrame({"operon_size_bin": labB, "n_operons": nB,
                              "pct_shared_in_2plus_organisms": pctB}),
                outdir / "fig03_recurrence_by_operon_size.tsv")


# ---------------------------------------------------------------------------
# fig04 -- distribution of each confidence component (C1..C4)
# ---------------------------------------------------------------------------
def fig04(genes: pd.DataFrame, outdir: Path, org_label: str) -> None:
    """Draws a histogram with median for each of C1-C4."""
    comps = [("C1", "c1_score"), ("C2", "c2_score_from_operon_probability"),
             ("C3", "c3_score"), ("C4", "c4_score")]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.0))
    axes = axes.ravel()
    rows = []
    for ax, (key, col) in zip(axes, comps):
        v = pd.to_numeric(genes[col], errors="coerce").dropna()
        ax.hist(v, bins=np.linspace(0, 1, 26), color=L.COMPONENT_COLOR[key],
                edgecolor="white", linewidth=0.5)
        med = float(v.median())
        ax.axvline(med, color=L.INK, lw=1.8, ls="--", label=f"median = {med:.2f}")
        ax.set_xlim(0, 1)
        ax.set_xlabel("score")
        ax.set_ylabel("number of genes")
        ax.legend(loc="best", fontsize=9, frameon=True, facecolor="white",
                  edgecolor="#bbbbbb", framealpha=0.92)
        L.set_title(ax, L.COMPONENT_LABEL[key])
        L.panel_letter(ax, chr(ord("a") + comps.index((key, col))))
        rows.append({"component": key, "n_genes": len(v), "mean": float(v.mean()),
                     "median": med, "min": float(v.min()), "max": float(v.max())})
    L.finish(fig, "Confidence components", organism=org_label)
    L.draw_sources_footer(fig, _RUN_ROOT, L.organism_source_lines(_RUN_ROOT, org_label))
    L.savefig(fig, outdir / "fig04_component_distributions.png")
    L.write_tsv(pd.DataFrame(rows), outdir / "fig04_component_summary.tsv")


# ---------------------------------------------------------------------------
# fig05 -- one real operon neighbourhood (gene track)
# ---------------------------------------------------------------------------
def _pick_window(genes: pd.DataFrame, span: int = 11) -> pd.DataFrame:
    """Returns a window of genes on the largest contig centred on a boundary between two operons at least 100 bp apart."""
    g = genes.dropna(subset=["gene_start"]).copy()
    if g.empty:
        return g
    contig = g["contig"].value_counts().index[0]
    sub = g[g["contig"] == contig].sort_values("gene_start").reset_index(drop=True)

    def oid(i):
        """Returns row i's operon id when it is in a real operon, else None."""
        v = sub.loc[i, "operon_id"]
        return v if (v and v != L._NOT_IN_OPERON and bool(sub.loc[i, "in_operon"])) else None

    sizes = sub["operon_id"].value_counts()
    best = None
    for i in range(len(sub) - 1):
        a, b = oid(i), oid(i + 1)
        if a and b and a != b:
            try:
                gap = int(sub.loc[i + 1, "gene_start"]) - int(sub.loc[i, "gene_end"]) - 1
            except (TypeError, ValueError):
                continue
            if gap >= 100 and sizes.get(a, 0) >= 2 and sizes.get(b, 0) >= 2:
                best = i
                break
    if best is None:
        return sub.iloc[: min(span, len(sub))]
    lo = max(0, best - span // 2)
    hi = min(len(sub), best + span // 2 + 1)
    return sub.iloc[lo:hi]


def fig05(genes: pd.DataFrame, outdir: Path, org_label: str) -> None:
    """Draws one real genome window with its operons via reportfig_lib.render_operon_page."""
    from matplotlib.patches import Patch
    win = _pick_window(genes)
    members = [{"start": r.gene_start, "end": r.gene_end, "strand": r.RAST_strand,
                "label": r.best_consensus_product_descriptor, "operon_id": r.operon_id,
                "source": getattr(r, "product_descriptor_source", ""),
                "needs_review": getattr(r, "needs_review", ""),
                "review_reason": getattr(r, "needs_review_reason", ""),
                "confidence": r.confidence_score, "c3": r.c3_score,
                "c2": getattr(r, "c2_score_from_operon_probability", None),
                "c1": getattr(r, "c1_score", None), "c4": getattr(r, "c4_score", None),
                "preliminary": r.preliminary_confidence_c1_c4,
                "operon_adjusted": r.final_confidence_operon_context,
                "operon_adjusted_hybrid": getattr(r, "final_confidence_operon_context_hybrid", None),
                "c3_hybrid": getattr(r, "c3_score_operon_context_hybrid", None),
                "tier": getattr(r, "confidence_tier", "")}
               for r in win.itertuples()]
    badges = {m["operon_id"]: m["operon_id"]
              for m in members if L.is_operon(m.get("operon_id"))}
    handles = [Patch(facecolor="#eef1f6", edgecolor="#c3ccd6", label="genes of one operon"),
               Patch(facecolor="#9aa0a6", edgecolor="#222222",
                     label="gene  (arrow points 5'→3' along its coding strand)")]
    # Shared renderer, same formatting as the atlas.
    span = f"{int(win.gene_start.min()):,}–{int(win.gene_end.max()):,} bp"
    L.render_operon_page(outdir / "fig05_operon_neighbourhood.png",
                         [{"members": members, "heading": span,
                           "detail": "neighbouring genes, each operon's members side by side"}],
                         org_label=org_label, suptitle="A real stretch of the genome",
                         run_root=_RUN_ROOT, note=L.OPERON_CORRECTION_NOTE,
                         footer_sources=L.organism_source_lines(_RUN_ROOT, org_label, coords=True),
                         per_row=10, min_span=10, badges=badges, legend_handles=handles)
    out = win[["feature_id", "operon_id", "operon_member_count", "gene_start",
               "gene_end", "RAST_strand", "best_consensus_product_descriptor",
               "confidence_score"]].copy()
    L.write_tsv(out, outdir / "fig05_operon_neighbourhood.tsv")


# ---------------------------------------------------------------------------
# fig06 / fig07 / fig08 -- galleries of the most-reproduced, most-unique and
# penalised operons: gene track plus full descriptors and confidence.
# ---------------------------------------------------------------------------
_GALLERY_N = 10
_GALLERY_CAP = 6   # members drawn per operon row


def _operon_gallery(ops_ranked: pd.DataFrame, outdir: Path, org_label: str,
                    fname: str, tsv: str, suptitle: str, mode: str = "reproduced",
                    note: str | None = None, pool_n: int | None = None) -> None:
    """Renders the top _GALLERY_N operons as one gallery page and writes its TSV."""
    top = ops_ranked.head(_GALLERY_N).reset_index(drop=True)
    if len(top) == 0:
        return
    blocks, rows = [], []
    for _, r in top.iterrows():
        size = int(r["size"])
        k = int(r.get("pangenome_organisms", 0) or 0)
        members = L.operon_to_members(r)                 # ALL members (no cap)
        ofn = f"{k} of {pool_n}" if pool_n else str(k)   # recurrence as k / pool size
        if mode == "unique":
            where = "seen in no other genome"
        elif mode == "penalized":
            # recurrence explains the penalty (low conservation, C3)
            where = f"lowered by weak conservation — operon seen in {ofn} genomes"
        else:
            where = f"conserved in {ofn} genomes" if k != 1 else f"in 1 of {pool_n or '?'} genomes"
        blocks.append({"members": members, "heading": r["operon_id"],
                       "detail": f"{size}-gene operon  |  {where}"})
        rows.append({"operon_id": r["operon_id"], "size": size,
                     "pangenome_organisms": k,
                     "operon_database_pool_organisms": pool_n,
                     "pangenome_occurrences": int(r.get("pangenome_occurrences", 0) or 0),
                     "operon_confidence_penalty": round(float(r.get("penalty", 0.0) or 0.0), 4),
                     "members_in_order": r["members_in_order"]})
    # Shared renderer, same formatting as the atlas.
    L.render_operon_page(outdir / fname, blocks, org_label=org_label, suptitle=suptitle,
                         run_root=_RUN_ROOT, note=note or L.OPERON_CORRECTION_NOTE,
                         footer_sources=L.organism_source_lines(_RUN_ROOT, org_label,
                                                                coords=True, operon_db=_OPERON_DB),
                         per_row=10, min_span=10)
    L.write_tsv(pd.DataFrame(rows), outdir / tsv)


def _operon_penalty(r) -> float:
    """Returns the summed downward correction Σ max(0, preliminary - final) over an operon's members."""
    tot = 0.0
    for p, f in zip(r.get("preliminaries", []) or [], r.get("confidences", []) or []):
        try:
            d = float(p) - float(f)
            if d > 0 and d == d:      # d==d guards NaN
                tot += d
        except (TypeError, ValueError):
            pass
    return tot


def fig06_07(operons: pd.DataFrame, recurrence: dict, outdir: Path, org_label: str,
             pool_n: int | None = None) -> None:
    """Draws the reproduced (fig06), unique (fig07) and penalised (fig08) operon galleries."""
    ops = operons.copy()
    ops["pangenome_organisms"] = ops["members_in_order"].map(
        lambda m: recurrence.get(m, {}).get("organism_count", 0))
    ops["pangenome_occurrences"] = ops["members_in_order"].map(
        lambda m: recurrence.get(m, {}).get("label_frequency", 0))
    # show operons of named-function genes only (same informative gate as scoring)
    ops["all_named"] = ops["members_in_order"].map(
        lambda m: L.operon_members_informative(m))
    ops["mostly_named"] = ops["members_in_order"].map(
        lambda m: L.operon_members_informative(m, min_informative=2))
    ops["penalty"] = ops.apply(_operon_penalty, axis=1)

    matched = ops[(ops["pangenome_organisms"] > 0) & ops["all_named"]]
    reproduced = matched.sort_values(
        ["pangenome_organisms", "pangenome_occurrences", "size"], ascending=False)
    _operon_gallery(reproduced, outdir, org_label,
                    "fig06_top_reproduced_operons.png",
                    "fig06_top_reproduced_operons.tsv",
                    "Named-function operons most often seen across the pangenome",
                    mode="reproduced", pool_n=pool_n)

    unique = ops[(ops["pangenome_organisms"] <= 1) & ops["mostly_named"]].copy()
    unique["pangenome_organisms"] = unique["pangenome_organisms"].replace(0, 1)
    unique = unique.sort_values("size", ascending=False)
    _operon_gallery(unique, outdir, org_label,
                    "fig07_most_unique_operons.png",
                    "fig07_most_unique_operons.tsv",
                    "Named-function operons seen in no other pangenome genome",
                    mode="unique", pool_n=pool_n)

    # fig08 -- operons where operon context PENALISED (lowered) member confidence
    penalized = ops[(ops["penalty"] > 0.01) & ops["mostly_named"]].sort_values(
        "penalty", ascending=False)
    _operon_gallery(penalized, outdir, org_label,
                    "fig08_penalized_operons.png",
                    "fig08_penalized_operons.tsv",
                    "Operons where operon context LOWERED gene confidence",
                    mode="penalized", note=L.OPERON_PENALTY_NOTE, pool_n=pool_n)


# ---------------------------------------------------------------------------
def main() -> None:
    """Loads the organism and its leave-one-out OCC pool, draws every figure and exits 2 if any fails."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--operon-db", default=str(L.DEFAULT_OPERON_DB))
    ap.add_argument("--output-dir", default=None,
                    help="default: <run>/<organism>/scoring/figures")
    args = ap.parse_args()

    run_root = Path(args.run_root)
    outdir = Path(args.output_dir) if args.output_dir else \
        run_root / args.organism / "scoring" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    global _RUN_ROOT, _OPERON_DB
    _RUN_ROOT, _OPERON_DB = run_root, Path(args.operon_db)
    L.apply_style()
    org_label = args.organism
    genes = L.load_organism_genes(run_root, args.organism)
    operons = L.build_operons(genes)
    # Recurrence and the pool caption use the OCC reference's organisms,
    # falling back to this run's organisms when the OCC is unreadable.
    restrict = L.load_occ_organisms() or set(L.discover_organisms(run_root))
    # Leave-one-out: the pool excludes this organism, as in its C3 scoring.
    restrict = set(restrict)
    restrict.discard(org_label)
    recurrence = L.load_operon_recurrence(Path(args.operon_db), restrict_to=restrict)
    # Provenance line: pool size and tallies from the OCC genome-stats sidecar.
    pool_list = sorted(restrict)
    pstats = L.aggregate_pool_stats(L.load_pool_stats(), pool_list)
    L.set_provenance(L.provenance_text(len(pool_list), pstats, leave_one_out=True))

    jobs = [
        ("fig01", lambda: fig01(genes, outdir, org_label)),
        ("fig02", lambda: fig02(genes, outdir, org_label)),
        ("fig03", lambda: fig03(genes, operons, recurrence, outdir, org_label)),
        ("fig04", lambda: fig04(genes, outdir, org_label)),
        ("fig05", lambda: fig05(genes, outdir, org_label)),
        ("fig06_07", lambda: fig06_07(operons, recurrence, outdir, org_label,
                                      pool_n=len(pool_list))),
    ]
    failures = []
    for name, fn in jobs:
        try:
            fn()
        except Exception:
            failures.append(name)
            print(f"[make_organism_report] {name} FAILED:\n{traceback.format_exc()}",
                  file=sys.stderr)
    L.write_sources_manifest(outdir, run_root, [args.organism], Path(args.operon_db))
    print(f"[make_organism_report] {args.organism}: {len(jobs) - len(failures)}/{len(jobs)} "
          f"figure groups ok → {outdir}")
    if failures:
        print(f"[make_organism_report] FAILURES: {failures}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
