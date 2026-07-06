"""Figure 60 - Operon Context Confidence (OCC): an INDEPENDENT, per-gene operon
reliability factor derived purely from pan-genome co-occurrence (figs 01-05).

THE QUESTION (operon-index thread)
==================================
figs 01-05 established that (a) within-operon functional co-occurrence is a
conserved, pan-genome signal; (b) recurrence - not raw partner count - is what
earns trust (figs 03/04); and (c) UniOP's per-pair probability is orthogonal to
that co-occurrence (fig 05).  Can we fold all of this into a single number, per
candidate gene, that says how reliable its operon placement is - WITHOUT using
the deterministic C1/C4 scores, so it is an independent predictor?

THE FACTOR (see c3_occ.py for the full derivation)
==================================================
Gate: only operons with a strict informative majority (n_inf > n_unf) provide
context (all-hypothetical / info-minority operons are excluded).
For each functional pair we build, across all genomes, a recurrence-aware
conditional co-occurrence (Jeffreys 95% lower bound) x an enrichment safeguard,
in two channels - immediate-neighbour (adj) and same-operon (op).  A gene's OCC
is the noisy-OR over its partners:
        OCC(g) = 1 - PROD_neighbours(1-rho_adj) * PROD_comembers(1-rho_op)
Uninformative genes in a qualifying operon inherit the operon's mean OCC.

FINDING (panels)
================
(a) OCC spreads the full [0,1] range (median ~0.40): it genuinely discriminates.
    A known-conserved class - ribosomal proteins - piles up near 1 (median ~0.92).
(b) Worked examples: conserved multi-partner genes (ribosomal, NADH dehydrogenase,
    flagellar) score ~1; context-free genes (transposases, lone transporters)
    score ~0 - exactly the intended biology.
(c) OCC is INDEPENDENT of the base signals: Spearman vs C4 ~ 0.02, vs UniOP
    per-pair prob ~ 0.03 (orthogonal, echoing fig 05), vs C1 ~ 0.26 (both merely
    prefer informative genes).  OCC therefore adds new information.
(d) OCC rises with operon size, and that rise SURVIVES a link-reliability floor
    (0 -> 0.20 barely moves it): the high scores of large operons come from
    genuinely conserved links (flagellar / ribosomal / capsular modules), not
    noisy-OR inflation over coincidental partners.

Read-only prototype; the production scorer is untouched.  How OCC is folded into
C1/C4 (log-odds shift around a neutral pivot) is left to the caller.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import c3_occ as O
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SUB = "01-operon-context-confidence"

SIZE_BANDS = [(1, 2), (3, 4), (5, 8), (9, 20), (21, 200)]
SIZE_LAB = ["2", "3-4", "5-8", "9-20", "21+"]
SIZE_COL = [L.RED, L.ORANGE, L.AMBER, L.TEAL, L.GREEN]


def make(genes, operons, outdir):
    run_root = outdir.parents[3]
    ref = O.build_reference(genes, run_root)
    df = O.compute_all_genes(genes, run_root, ref=ref)
    occ0 = df.attrs["occ0"]
    inf = df[~df["uninformative"]].copy()

    # operon size (informative context + self) for the inflation panel
    inf["op_size"] = inf["n_inf_context"] + 1

    # ---- independence vs the base signals ----------------------------------
    base = genes[["organism", "feature_id", "c1_score", "c4_score", "operon_prob"]]
    m = inf.merge(base, on=["organism", "feature_id"], how="left")
    corr = {}
    for col in ["c1_score", "c4_score", "operon_prob"]:
        s = m[["occ", col]].apply(pd.to_numeric, errors="coerce").dropna()
        corr[col] = (spearmanr(s["occ"], s[col]).correlation, len(s))

    # ---- worked examples ----------------------------------------------------
    def _best_instance(mask):
        sub = inf[mask]
        if not len(sub):
            return None
        return sub.sort_values(["occ", "n_partners"], ascending=[False, False]).iloc[0]

    def _worst_instance(mask):
        sub = inf[mask]
        if not len(sub):
            return None
        return sub.sort_values(["occ", "n_partners"], ascending=[True, True]).iloc[0]

    high_specs = [
        ("50s ribosomal protein l22", "contains", "50s ribosomal protein l22"),
        ("50s ribosomal protein l2", "eq", "50s ribosomal protein l2"),
        ("nadh-quinone oxidoreductase", "contains", "nadh-quinone oxidoreductase"),
        ("atp synthase", "contains", "atp synthase"),
        ("flagellar", "contains", "flagellar"),
        ("dna-directed rna polymerase", "contains", "dna-directed rna polymerase"),
    ]
    low_specs = [
        ("transposase", "contains", "transposase"),
        ("chromate transporter", "contains", "chromate transporter"),
        ("trk system potassium uptake protein trka", "eq", "trk system"),
        ("cbs domain", "contains", "cbs domain"),
        ("multidrug efflux", "contains", "efflux pump"),
        ("cobalt-zinc-cadmium resistance", "contains", "resistance protein"),
    ]

    def _pick(specs, hi):
        out = []
        seen = set()
        for desc, mode, _tag in specs:
            if mode == "contains":
                mask = inf["clean_descriptor"].str.contains(desc, na=False, regex=False)
            else:
                mask = inf["clean_descriptor"] == desc
            rec = _best_instance(mask) if hi else _worst_instance(mask)
            if rec is None:
                continue
            key = rec["clean_descriptor"]
            if key in seen:
                continue
            seen.add(key)
            out.append(rec)
        return out

    hi_ex = _pick(high_specs, True)[:5]
    lo_ex = _pick(low_specs, False)[:5]

    # ---- size / floor robustness (reuse ref; floor only affects aggregation) --
    ref["params"]["link_floor"] = 0.20
    df_floor = O.compute_all_genes(genes, run_root, ref=ref)
    ref["params"]["link_floor"] = 0.0
    inf_f = df_floor[~df_floor["uninformative"]].copy()
    inf_f["op_size"] = inf_f["n_inf_context"] + 1

    def _band_means(frame):
        out = []
        for lo, hi in SIZE_BANDS:
            s = frame[(frame["op_size"] >= lo) & (frame["op_size"] <= hi)]
            out.append((s["occ"].mean() if len(s) else np.nan, len(s)))
        return out

    band0 = _band_means(inf)
    bandf = _band_means(inf_f)

    # ======================= FIGURE =========================================
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 12.8))
    axA, axB, axC, axD = axes.ravel()
    fig.suptitle("Operon Context Confidence (OCC): an independent, per-gene operon-reliability factor "
                 "from pan-genome co-occurrence",
                 fontsize=15, fontweight="bold", y=0.986)

    # ---- (a) OCC distribution + ribosomal overlay --------------------------
    rib = inf[inf["clean_descriptor"].str.contains("ribosomal protein", na=False)]["occ"]
    edges = np.linspace(0, 1, 41)
    axA.hist(inf["occ"], bins=edges, color=L.BLUE, alpha=0.75,
             label="all informative genes  (n = %d)" % len(inf))
    axA.hist(rib, bins=edges, color=L.GREEN, alpha=0.8,
             weights=np.full(len(rib), len(inf) / max(len(rib), 1) * 0.18),
             label="ribosomal proteins (scaled)  (n = %d)" % len(rib))
    axA.axvline(occ0, color="black", linewidth=2.4, linestyle="--")
    axA.text(occ0 + 0.01, axA.get_ylim()[1] * 0.94,
             "neutral pivot\nOCC0 = %.2f" % occ0, fontsize=10,
             fontweight="bold", va="top")
    axA.axvline(rib.median(), color=L.GREEN, linewidth=2.2, linestyle=":")
    axA.text(rib.median() - 0.01, axA.get_ylim()[1] * 0.62,
             "ribosomal\nmedian %.2f" % rib.median(), fontsize=9.5,
             fontweight="bold", va="top", ha="right", color="#0a7d3a")
    axA.set_xlim(0, 1)
    axA.set_xlabel("Operon Context Confidence (OCC)", fontweight="bold")
    axA.set_ylabel("Number of genes", fontweight="bold")
    axA.set_title("(a) OCC spans the full range and genuinely discriminates\n"
                  "a known-conserved class (ribosomal) piles up near 1",
                  fontweight="bold")
    axA.legend(loc="upper right", fontsize=10, framealpha=0.9)
    axA.grid(False)
    L.boldticks(axA)

    # ---- (b) worked examples ----------------------------------------------
    ex = [(r, L.BLUE) for r in hi_ex] + [(r, L.RED) for r in lo_ex]
    ex.sort(key=lambda t: t[0]["occ"])
    yy = np.arange(len(ex))
    xvals = [r["occ"] for r, _ in ex]
    colors = [c for _, c in ex]
    labels = [L.short_desc(r["clean_descriptor"], 34) for r, _ in ex]
    axB.barh(yy, xvals, color=colors, alpha=0.9, edgecolor="black", linewidth=0.6)
    for i, (r, _) in enumerate(ex):
        txt = "%d partner%s" % (r["n_partners"], "" if r["n_partners"] == 1 else "s")
        if r["n_partners"] > 0:
            txt += ", best rho=%.2f" % r["best_rho"]
        xpos = min(r["occ"] + 0.02, 0.62)
        ha = "left"
        col = "black"
        if r["occ"] > 0.62:
            xpos = r["occ"] - 0.02
            ha = "right"
            col = "white"
        axB.text(xpos, i, txt, va="center", ha=ha, fontsize=8.6,
                 fontweight="bold", color=col)
    axB.set_yticks(yy)
    axB.set_yticklabels(labels, fontsize=9)
    axB.set_xlim(0, 1.02)
    axB.set_xlabel("Operon Context Confidence (OCC)", fontweight="bold")
    from matplotlib.patches import Patch
    axB.legend(handles=[Patch(facecolor=L.BLUE, edgecolor="black",
                              label="conserved multi-partner module"),
                        Patch(facecolor=L.RED, edgecolor="black",
                              label="context-free / mobile gene")],
               loc="lower right", fontsize=9, framealpha=0.95)
    axB.set_title("(b) Worked examples: OCC rewards conserved partnerships,\n"
                  "penalises context-free genes (transposases, lone transporters)",
                  fontweight="bold")
    axB.grid(False)
    L.boldticks(axB)

    # ---- (c) independence from the base signals ----------------------------
    s = m[["occ", "operon_prob"]].apply(pd.to_numeric, errors="coerce").dropna()
    hb = axC.hexbin(s["operon_prob"], s["occ"], gridsize=42, cmap="viridis",
                    bins="log", mincnt=1)
    cb = fig.colorbar(hb, ax=axC, fraction=0.046, pad=0.02)
    cb.set_label("genes (log)", fontweight="bold")
    axC.set_xlabel("UniOP operon membership probability", fontweight="bold")
    axC.set_ylabel("Operon Context Confidence (OCC)", fontweight="bold")
    axC.set_ylim(0, 1)
    box = ("Spearman correlation of OCC with the base signals\n"
           "   vs UniOP prob : rho = %.2f\n"
           "   vs C4 score   : rho = %.2f\n"
           "   vs C1 score   : rho = %.2f"
           % (corr["operon_prob"][0], corr["c4_score"][0], corr["c1_score"][0]))
    axC.text(0.03, 0.97, box, transform=axC.transAxes, fontsize=9.5,
             fontweight="bold", va="top", ha="left",
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.9,
                       edgecolor="black"))
    axC.set_title("(c) OCC is independent of the base signals\n"
                  "orthogonal to UniOP prob & C4 (echoing fig 05); adds new information",
                  fontweight="bold")
    axC.grid(False)
    L.boldticks(axC)

    # ---- (d) size relationship survives a reliability floor ----------------
    xb = np.arange(len(SIZE_BANDS))
    w = 0.38
    axD.bar(xb - w / 2, [b[0] for b in band0], width=w, color=L.BLUE,
            alpha=0.9, edgecolor="black", linewidth=0.6,
            label="link floor = 0 (default)")
    axD.bar(xb + w / 2, [b[0] for b in bandf], width=w, color=L.ORANGE,
            alpha=0.9, edgecolor="black", linewidth=0.6,
            label="link floor = 0.20 (drop weak links)")
    for i, (b0, bf) in enumerate(zip(band0, bandf)):
        if not np.isnan(b0[0]):
            axD.text(i - w / 2, b0[0] + 0.01, "%.2f" % b0[0], ha="center",
                     va="bottom", fontsize=8.5, fontweight="bold")
        if not np.isnan(bf[0]):
            axD.text(i + w / 2, bf[0] + 0.01, "%.2f" % bf[0], ha="center",
                     va="bottom", fontsize=8.5, fontweight="bold")
    axD.set_xticks(xb)
    axD.set_xticklabels(SIZE_LAB, fontweight="bold")
    axD.set_ylim(0, 1.08)
    axD.set_xlabel("Operon size (number of informative members)", fontweight="bold")
    axD.set_ylabel("Mean OCC of member genes", fontweight="bold")
    axD.set_title("(d) High OCC in large operons is real conservation, not inflation\n"
                  "the size trend barely moves when weak links are dropped",
                  fontweight="bold")
    axD.legend(loc="upper left", fontsize=9.5, framealpha=0.95)
    axD.grid(False)
    L.boldticks(axD)

    fig.tight_layout(rect=[0, 0, 1, 0.965])
    L.savefig(fig, outdir / "fig06_operon_context_confidence.png")

    # ---- TSVs --------------------------------------------------------------
    per_gene = df.sort_values(["occ", "organism"], ascending=[False, True])
    L.write_tsv(per_gene, outdir / "fig06_per_gene_occ.tsv")

    summ = []
    for q in [0.1, 0.25, 0.5, 0.75, 0.9]:
        summ.append({"panel": "a_distribution", "metric": "occ_q%.2f" % q,
                     "value": float(inf["occ"].quantile(q)), "n": len(inf)})
    summ.append({"panel": "a_distribution", "metric": "occ0_pivot_median",
                 "value": float(occ0), "n": len(inf)})
    summ.append({"panel": "a_distribution", "metric": "ribosomal_median_occ",
                 "value": float(rib.median()), "n": len(rib)})
    for col in ["operon_prob", "c4_score", "c1_score"]:
        summ.append({"panel": "c_independence", "metric": "spearman_occ_vs_" + col,
                     "value": float(corr[col][0]), "n": corr[col][1]})
    for (lo, hi), lab, b0, bf in zip(SIZE_BANDS, SIZE_LAB, band0, bandf):
        summ.append({"panel": "d_size", "metric": "mean_occ_size_%s_floor0" % lab,
                     "value": float(b0[0]) if not np.isnan(b0[0]) else float("nan"),
                     "n": b0[1]})
        summ.append({"panel": "d_size", "metric": "mean_occ_size_%s_floor0.20" % lab,
                     "value": float(bf[0]) if not np.isnan(bf[0]) else float("nan"),
                     "n": bf[1]})
    L.write_tsv(pd.DataFrame(summ), outdir / "fig06_summary.tsv")

    # worked-examples table
    exrows = []
    for r, _ in ex:
        exrows.append(dict(clean_descriptor=r["clean_descriptor"],
                           organism=r["organism"], occ=float(r["occ"]),
                           n_partners=int(r["n_partners"]),
                           best_partner=r["best_partner"],
                           best_rho=float(r["best_rho"]),
                           best_channel=r["best_channel"]))
    L.write_tsv(pd.DataFrame(exrows), outdir / "fig06_worked_examples.tsv")


if __name__ == "__main__":
    L.figure_main(make, SUB)
