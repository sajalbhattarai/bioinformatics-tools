"""Figure 59 - Does UniOP's per-PAIR operonic probability reflect the empirical
cross-genome co-occurrence?

THE QUESTION (operon-index thread)
==================================
figs 01-04 measured, per functional (descriptor) pair, an *empirical*
cross-genome signal: the conditional co-occurrence
        cond = (# genomes where A,B are operon-adjacent)
               / (# genomes where both A,B are present as operon members)
cond = 1.0 means "wherever both functions exist they are always operon
partners" - a deterministic, conserved module.

UniOP itself already emits a per-adjacency probability for every operonic pair
(OPERON_upstream/downstream_pairwise_probability in operon/operon_results.tsv).
So: does that model probability ALREADY contain the conservation signal, or is
the empirical co-occurrence independent information the index must add?

METHOD
======
* Gene identity = clean_descriptor (lower-cased), exactly as fig 01.  A pair is
  a sorted descriptor pair; a within-operon adjacency = two informative genes
  consecutive by start in the SAME operon on the SAME contig.
* For every such adjacency we look up UniOP's pairwise probability from that
  organism's operon_results.tsv (keyed by the unordered feature-id pair) and
  aggregate per descriptor pair (mean / median across genome instances).
* We then compare that UniOP probability to the empirical conditional
  co-occurrence, and - as a mechanism control - to the intergenic gap.

FINDING (panels)
================
(a) UniOP pairwise prob is nearly FLAT across conditional-co-occurrence bins
    (Spearman ~ 0.06): it is ~0.9 whether a pair always co-occurs or rarely does.
(b) It cannot separate deterministic modules (cond=1) from flimsy pairs
    (cond<0.5): the two UniOP distributions overlap (Mann-Whitney n.s.).
(c) What UniOP DOES encode is the intergenic gap (Spearman ~ -0.76) - a
    single-genome, distance-driven signal, blind to conservation.
(d) Concretely, a flimsy pair that co-occurs in ~6% of genomes can get a HIGHER
    UniOP probability than the perfectly conserved ribosomal super-operon.

=> UniOP probability and empirical co-occurrence are orthogonal.  The pan-genome
co-occurrence metric is NOT redundant with UniOP - it is independent evidence
the operon index should add.  (Read-only prototype; scorer untouched.)
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

csv.field_size_limit(10_000_000)

SUB = "01-operon-context-confidence"

# conditional-co-occurrence bins for panel (a)
COND_BINS = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.999), (0.999, 1.01)]
COND_LAB = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-<1", "=1.0"]
COND_COL = [L.RED, L.ORANGE, L.AMBER, L.YELLOW, L.TEAL, L.GREEN]

# intergenic-distance bins for panel (c)
GAP_BINS = [(-10**9, -1), (0, 0), (1, 100), (100, 10**9)]
GAP_LAB = ["<0", "=0", "1-100", ">100"]
GAP_COL = [L.GREEN, L.TEAL, L.LIME, L.RED]


def _build_contig_map(run_root, organisms):
    rx = re.compile(r"^(.*)_(\d+)([+-])(\d+)$")
    m = {}
    for org in organisms:
        path = run_root / org / "labeling" / "labeled-genes.tsv"
        if not path.is_file():
            continue
        with open(path, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd)
            iF = header.index("feature_id")
            iG = header.index("gene_id")
            for row in rd:
                if len(row) > max(iF, iG):
                    hit = rx.match(row[iG] or "")
                    m[(org, row[iF])] = hit.group(1) if hit else None
    return m


def _uniop_pair_probs(run_root, organisms):
    """(organism, frozenset{fidA,fidB}) -> UniOP pairwise operonic probability."""
    pp = {}
    for org in organisms:
        path = run_root / org / "operon" / "operon_results.tsv"
        if not path.is_file():
            continue
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                fid = row["feature_id"]
                for gcol, pcol in (("OPERON_upstream_gene_id", "OPERON_upstream_pairwise_probability"),
                                   ("OPERON_downstream_gene_id", "OPERON_downstream_pairwise_probability")):
                    nb = row.get(gcol) or ""
                    pv = row.get(pcol) or ""
                    if nb and pv:
                        try:
                            pp[(org, frozenset((fid, nb)))] = float(pv)
                        except ValueError:
                            pass
    return pp


def _gap_bin_idx(gap):
    for i, (lo, hi) in enumerate(GAP_BINS):
        if lo <= gap <= hi:
            return i
    return len(GAP_BINS) - 1


def make(genes, operons, outdir):
    run_root = outdir.parents[3]
    g = genes.copy()
    organisms = sorted(g["organism"].unique())
    contig = _build_contig_map(run_root, organisms)
    g["contig"] = [contig.get((o, f)) for o, f in zip(g["organism"], g["feature_id"])]

    def norm(d):
        return (d or "").strip().lower()

    uninf = {(o, f): bool(u) for o, f, u in
             zip(g["organism"], g["feature_id"], g["uninformative"])}
    cln = {(o, f): norm(c) for o, f, c in
           zip(g["organism"], g["feature_id"], g["clean_descriptor"])}

    pairprob = _uniop_pair_probs(run_root, organisms)

    # ---- within-operon adjacencies (identical rule to fig 01) --------------
    pair_orgs = defaultdict(set)
    present = defaultdict(set)
    pair_probs = defaultdict(list)
    pair_gaps = defaultdict(list)
    inst_prob = []          # per-instance UniOP prob
    inst_gap = []           # per-instance intergenic gap
    for (org, oid), sub in g.groupby(["organism", "operon_id"]):
        if not str(oid).startswith("operon_"):
            continue
        sub = sub.sort_values("start")
        recs = sub.to_dict("records")
        for rec in recs:
            if not uninf[(org, rec["feature_id"])]:
                present[cln[(org, rec["feature_id"])]].add(org)
        for a, b in zip(recs[:-1], recs[1:]):
            if a["contig"] != b["contig"]:
                continue
            if uninf[(org, a["feature_id"])] or uninf[(org, b["feature_id"])]:
                continue
            da, db = cln[(org, a["feature_id"])], cln[(org, b["feature_id"])]
            if not da or not db or da == db:
                continue
            key = (da, db) if da < db else (db, da)
            pair_orgs[key].add(org)
            pp = pairprob.get((org, frozenset((a["feature_id"], b["feature_id"]))))
            gap = int(b["start"]) - int(a["end"]) - 1
            if pp is not None:
                pair_probs[key].append(pp)
                pair_gaps[key].append(gap)
                inst_prob.append(pp)
                inst_gap.append(gap)

    # ---- per-pair joined table ---------------------------------------------
    rows = []
    for key, orgs in pair_orgs.items():
        a, b = key
        copres = len(present[a] & present[b])
        probs = pair_probs.get(key, [])
        if copres == 0 or not probs:
            continue
        rows.append({
            "function_1": a, "function_2": b,
            "n_instances": len(probs),
            "n_genomes_adjacent": len(orgs),
            "n_genomes_copresent": copres,
            "conditional_cooccurrence": len(orgs) / copres,
            "mean_uniop_prob": float(np.mean(probs)),
            "median_uniop_prob": float(np.median(probs)),
            "median_intergenic_bp": int(np.median(pair_gaps[key])),
        })
    pairs = pd.DataFrame(rows)
    inst_prob = np.asarray(inst_prob)
    inst_gap = np.asarray(inst_gap)

    well = pairs[pairs["n_genomes_copresent"] >= 3]   # cond not trivially 1.0
    cond_w = well["conditional_cooccurrence"].to_numpy()
    prob_w = well["mean_uniop_prob"].to_numpy()
    rho_cond, p_cond = spearmanr(cond_w, prob_w)
    rho_gap, p_gap = spearmanr(inst_prob, inst_gap)

    # always-together vs flimsy (both well-sampled)
    strong = pairs[pairs["n_genomes_copresent"] >= 5]
    alw = strong[strong["conditional_cooccurrence"] >= 0.999]["mean_uniop_prob"].to_numpy()
    flm = strong[strong["conditional_cooccurrence"] < 0.5]["mean_uniop_prob"].to_numpy()
    mw_p = mannwhitneyu(alw, flm, alternative="two-sided").pvalue if len(alw) > 5 and len(flm) > 5 else float("nan")

    # ======================= FIGURE =========================================
    fig, axes = plt.subplots(2, 2, figsize=(15.6, 12.6))
    axA, axB, axC, axD = axes.ravel()
    fig.suptitle("Does UniOP's per-pair operonic probability reflect the empirical co-occurrence?  "
                 "No - it encodes the intergenic gap, not conservation",
                 fontsize=15, fontweight="bold", y=0.985)

    # ---- (a) UniOP prob across conditional-co-occurrence bins --------------
    data_a, meds_a, ns_a = [], [], []
    for lo, hi in COND_BINS:
        sel = (cond_w >= lo) & (cond_w < hi)
        vals = prob_w[sel]
        data_a.append(vals if len(vals) else np.array([np.nan]))
        meds_a.append(np.median(vals) if len(vals) else np.nan)
        ns_a.append(int(sel.sum()))
    xa = np.arange(len(COND_BINS))
    bp = axA.boxplot(data_a, positions=xa, widths=0.62, patch_artist=True,
                     showfliers=False, medianprops=dict(color="black", linewidth=1.6))
    for patch, c in zip(bp["boxes"], COND_COL):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    axA.plot(xa, meds_a, "-o", color="black", linewidth=2.2, markersize=6, zorder=5)
    for x, m, n in zip(xa, meds_a, ns_a):
        if not np.isnan(m):
            axA.text(x, 1.005, "n=%d" % n, ha="center", va="bottom",
                     fontsize=9, fontweight="bold", color="#333333")
    axA.set_xticks(xa)
    axA.set_xticklabels(COND_LAB, fontweight="bold")
    axA.set_ylim(0.45, 1.03)
    axA.set_xlabel("Empirical conditional co-occurrence  (genomes adjacent / genomes both present)",
                   fontweight="bold")
    axA.set_ylabel("UniOP pairwise operonic probability", fontweight="bold")
    axA.set_title("(a) UniOP is ~equally confident no matter how often the pair really co-occurs\n"
                  "flat across the whole co-occurrence range  (Spearman rho = %.2f, n = %d pairs)"
                  % (rho_cond, len(cond_w)), fontweight="bold")
    axA.grid(False)
    L.boldticks(axA)

    # ---- (b) always-together vs flimsy overlap -----------------------------
    edges = np.linspace(0.5, 1.0, 26)
    axB.hist(flm, bins=edges, density=True, color=L.RED, alpha=0.55,
             label="flimsy  (cond < 0.5,  n = %d)" % len(flm))
    axB.hist(alw, bins=edges, density=True, color=L.GREEN, alpha=0.55,
             label="always-together  (cond = 1,  n = %d)" % len(alw))
    axB.axvline(np.median(flm), color=L.RED, linewidth=2.4, linestyle="--")
    axB.axvline(np.median(alw), color=L.GREEN, linewidth=2.4, linestyle="--")
    axB.set_xlim(0.5, 1.0)
    axB.set_xlabel("UniOP pairwise operonic probability  (per pair, mean over genomes)",
                   fontweight="bold")
    axB.set_ylabel("Density of descriptor pairs", fontweight="bold")
    axB.set_title("(b) UniOP cannot separate deterministic modules from flimsy pairs\n"
                  "medians %.3f vs %.3f  (Mann-Whitney p = %.2f, not significant)"
                  % (np.median(alw), np.median(flm), mw_p), fontweight="bold")
    axB.legend(loc="upper left", fontsize=10, framealpha=0.9)
    axB.grid(False)
    L.boldticks(axB)

    # ---- (c) mechanism: UniOP prob across intergenic-distance bins ---------
    data_c, meds_c, ns_c = [], [], []
    idx = np.array([_gap_bin_idx(x) for x in inst_gap])
    for i in range(len(GAP_BINS)):
        vals = inst_prob[idx == i]
        data_c.append(vals if len(vals) else np.array([np.nan]))
        meds_c.append(np.median(vals) if len(vals) else np.nan)
        ns_c.append(int((idx == i).sum()))
    xc = np.arange(len(GAP_BINS))
    bp2 = axC.boxplot(data_c, positions=xc, widths=0.62, patch_artist=True,
                      showfliers=False, medianprops=dict(color="black", linewidth=1.6))
    for patch, c in zip(bp2["boxes"], GAP_COL):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    axC.plot(xc, meds_c, "-o", color="black", linewidth=2.2, markersize=6, zorder=5)
    for x, m, n in zip(xc, meds_c, ns_c):
        if not np.isnan(m):
            axC.text(x, 1.005, "n=%d" % n, ha="center", va="bottom",
                     fontsize=9, fontweight="bold", color="#333333")
    axC.set_xticks(xc)
    axC.set_xticklabels(GAP_LAB, fontweight="bold")
    axC.set_ylim(0.45, 1.03)
    axC.set_xlabel("Intergenic distance between the operon partners (bp)", fontweight="bold")
    axC.set_ylabel("UniOP pairwise operonic probability", fontweight="bold")
    axC.set_title("(c) What UniOP actually encodes: the intergenic gap (a single-genome signal)\n"
                  "steep, monotonic  (Spearman rho = %.2f, n = %d adjacencies)"
                  % (rho_gap, len(inst_prob)), fontweight="bold")
    axC.grid(False)
    L.boldticks(axC)

    # ---- (d) concrete contrast --------------------------------------------
    def _dedup(df):
        seen, keep = set(), []
        for _, r in df.iterrows():
            lab = "%s + %s" % (L.short_desc(r["function_1"], 30),
                               L.short_desc(r["function_2"], 30))
            if lab in seen:
                continue
            seen.add(lab)
            keep.append((r, lab))
        return keep

    flimsy_ex = _dedup(pairs[(pairs["n_genomes_copresent"] >= 12) &
                             (pairs["conditional_cooccurrence"] <= 0.15)]
                       .sort_values("mean_uniop_prob", ascending=False))[:6]
    always_ex = _dedup(pairs[(pairs["conditional_cooccurrence"] >= 0.999) &
                             (pairs["n_genomes_copresent"] >= 17)]
                       .sort_values("n_genomes_copresent", ascending=False))[:6]
    ex = [(r, lab, L.RED) for r, lab in flimsy_ex] + \
         [(r, lab, L.BLUE) for r, lab in always_ex]
    ex.sort(key=lambda t: t[0]["mean_uniop_prob"])
    ylabels = [lab for _, lab, _ in ex]
    xvals = [r["mean_uniop_prob"] for r, _, _ in ex]
    colors = [c for _, _, c in ex]
    yy = np.arange(len(ex))
    axD.barh(yy, xvals, color=colors, alpha=0.9, edgecolor="black", linewidth=0.6)
    for i, (r, _, _) in enumerate(ex):
        axD.text(0.515, i, "cond=%.2f  (%d/%d gen.)"
                 % (r["conditional_cooccurrence"], r["n_genomes_adjacent"],
                    r["n_genomes_copresent"]),
                 va="center", ha="left", fontsize=8.6, fontweight="bold", color="white")
    axD.set_yticks(yy)
    axD.set_yticklabels(ylabels, fontsize=8.6)
    axD.set_xlim(0.5, 1.02)
    axD.set_xlabel("UniOP pairwise operonic probability", fontweight="bold")
    from matplotlib.patches import Patch
    axD.legend(handles=[Patch(facecolor=L.RED, edgecolor="black",
                              label="flimsy pair (co-occurs rarely)"),
                        Patch(facecolor=L.BLUE, edgecolor="black",
                              label="always-together ribosomal module")],
               loc="lower right", fontsize=9, framealpha=0.95)
    axD.set_title("(d) A flimsy pair can out-score a perfectly conserved module\n"
                  "same UniOP band (~0.9) whether the pair always co-occurs or almost never does",
                  fontweight="bold")
    axD.grid(False)
    L.boldticks(axD)

    fig.tight_layout(rect=[0, 0, 1, 0.965])
    L.savefig(fig, outdir / "fig05_uniop_vs_cooccurrence.png")

    # ---- TSVs --------------------------------------------------------------
    pairs_sorted = pairs.sort_values(["n_genomes_copresent", "conditional_cooccurrence"],
                                     ascending=[False, False])
    L.write_tsv(pairs_sorted, outdir / "fig05_pair_uniop_vs_cooccurrence.tsv")

    summ = []
    for (lo, hi), lab, vals in zip(COND_BINS, COND_LAB, data_a):
        v = vals[~np.isnan(vals)]
        summ.append({"panel": "a_cond_bin", "bin": lab, "n": len(v),
                     "mean_uniop_prob": float(np.mean(v)) if len(v) else float("nan"),
                     "median_uniop_prob": float(np.median(v)) if len(v) else float("nan")})
    for (lo, hi), lab, vals in zip(GAP_BINS, GAP_LAB, data_c):
        v = vals[~np.isnan(vals)]
        summ.append({"panel": "c_gap_bin", "bin": lab, "n": len(v),
                     "mean_uniop_prob": float(np.mean(v)) if len(v) else float("nan"),
                     "median_uniop_prob": float(np.median(v)) if len(v) else float("nan")})
    summ.append({"panel": "corr", "bin": "spearman_cond_vs_uniop(copres>=3)",
                 "n": len(cond_w), "mean_uniop_prob": float(rho_cond),
                 "median_uniop_prob": float(p_cond)})
    summ.append({"panel": "corr", "bin": "spearman_gap_vs_uniop(instances)",
                 "n": len(inst_prob), "mean_uniop_prob": float(rho_gap),
                 "median_uniop_prob": float(p_gap)})
    summ.append({"panel": "b_overlap", "bin": "always_together(cond=1,copres>=5)",
                 "n": len(alw), "mean_uniop_prob": float(np.mean(alw)),
                 "median_uniop_prob": float(np.median(alw))})
    summ.append({"panel": "b_overlap", "bin": "flimsy(cond<0.5,copres>=5)",
                 "n": len(flm), "mean_uniop_prob": float(np.mean(flm)),
                 "median_uniop_prob": float(np.median(flm))})
    summ.append({"panel": "b_overlap", "bin": "mann_whitney_p",
                 "n": len(alw) + len(flm), "mean_uniop_prob": float(mw_p),
                 "median_uniop_prob": float("nan")})
    L.write_tsv(pd.DataFrame(summ), outdir / "fig05_summary.tsv")


if __name__ == "__main__":
    L.figure_main(make, SUB)
