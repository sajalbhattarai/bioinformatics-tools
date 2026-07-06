#!/usr/bin/env python3
"""Figure 55 - Intergenic-distance bins of operon co-occurrence + by-chance test.

THE QUESTION (from the operon-index thread)
--------------------------------------------
For two genes that sit together INSIDE the same operon, how far apart are they,
and given that spacing, do they ALWAYS co-occur - or could that partnership be a
coincidence?  A model predicting two adjacent genes are "in an operon" is only
useful if their two FUNCTIONS genuinely belong together, recurrently, across
genomes.  This figure measures exactly that, binned by intergenic distance.

DEFINITIONS
-----------
* Gene identity  = clean_descriptor (the functional name).  The SAME descriptor
  in different genomes is the SAME gene - co-occurrence is counted over
  descriptor pairs, never feature_ids.
* Within-operon adjacency = two genes that are consecutive (sorted by start) in
  the SAME operon on the SAME contig.
* Intergenic distance = downstream.start - upstream.end - 1
        < 0  overlapping ORFs      = 0  abutting      > 0  a real gap (bp)
  Bins: <0, 0, 1-100, 100-200, 200-300, 300-400, 400-500, >500.

BY-CHANCE MODEL (informative-informative pairs only; a hypothetical partner has
no function to co-occur)
-----------------------------------------------------------------------------
Treat every within-operon informative adjacency as an edge of a multigraph on
functions.  M = total edges; deg(A) = edges touching function A.  For a pair
(A,B) observed k_inst times, the configuration-model expectation is
        E = deg(A) * deg(B) / (2M)
        lift    = k_inst / E                     (>1 = enriched over chance)
        p_chance = P(X >= k_inst | Poisson(E))   (the "by chance" probability)
"Do they always co-occur?" = conditional co-occurrence
        cond = (# genomes where A,B are operon-adjacent)
               / (# genomes where both A,B are present as operon members)
cond = 1.0 means: wherever both functions exist, they are always operon partners.

Panels
------
(a) how many within-operon adjacencies fall in each intergenic bin (the spacing
    of operon partners) - operon-internal genes are overwhelmingly overlapping
    or <100 bp apart.
(b) "do they always co-occur?" - distribution of conditional co-occurrence for
    informative descriptor pairs; a large mode at 1.0 = deterministic partners.
(c) by-chance probability per bin - distribution of -log10 p_chance; even the
    least-conserved partners beat chance, conserved modules reach ~1e-50.
(d) the strongest conserved modules (most genomes co-adjacent), annotated with
    their spacing and by-chance probability.

Read-only analysis; does NOT modify the scoring pipeline.
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import poisson

csv.field_size_limit(10_000_000)

SUB = "01-operon-context-confidence"

# intergenic-distance bins requested by the user, + a >500 catch-all
BINS = [(-10**9, -1, "<0"), (0, 0, "0"), (1, 100, "1-100"),
        (101, 200, "100-200"), (201, 300, "200-300"), (301, 400, "300-400"),
        (401, 500, "400-500"), (501, 10**9, ">500")]
BIN_NAMES = [b[2] for b in BINS]
# green (tight / overlapping = strong coupling) -> red (far apart = weak)
BIN_COLORS = [L.GREEN, L.TEAL, L.LIME, L.YELLOW, L.AMBER, L.ORANGE, L.RED, "#8a1220"]


def _bin_of(gap):
    for lo, hi, nm in BINS:
        if lo <= gap <= hi:
            return nm
    return ">500"


def _build_contig_map(run_root, organisms):
    """(organism, feature_id) -> contig id, parsed from gene_id in labeled-genes."""
    rx = __import__("re").compile(r"^(.*)_(\d+)([+-])(\d+)$")
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

    # ---- within-operon adjacencies -----------------------------------------
    inst_bins = defaultdict(int)         # bin -> instance count (ALL adjacencies)
    inst_bins_info = defaultdict(int)    # bin -> instance count (info-info only)
    pair_inst = defaultdict(int)         # (dA,dB) -> #adjacency instances
    pair_orgs = defaultdict(set)         # (dA,dB) -> set(org adjacent)
    pair_gaps = defaultdict(list)        # (dA,dB) -> [gaps]
    deg = defaultdict(int)               # function -> #informative adjacencies
    M = 0
    present = defaultdict(set)           # function -> set(org) operon-present
    all_gaps = []

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
            gap = int(b["start"]) - int(a["end"]) - 1
            nm = _bin_of(gap)
            inst_bins[nm] += 1
            all_gaps.append(gap)
            ua = uninf[(org, a["feature_id"])]
            ub = uninf[(org, b["feature_id"])]
            if ua or ub:
                continue
            da, db = cln[(org, a["feature_id"])], cln[(org, b["feature_id"])]
            if not da or not db or da == db:
                continue
            inst_bins_info[nm] += 1
            key = (da, db) if da < db else (db, da)
            pair_inst[key] += 1
            pair_orgs[key].add(org)
            pair_gaps[key].append(gap)
            deg[da] += 1
            deg[db] += 1
            M += 1

    # ---- per descriptor-pair statistics (by-chance + always-co-occur) ------
    rows = []
    for key, inst in pair_inst.items():
        a, b = key
        k = len(pair_orgs[key])
        gaps = pair_gaps[key]
        med_gap = int(np.median(gaps))
        E = deg[a] * deg[b] / (2.0 * M) if M else 0.0
        lift = inst / E if E > 0 else float("inf")
        p_chance = float(poisson.sf(inst - 1, E)) if E > 0 else 0.0
        copres = len(present[a] & present[b])
        cond = k / copres if copres > 0 else float("nan")
        rows.append({
            "function_1": a, "function_2": b, "n_instances": inst,
            "n_genomes_adjacent": k, "n_genomes_copresent": copres,
            "conditional_cooccurrence": cond, "median_intergenic_bp": med_gap,
            "intergenic_bin": _bin_of(med_gap), "expected_by_chance": E,
            "lift": lift, "p_by_chance": p_chance,
            "neg_log10_p": -math.log10(max(p_chance, 1e-300)),
        })
    pairs = pd.DataFrame(rows)

    # ======================= FIGURE =========================================
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 12.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- (a) intergenic-bin histogram of within-operon adjacencies ---------
    tot_all = sum(inst_bins.values()) or 1
    all_counts = [inst_bins.get(nm, 0) for nm in BIN_NAMES]
    info_counts = [inst_bins_info.get(nm, 0) for nm in BIN_NAMES]
    x = np.arange(len(BIN_NAMES))
    bars = axA.bar(x, all_counts, color=BIN_COLORS, edgecolor="black",
                   linewidth=0.6, width=0.74)
    for xi, c in zip(x, all_counts):
        if c:
            axA.annotate("%d\n%.1f%%" % (c, 100 * c / tot_all), (xi, c),
                         ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    axA.set_xticks(x)
    axA.set_xticklabels(BIN_NAMES, rotation=25, ha="right")
    axA.set_ylim(0, max(all_counts) * 1.16 if max(all_counts) else 1)
    tight = 100 * (inst_bins.get("<0", 0) + inst_bins.get("0", 0)
                   + inst_bins.get("1-100", 0)) / tot_all
    axA.set_xlabel("Intergenic distance between operon partners (bp)",
                   fontweight="bold")
    axA.set_ylabel("Within-operon adjacencies (instances)", fontweight="bold")
    axA.set_title("(a) How close are operon partners?  %d adjacencies\n"
                  "%.1f%% overlap or sit <=100 bp apart" % (tot_all, tight))
    L.boldticks(axA)
    axA.grid(False)

    # ---- (b) do they always co-occur?  conditional co-occurrence -----------
    cond = pairs["conditional_cooccurrence"].dropna().to_numpy()
    multi = pairs[pairs["n_genomes_copresent"] >= 2]["conditional_cooccurrence"].dropna().to_numpy()
    axB.hist(cond, bins=np.linspace(0, 1, 21), color=L.BLUE, edgecolor="black",
             linewidth=0.3, alpha=0.55,
             label="all informative pairs  (n=%d)" % len(cond))
    if len(multi):
        axB.hist(multi, bins=np.linspace(0, 1, 21), color=L.PURPLE,
                 edgecolor="black", linewidth=0.3, alpha=0.75,
                 label="present in >=2 genomes  (n=%d)" % len(multi))
    always = 100 * float((cond >= 0.999).mean()) if len(cond) else 0.0
    axB.axvline(1.0, color=L.RED, lw=1.4, ls="--")
    axB.set_xlabel("Conditional co-occurrence  =  genomes adjacent / genomes both present",
                   fontweight="bold")
    axB.set_ylabel("Informative descriptor pairs", fontweight="bold")
    axB.set_title("(b) Do operon partners ALWAYS co-occur?\n"
                  "%.0f%% of pairs are adjacent in EVERY genome where both exist"
                  % always)
    axB.legend(frameon=True, loc="upper center", fontsize=9)
    L.boldticks(axB)
    axB.grid(False)

    # ---- (c) by-chance probability per bin ---------------------------------
    # well-populated bins get their own box; sparse far bins are pooled as >=100
    near = ["<0", "0", "1-100"]
    box_data, box_labels, box_cols = [], [], []
    for nm, col in zip(BIN_NAMES, BIN_COLORS):
        if nm not in near:
            continue
        vals = pairs[pairs["intergenic_bin"] == nm]["neg_log10_p"].to_numpy()
        if len(vals):
            box_data.append(vals)
            box_labels.append("%s\n(n=%d)" % (nm, len(vals)))
            box_cols.append(col)
    far = pairs[~pairs["intergenic_bin"].isin(near)]["neg_log10_p"].to_numpy()
    if len(far):
        box_data.append(far)
        box_labels.append(">=100\n(n=%d)" % len(far))
        box_cols.append(L.RED)
    if box_data:
        bp = axC.boxplot(box_data, patch_artist=True, showfliers=True,
                         widths=0.6, flierprops=dict(marker="o", markersize=2.4,
                         markerfacecolor="black", alpha=0.35, markeredgecolor="none"))
        for patch, col in zip(bp["boxes"], box_cols):
            patch.set_facecolor(col)
            patch.set_alpha(0.75)
            patch.set_edgecolor("black")
        for med in bp["medians"]:
            med.set_color("black")
            med.set_linewidth(1.6)
        axC.set_xticklabels(box_labels, fontsize=9)
    axC.axhline(-math.log10(0.05), color="black", lw=1.2, ls="--",
                label="p = 0.05")
    axC.set_xlabel("Intergenic distance bin", fontweight="bold")
    axC.set_ylabel("By-chance improbability  -log10 p", fontweight="bold")
    axC.set_title("(c) Could the partnership be chance?\n"
                  "overlapping / close pairs are astronomically non-random")
    axC.legend(frameon=True, loc="upper right", fontsize=9)
    L.boldticks(axC)
    axC.grid(False)

    # ---- (d) strongest conserved modules -----------------------------------
    top = pairs.sort_values(["n_genomes_adjacent", "lift"],
                            ascending=False).head(12).iloc[::-1]
    if len(top):
        labels = ["%s + %s" % (L.short_desc(a, 20), L.short_desc(b, 20))
                  for a, b in zip(top["function_1"], top["function_2"])]
        y = np.arange(len(labels))
        cols = [BIN_COLORS[BIN_NAMES.index(bn)] for bn in top["intergenic_bin"]]
        axD.barh(y, top["n_genomes_adjacent"], color=cols, edgecolor="black",
                 linewidth=0.6, height=0.72)
        for yi, (kk, pv, gp) in enumerate(zip(top["n_genomes_adjacent"],
                                              top["p_by_chance"],
                                              top["median_intergenic_bp"])):
            axD.annotate("%d gen.  %dbp  p=%.0e" % (kk, gp, pv),
                         (kk, yi), xytext=(4, 0), textcoords="offset points",
                         va="center", ha="left", fontsize=7.5, fontweight="bold")
        axD.set_yticks(y)
        axD.set_yticklabels(labels, fontsize=8)
        axD.set_xlim(0, top["n_genomes_adjacent"].max() * 1.42)
        axD.set_xlabel("Genomes where the two functions are operon-adjacent",
                       fontweight="bold")
        axD.set_title("(d) Strongest conserved operon modules\n"
                      "bar colour = intergenic bin (green = overlap/close)")
    L.boldticks(axD)
    axD.grid(False)

    fig.suptitle("Intergenic spacing of operon partners and the probability their "
                 "co-occurrence is chance", fontsize=13, fontweight="bold",
                 y=1.004)
    fig.tight_layout(h_pad=2.6, w_pad=3.0)

    # ---- TSVs --------------------------------------------------------------
    L.write_tsv(pairs.sort_values(["n_genomes_adjacent", "lift"], ascending=False),
                outdir / "fig01_pair_intergenic_cooccurrence.tsv")

    bin_rows = []
    for nm in BIN_NAMES:
        sub = pairs[pairs["intergenic_bin"] == nm]
        cvals = sub["conditional_cooccurrence"].dropna().to_numpy()
        bin_rows.append({
            "intergenic_bin": nm,
            "n_adjacency_instances_all": int(inst_bins.get(nm, 0)),
            "pct_adjacencies_all": round(100 * inst_bins.get(nm, 0) / tot_all, 3),
            "n_adjacency_instances_info": int(inst_bins_info.get(nm, 0)),
            "n_descriptor_pairs": int(len(sub)),
            "median_lift": float(sub["lift"].median()) if len(sub) else float("nan"),
            "median_neg_log10_p": float(sub["neg_log10_p"].median()) if len(sub) else float("nan"),
            "pct_always_cooccur": round(100 * float((cvals >= 0.999).mean()), 2) if len(cvals) else float("nan"),
        })
    L.write_tsv(pd.DataFrame(bin_rows), outdir / "fig01_intergenic_bin_summary.tsv")

    agaps = np.array(all_gaps)
    summ = pd.DataFrame([{
        "n_within_operon_adjacencies": int(tot_all),
        "n_info_info_adjacencies": int(M),
        "n_descriptor_pairs": int(len(pairs)),
        "median_gap_bp": float(np.median(agaps)) if len(agaps) else float("nan"),
        "mean_gap_bp": float(np.mean(agaps)) if len(agaps) else float("nan"),
        "pct_overlap_or_abut": round(100 * float((agaps <= 0).mean()), 2) if len(agaps) else float("nan"),
        "pct_within_100bp": round(100 * float((agaps <= 100).mean()), 2) if len(agaps) else float("nan"),
        "pct_pairs_always_cooccur": round(100 * float((pairs["conditional_cooccurrence"].dropna() >= 0.999).mean()), 2) if len(pairs) else float("nan"),
        "median_lift_all_pairs": float(pairs["lift"].median()) if len(pairs) else float("nan"),
        "median_neg_log10_p_all_pairs": float(pairs["neg_log10_p"].median()) if len(pairs) else float("nan"),
    }])
    L.write_tsv(summ, outdir / "fig01_summary.tsv")

    L.savefig(fig, outdir / "fig01_intergenic_cooccurrence_bins.png")


if __name__ == "__main__":
    L.figure_main(make, SUB)
