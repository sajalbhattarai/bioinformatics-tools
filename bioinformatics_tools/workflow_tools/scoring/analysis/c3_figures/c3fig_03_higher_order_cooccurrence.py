#!/usr/bin/env python3
"""Figure 57 - Higher-order operon modules: co-occurrence & the trust-decay curve.

THE QUESTION (scaling the pairwise operon co-occurrence up)
-----------------------------------------------------------
fig 01 measured how reliably TWO genes co-occur in an operon.  The natural
question is: what about THREE genes together, four, five, ...?  Which SETS of
genes always appear in one operon, how are they spaced, and - crucially - how
does the trustworthiness of the "they always operon-together" claim decay as the
set grows?  That decay curve tells us how far the partnership count can be pushed
before it stops being reliable.

DEFINITIONS  (identical spirit to fig 01; gene identity = clean_descriptor)
---------------------------------------------------------------------------
* k-module = a set of k informative descriptors that occur as k CONSECUTIVE
  operon members (sorted by start, same contig, same operon) in >=1 genome.
  Its canonical key is the sorted descriptor tuple (orientation-independent).
* genomes_together(S) = # genomes with such a consecutive run of exactly S.
* co-present(S) = # genomes where ALL k descriptors are present (operonic OR
  not - a member that dropped to a singleton is still "present").
* conditional co-occurrence  =  genomes_together / co-present   in [0,1]
      1.0  ==  wherever all k functions exist, they are always one operon run.
* weakest-link intergenic gap = MAX of the k-1 internal gaps of a run (the gap
  most likely to break the module).

Panels
------
(a) The pool of trustable modules THINS with k: distinct modules, multi-genome
    (>=2) modules, and "always-together" modules (cond = 1 across >=3 genomes).
(b) The weakest internal gap WIDENS with k: per-module median max-internal-gap.
(c) THE TRUST-DECAY CURVE: flagship reliability - the conditional co-occurrence
    of the most-recurrent size-k modules falls as k grows (even the best-
    conserved module assembles fully in fewer genomes).
(d) The flagship module at each size k (the gene set, its genome count, its
    conditional co-occurrence and spacing) - the ribosomal super-operon extended.

Read-only analysis; does NOT modify the scoring pipeline.
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

csv.field_size_limit(10_000_000)

SUB = "01-operon-context-confidence"
KS = [2, 3, 4, 5, 6, 7]
# small set = cool / trustable  ->  large set = warm / fragile
KCOLOR = {2: L.BLUE, 3: L.CYAN, 4: L.GREEN, 5: L.AMBER, 6: L.ORANGE, 7: L.RED}


def _norm(d):
    return (d or "").strip().lower()


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
            iF, iG = header.index("feature_id"), header.index("gene_id")
            for row in rd:
                if len(row) > max(iF, iG):
                    hit = rx.match(row[iG] or "")
                    m[(org, row[iF])] = hit.group(1) if hit else None
    return m


def compute_modules(genes, run_root):
    """Enumerate consecutive informative k-runs within operons.

    Returns:
        present : descriptor -> set(genome) present (informative genes, any locus)
        mods    : {k: {key: [n_instances, set(genomes_together), [max_gaps]]}}
    """
    g = genes.copy()
    organisms = sorted(g["organism"].unique())
    contig = _build_contig_map(run_root, organisms)
    g["contig"] = [contig.get((o, f)) for o, f in zip(g["organism"], g["feature_id"])]

    present = defaultdict(set)
    for o, cd, un in zip(g["organism"], g["clean_descriptor"], g["uninformative"]):
        if not bool(un):
            present[_norm(cd)].add(o)

    mods = {k: defaultdict(lambda: [0, set(), []]) for k in KS}
    for (org, oid), sub in g.groupby(["organism", "operon_id"]):
        if not str(oid).startswith("operon_"):
            continue
        r = sub.sort_values("start").to_dict("records")
        n = len(r)
        for k in KS:
            for i in range(n - k + 1):
                win = r[i:i + k]
                if any(bool(w["uninformative"]) for w in win):
                    continue
                if len(set(w["contig"] for w in win)) > 1:
                    continue
                ds = [_norm(w["clean_descriptor"]) for w in win]
                if any(not d for d in ds):
                    continue
                key = tuple(sorted(ds))
                slot = mods[k][key]
                slot[0] += 1
                slot[1].add(org)
                slot[2].append(max(int(win[j + 1]["start"]) - int(win[j]["end"]) - 1
                                   for j in range(k - 1)))
    return present, mods


def _module_frame(present, mods_k):
    rows = []
    for key, (inst, orgs, gaps) in mods_k.items():
        members = set(key)
        cp = len(set.intersection(*[present[d] for d in members])) if members else 0
        gt = len(orgs)
        rows.append({
            "key": key, "n_instances": inst, "genomes_together": gt,
            "co_present": cp, "conditional": (gt / cp) if cp else float("nan"),
            "median_max_gap": int(np.median(gaps)) if gaps else 0,
        })
    return pd.DataFrame(rows)


def make(genes, operons, outdir):
    run_root = outdir.parents[3]
    present, mods = compute_modules(genes, run_root)
    frames = {k: _module_frame(present, mods[k]) for k in KS}

    # ---- per-k aggregates --------------------------------------------------
    agg = []
    for k in KS:
        f = frames[k]
        multi = f[f["co_present"] >= 2]
        always = f[(f["co_present"] >= 3) & (f["conditional"] >= 0.999)]
        # flagship reliability: conditional of the most-recurrent modules
        recur = multi.sort_values("genomes_together", ascending=False)
        top10 = recur.head(10)["conditional"].to_numpy()
        agg.append({
            "k": k, "n_modules": len(f),
            "n_multi_genome": len(multi), "n_always_together": len(always),
            "n_robust_ge5": int((f["genomes_together"] >= 5).sum()),
            "n_robust_ge10": int((f["genomes_together"] >= 10).sum()),
            "flagship_cond": float(recur.head(1)["conditional"].iloc[0]) if len(recur) else float("nan"),
            "flagship_top10_mean": float(np.mean(top10)) if len(top10) else float("nan"),
            "flagship_top10_min": float(np.min(top10)) if len(top10) else float("nan"),
            "flagship_top10_max": float(np.max(top10)) if len(top10) else float("nan"),
            "median_max_gap": float(f["median_max_gap"].median()) if len(f) else float("nan"),
        })
    A = pd.DataFrame(agg)

    fig, axes = plt.subplots(2, 2, figsize=(15.8, 12.6))
    axA, axB, axC, axD = axes.ravel()

    # ---- (a) thinning of the trustable pool -------------------------------
    x = np.arange(len(KS))
    w = 0.26
    axA.bar(x - w, A["n_modules"], w, color="#9e9e9e", edgecolor="black",
            linewidth=0.5, label="distinct modules")
    axA.bar(x, A["n_multi_genome"], w, color=L.BLUE, edgecolor="black",
            linewidth=0.5, label="in >=2 genomes")
    axA.bar(x + w, A["n_always_together"], w, color=L.RED, edgecolor="black",
            linewidth=0.5, label="always together (cond=1, >=3 genomes)")
    axA.set_yscale("log")
    axA.set_xticks(x)
    axA.set_xticklabels(["%d" % k for k in KS])
    for xi, n in zip(x, A["n_always_together"]):
        axA.annotate("%d" % n, (xi + w, n), ha="center", va="bottom",
                     fontsize=8, fontweight="bold", color=L.RED)
    axA.set_xlabel("Module size k (genes)", fontweight="bold")
    axA.set_ylabel("Number of modules (log)", fontweight="bold")
    axA.set_title("(a) The pool of trustable modules thins with size\n"
                  "always-together sets: %d pairs -> %d of size 7"
                  % (A["n_always_together"].iloc[0], A["n_always_together"].iloc[-1]))
    axA.legend(frameon=True, loc="upper right", fontsize=8.5)
    L.boldticks(axA)
    axA.grid(False)

    # ---- (b) weakest-link gap widens --------------------------------------
    box_data = [frames[k]["median_max_gap"].to_numpy() for k in KS]
    bp = axB.boxplot(box_data, patch_artist=True, showfliers=False, widths=0.62)
    for patch, k in zip(bp["boxes"], KS):
        patch.set_facecolor(KCOLOR[k])
        patch.set_alpha(0.78)
        patch.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black")
        med.set_linewidth(1.6)
    axB.plot(np.arange(1, len(KS) + 1), A["median_max_gap"], color="black",
             lw=1.6, marker="o", markersize=5, markerfacecolor="white", zorder=5)
    axB.set_xticklabels(["%d" % k for k in KS])
    axB.set_xlabel("Module size k (genes)", fontweight="bold")
    axB.set_ylabel("Weakest-link internal gap  (max of k-1 gaps, bp)",
                   fontweight="bold")
    axB.set_title("(b) The weakest internal gap widens with size\n"
                  "median %d bp (pairs) -> %d bp (size 7)"
                  % (A["median_max_gap"].iloc[0], A["median_max_gap"].iloc[-1]))
    L.boldticks(axB)
    axB.grid(False)

    # ---- (c) THE TRUST-DECAY CURVE ----------------------------------------
    kx = np.array(KS)
    axC.fill_between(kx, A["flagship_top10_min"], A["flagship_top10_max"],
                     color=L.BLUE, alpha=0.16,
                     label="top-10 recurrent modules (spread)")
    axC.plot(kx, A["flagship_top10_mean"], color=L.BLUE, lw=2.4, marker="o",
             markersize=7, markerfacecolor="white", markeredgewidth=1.6,
             label="top-10 recurrent modules (mean)")
    axC.plot(kx, A["flagship_cond"], color=L.RED, lw=2.0, marker="s",
             markersize=6, ls="--", label="single most-recurrent module")
    for xi, yv in zip(kx, A["flagship_top10_mean"]):
        axC.annotate("%.2f" % yv, (xi, yv), xytext=(0, 8),
                     textcoords="offset points", ha="center", fontsize=8.5,
                     fontweight="bold", color=L.BLUE)
    axC.set_ylim(0, 1.05)
    axC.set_xticks(kx)
    axC.set_xlabel("Module size k (genes)", fontweight="bold")
    axC.set_ylabel("Conditional co-occurrence  (fully assembled where present)",
                   fontweight="bold")
    axC.set_title("(c) Trust-decay curve: even the best-conserved module\n"
                  "assembles fully in fewer genomes as it grows")
    axC.legend(frameon=True, loc="lower left", fontsize=8.8)
    L.boldticks(axC)
    axC.grid(False)

    # ---- (d) the flagship module at each size -----------------------------
    flag = []
    for k in KS:
        multi = frames[k][frames[k]["co_present"] >= 2]
        if not len(multi):
            continue
        best = multi.sort_values(["genomes_together", "conditional"],
                                 ascending=False).iloc[0]
        flag.append((k, best))
    y = np.arange(len(flag))
    gts = [b["genomes_together"] for _, b in flag]
    cols = [KCOLOR[k] for k, _ in flag]
    axD.barh(y, gts, color=cols, edgecolor="black", linewidth=0.6, height=0.66)
    for yi, (k, b) in zip(y, flag):
        members = [L.short_desc(d, 15) for d in b["key"]]
        label = " + ".join(members)
        if len(label) > 62:
            label = label[:61] + "\u2026"
        axD.annotate("k=%d  %d gen.  cond=%.2f  gap=%dbp\n%s"
                     % (k, b["genomes_together"], b["conditional"],
                        b["median_max_gap"], label),
                     (b["genomes_together"], yi), xytext=(5, 0),
                     textcoords="offset points", va="center", ha="left",
                     fontsize=7.0, fontweight="bold")
    axD.set_yticks(y)
    axD.set_yticklabels(["k=%d" % k for k, _ in flag])
    axD.set_xlim(0, max(gts) * 1.75 if gts else 1)
    axD.set_xlabel("Genomes where the whole set is one operon run",
                   fontweight="bold")
    axD.set_title("(d) The flagship module at each size\n"
                  "(the ribosomal super-operon, extended gene by gene)")
    L.boldticks(axD)
    axD.grid(False)

    fig.suptitle("Higher-order operon modules: which gene sets always co-occur, "
                 "and how the partnership's trust decays with set size",
                 fontsize=13, fontweight="bold", y=1.004)
    fig.tight_layout(h_pad=2.6, w_pad=3.0)

    # ---- TSVs --------------------------------------------------------------
    L.write_tsv(A, outdir / "fig03_ksize_summary.tsv")
    allrows = []
    for k in KS:
        f = frames[k].copy()
        f.insert(0, "k", k)
        f["module"] = [" | ".join(key) for key in f["key"]]
        f = f.drop(columns=["key"])
        allrows.append(f[f["co_present"] >= 2])
    modules = pd.concat(allrows, ignore_index=True).sort_values(
        ["k", "genomes_together", "conditional"], ascending=[True, False, False])
    L.write_tsv(modules, outdir / "fig03_modules_by_size.tsv")

    always_rows = modules[(modules["co_present"] >= 3) & (modules["conditional"] >= 0.999)]
    L.write_tsv(always_rows, outdir / "fig03_always_together_modules.tsv")

    L.savefig(fig, outdir / "fig03_higher_order_cooccurrence.png")


if __name__ == "__main__":
    L.figure_main(make, SUB)
