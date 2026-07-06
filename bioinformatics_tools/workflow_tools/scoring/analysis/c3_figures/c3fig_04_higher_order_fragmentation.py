#!/usr/bin/env python3
"""Figure 58 - Higher-order fragmentation: does a growing module stay whole?

THE QUESTION (scaling the operon-boundary / partner-switching analysis up)
--------------------------------------------------------------------------
fig 02 asked, for a PAIR, whether a member is separated into another operon or
left a singleton.  Here we ask the same for k-gene modules: when all k members
of a conserved module are present in a genome, how often is the WHOLE set really
one operon - and when it is not, does a member separate into a DIFFERENT operon
or drop out to a SINGLETON?  This is the fragmentation side of the trust question:
how much can we trust that a counted k-gene partnership is a real, intact operon?

DEFINITIONS  (gene identity = clean_descriptor; consistent with figs 01-03)
---------------------------------------------------------------------------
For a k-module S (k consecutive informative operon members in >=1 genome) and a
genome where ALL k descriptors are present, classify the instance:
    intact       the whole set is one consecutive operon run
    rearranged   all k share a single operon but are not a consecutive run
    split        all k are operonic but occupy >=2 DIFFERENT operons
    >=1 singleton  at least one member has NO operonic copy (dropped out)
(priority when ambiguous: a singleton dropout is the most severe, then split,
then rearranged.)  "separated" = split OR singleton (a member truly leaves the
module's operon); rearranged stays inside one operon.

Panels
------
(a) Fate of every co-present instance, per module size k: intact vs rearranged
    vs split vs >=1-singleton (stacked).  Full assembly is the exception.
(b) THE FRAGMENTATION CURVE: fraction of instances where >=1 member separates
    into another operon (split) or drops to a singleton, vs k - with the intact
    fraction for reference.
(c) Recurrent modules are far more trustable: fate breakdown for modules seen
    together in >=5 genomes vs all modules - the count IS reliable for real
    modules, not for the incidental long tail.
(d) A worked example: one conserved multi-gene module drawn across all 21
    genomes, each cell coloured by its fate - intact in some, split or
    singleton-fragmented in others (partner switching at higher order).

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
from matplotlib.patches import Rectangle

csv.field_size_limit(10_000_000)

SUB = "01-operon-context-confidence"
KS = [2, 3, 4, 5, 6, 7]

FATE_ORDER = ["intact", "rearranged", "split", "singleton"]
FATE_COLOR = {"intact": L.BLUE, "rearranged": L.CYAN, "split": L.ORANGE,
              "singleton": L.RED, "absent": "#e4e4e4"}
FATE_LABEL = {"intact": "intact\n(one operon run)",
              "rearranged": "rearranged\n(same operon,\nnon-adjacent)",
              "split": "split\n(>=2 operons)",
              "singleton": ">=1 singleton\n(member dropped out)"}


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


def compute(genes, run_root):
    g = genes.copy()
    organisms = sorted(g["organism"].unique())
    contig = _build_contig_map(run_root, organisms)
    g["contig"] = [contig.get((o, f)) for o, f in zip(g["organism"], g["feature_id"])]

    present = defaultdict(set)
    desc_ops = defaultdict(lambda: defaultdict(set))    # org -> d -> {operon_id}
    for o, cd, un, oid, io, mc in zip(g["organism"], g["clean_descriptor"],
                                      g["uninformative"], g["operon_id"],
                                      g["in_operon"], g["member_count"]):
        if bool(un):
            continue
        d = _norm(cd)
        present[d].add(o)
        if bool(io) and int(mc) >= 2 and str(oid).startswith("operon_"):
            desc_ops[o][d].add(oid)

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
    return organisms, present, desc_ops, mods


def _fate(key, org, present, desc_ops, together):
    members = set(key)
    if not all(org in present[d] for d in members):
        return "absent"
    if org in together:
        return "intact"
    if any(len(desc_ops[org].get(d, set())) == 0 for d in members):
        return "singleton"
    common = set.intersection(*[desc_ops[org].get(d, set()) for d in members])
    return "rearranged" if common else "split"


def make(genes, operons, outdir):
    run_root = outdir.parents[3]
    organisms, present, desc_ops, mods = compute(genes, run_root)

    # ---- fate tallies per k (all co-present instances) --------------------
    tally = {k: defaultdict(int) for k in KS}
    tally_recur = {k: defaultdict(int) for k in KS}   # modules together in >=5 genomes
    for k in KS:
        for key, (inst, together, gaps) in mods[k].items():
            members = set(key)
            cp = set.intersection(*[present[d] for d in members]) if members else set()
            if len(cp) < 2:
                continue
            recur = len(together) >= 5
            for org in cp:
                fate = _fate(key, org, present, desc_ops, together)
                tally[k][fate] += 1
                if recur:
                    tally_recur[k][fate] += 1

    def fractions(t):
        tot = sum(t[f] for f in FATE_ORDER)
        return {f: (100 * t[f] / tot if tot else 0.0) for f in FATE_ORDER}, tot

    fr = {k: fractions(tally[k]) for k in KS}

    fig, axes = plt.subplots(2, 2, figsize=(15.8, 12.8))
    axA, axB, axC, axD = axes.ravel()

    # ---- (a) stacked fates per k ------------------------------------------
    x = np.arange(len(KS))
    bottom = np.zeros(len(KS))
    for fate in FATE_ORDER:
        vals = np.array([fr[k][0][fate] for k in KS])
        axA.bar(x, vals, 0.66, bottom=bottom, color=FATE_COLOR[fate],
                edgecolor="black", linewidth=0.5, label=FATE_LABEL[fate])
        for xi, v, b in zip(x, vals, bottom):
            if v >= 4:
                axA.annotate("%.0f%%" % v, (xi, b + v / 2), ha="center",
                             va="center", fontsize=8, fontweight="bold",
                             color="white" if fate != "rearranged" else "black")
        bottom += vals
    axA.set_xticks(x)
    axA.set_xticklabels(["%d" % k for k in KS])
    axA.set_ylim(0, 100)
    axA.set_xlabel("Module size k (genes)", fontweight="bold")
    axA.set_ylabel("Co-present instances  %", fontweight="bold")
    axA.set_title("(a) Fate of a module across its co-present genomes\n"
                  "full assembly is the exception at every size")
    axA.legend(frameon=True, loc="upper center", bbox_to_anchor=(0.5, -0.12),
               ncol=4, fontsize=7.6)
    L.boldticks(axA)
    axA.grid(False)

    # ---- (b) fragmentation curve ------------------------------------------
    kx = np.array(KS)
    intact = np.array([fr[k][0]["intact"] for k in KS])
    split = np.array([fr[k][0]["split"] for k in KS])
    singl = np.array([fr[k][0]["singleton"] for k in KS])
    separated = split + singl
    axB.plot(kx, separated, color=L.RED, lw=2.4, marker="o", markersize=7,
             markerfacecolor="white", markeredgewidth=1.6,
             label=">=1 member separated (split or singleton)")
    axB.plot(kx, split, color=L.ORANGE, lw=1.9, marker="^", markersize=6,
             label="  ...into another operon (split)")
    axB.plot(kx, singl, color=L.PURPLE, lw=1.9, marker="v", markersize=6,
             label="  ...dropped to a singleton")
    axB.plot(kx, intact, color=L.BLUE, lw=1.9, marker="s", markersize=6, ls="--",
             label="intact (whole set = one operon)")
    for xi, yv in zip(kx, separated):
        axB.annotate("%.0f%%" % yv, (xi, yv), xytext=(0, 8),
                     textcoords="offset points", ha="center", fontsize=8.5,
                     fontweight="bold", color=L.RED)
    axB.set_ylim(0, 100)
    axB.set_xticks(kx)
    axB.set_xlabel("Module size k (genes)", fontweight="bold")
    axB.set_ylabel("Co-present instances  %", fontweight="bold")
    axB.set_title("(b) Fragmentation curve: how often at least one member\n"
                  "leaves the module's operon, vs module size")
    axB.legend(frameon=True, loc="center right", fontsize=8.2)
    L.boldticks(axB)
    axB.grid(False)

    # ---- (c) recurrent modules vs all -------------------------------------
    fr_rec = {k: fractions(tally_recur[k]) for k in KS}
    all_intact = np.array([fr[k][0]["intact"] for k in KS])
    rec_intact = np.array([fr_rec[k][0]["intact"] for k in KS])
    wg = 0.38
    axC.bar(x - wg / 2, all_intact, wg, color="#b0b0b0", edgecolor="black",
            linewidth=0.5, label="all co-present modules")
    axC.bar(x + wg / 2, rec_intact, wg, color=L.GREEN, edgecolor="black",
            linewidth=0.5, label="recurrent modules (together in >=5 genomes)")
    for xi, v in zip(x - wg / 2, all_intact):
        axC.annotate("%.0f%%" % v, (xi, v), ha="center", va="bottom", fontsize=7.6,
                     fontweight="bold")
    for xi, v in zip(x + wg / 2, rec_intact):
        axC.annotate("%.0f%%" % v, (xi, v), ha="center", va="bottom", fontsize=7.6,
                     fontweight="bold", color=L.GREEN)
    axC.set_xticks(x)
    axC.set_xticklabels(["%d" % k for k in KS])
    axC.set_ylim(0, 100)
    axC.set_xlabel("Module size k (genes)", fontweight="bold")
    axC.set_ylabel("Instances fully intact  %", fontweight="bold")
    axC.set_title("(c) The partnership count is trustable for REAL modules\n"
                  "recurrent modules stay intact far more often than the tail")
    axC.legend(frameon=True, loc="upper right", fontsize=8.2)
    L.boldticks(axC)
    axC.grid(False)

    # ---- (d) worked example across all genomes ----------------------------
    _draw_example(axD, organisms, present, desc_ops, mods)

    fig.suptitle("Higher-order fragmentation: a counted k-gene partnership is a "
                 "reliable operon only for small, recurrent modules",
                 fontsize=13, fontweight="bold", y=1.005)
    fig.tight_layout(h_pad=3.1, w_pad=3.0)

    # ---- TSVs --------------------------------------------------------------
    rows = []
    for k in KS:
        fs, tot = fr[k]
        fsr, totr = fr_rec[k]
        rows.append({
            "k": k, "n_instances": tot,
            "pct_intact": round(fs["intact"], 2),
            "pct_rearranged": round(fs["rearranged"], 2),
            "pct_split_other_operon": round(fs["split"], 2),
            "pct_singleton_dropout": round(fs["singleton"], 2),
            "pct_separated": round(fs["split"] + fs["singleton"], 2),
            "n_instances_recurrent": totr,
            "pct_intact_recurrent": round(fsr["intact"], 2),
        })
    L.write_tsv(pd.DataFrame(rows), outdir / "fig04_fragmentation_by_size.tsv")

    L.savefig(fig, outdir / "fig04_higher_order_fragmentation.png")


def _pick_example(present, desc_ops, mods):
    """A conserved size-4/5 module whose genomes show the most fate variety."""
    best = None
    for k in (5, 4):
        for key, (inst, together, gaps) in mods[k].items():
            members = set(key)
            cp = set.intersection(*[present[d] for d in members]) if members else set()
            if len(cp) < 12 or len(together) < 6:
                continue
            fates = defaultdict(int)
            for org in cp:
                fates[_fate(key, org, present, desc_ops, together)] += 1
            variety = sum(1 for f in ("intact", "split", "singleton") if fates[f] > 0)
            score = (variety, len(cp), min(fates["split"] + fates["singleton"],
                                           fates["intact"]))
            if best is None or score > best[0]:
                best = (score, key, together)
    if best is None:
        return None, None
    return best[1], best[2]


def _draw_example(ax, organisms, present, desc_ops, mods):
    key, together = _pick_example(present, desc_ops, mods)
    if key is None:
        ax.axis("off")
        return
    orgs = sorted(organisms)
    ncol = 7
    nrow = int(np.ceil(len(orgs) / ncol))
    for idx, org in enumerate(orgs):
        rr = nrow - 1 - (idx // ncol)
        cc = idx % ncol
        fate = _fate(key, org, present, desc_ops, together)
        ax.add_patch(Rectangle((cc, rr), 0.92, 0.92, facecolor=FATE_COLOR[fate],
                               edgecolor="black", linewidth=0.7))
        ax.text(cc + 0.46, rr + 0.46, L.short_label(org), ha="center", va="center",
                fontsize=6.6, fontweight="bold", fontstyle="italic",
                color="white" if fate in ("intact", "split", "singleton") else "black")
    ax.set_xlim(-0.1, ncol + 0.05)
    ax.set_ylim(-0.15, nrow + 0.05)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(False)
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=FATE_COLOR[f], edgecolor="black", label=lab)
               for f, lab in [("intact", "intact"), ("rearranged", "rearranged"),
                              ("split", "split (other operon)"),
                              ("singleton", "singleton"), ("absent", "not co-present")]]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.16),
              ncol=5, frameon=True, fontsize=7.6)
    mod = " + ".join(L.short_desc(d, 16) for d in key)
    if len(mod) > 74:
        mod = mod[:73] + "\u2026"
    n_int = sum(1 for o in orgs if _fate(key, o, present, desc_ops, together) == "intact")
    n_cp = sum(1 for o in orgs if all(o in present[d] for d in set(key)))
    ax.set_title("(d) One %d-gene module across 21 genomes  (intact in %d/%d present)\n%s"
                 % (len(key), n_int, n_cp, mod), fontsize=9.2, pad=10)
    ax.grid(False)


if __name__ == "__main__":
    L.figure_main(make, SUB)
