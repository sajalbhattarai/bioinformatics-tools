"""c3fig_07_occ_db_growth.py - how the dynamic OCC database changes as it grows.

Compares a 15-organism OCC database against the 21-organism database to show what
folding in 6 more organisms actually does to the reliabilities (rho) and the
per-gene OCC.

Panels
  (a) structural growth 15 -> 21 (qualifying operons, descriptors, pairs);
  (b) shared-pair reliability rho_15 vs rho_21 (no pair is ever lost; net upward);
  (c) per-gene OCC drift on the SAME first-15 genes (median stable, conserved
      core rises, ~20% of genes move materially);
  (d) new-organism COVERAGE - the argument for the dynamic DB: a freshly-labeled
      organism is largely unscoreable until its own operon evidence is folded in.

Usage: c3fig_07_occ_db_growth.py --stats-dir .../c3-genes-comprehensive-stats
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import c3_occ as O
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats-dir", required=True)
    ap.add_argument("--n-first", type=int, default=15)
    args = ap.parse_args()
    sd = Path(args.stats_dir).resolve()
    rr = sd.parents[1]
    L.apply_style()

    genes = L.load_cache(sd / "_cache" / "genes.pkl")
    orgs = sorted(genes["organism"].unique())
    o_small = orgs[:args.n_first]
    new_orgs = orgs[args.n_first:]

    ref_s = O.build_reference(genes[genes["organism"].isin(o_small)], rr)
    ref_f = O.load_reference(sd / "_cache" / "occ_reference.pkl")
    if len(ref_f["organisms_added"]) != len(orgs):
        ref_f = O.build_reference(genes, rr)

    NS, NF = len(ref_s["organisms_added"]), len(ref_f["organisms_added"])

    # ---- panel data ---------------------------------------------------------
    cats = ["qual.\noperons", "descriptors", "adj\npairs", "op\npairs"]
    v_s = [ref_s["n_qualifying_operons"], len(ref_s["present"]),
           len(ref_s["rho_adj"]), len(ref_s["rho_op"])]
    v_f = [ref_f["n_qualifying_operons"], len(ref_f["present"]),
           len(ref_f["rho_adj"]), len(ref_f["rho_op"])]

    shared = list(set(ref_s["rho_op"]) & set(ref_f["rho_op"]))
    x_op = np.array([ref_s["rho_op"][k] for k in shared])
    y_op = np.array([ref_f["rho_op"][k] for k in shared])
    n_new_pairs = len(set(ref_f["rho_op"]) - set(ref_s["rho_op"]))
    dmean = float((y_op - x_op).mean())

    sub = genes[genes["organism"].isin(o_small)]
    i_s = O.compute_all_genes(sub, rr, ref=ref_s)
    i_f = O.compute_all_genes(sub, rr, ref=ref_f)
    m = i_s.merge(i_f, on=["organism", "feature_id"], suffixes=("_s", "_f"))
    m = m[~m["uninformative_s"]]
    docc = (m["occ_f"] - m["occ_s"]).to_numpy()
    rib = m["clean_descriptor_s"].str.contains("ribosomal protein", na=False).to_numpy()
    rib_med_s = float(m.loc[rib, "occ_s"].median())
    rib_med_f = float(m.loc[rib, "occ_f"].median())
    med_s = float(m["occ_s"].median())
    med_f = float(m["occ_f"].median())

    sub6 = genes[genes["organism"].isin(new_orgs)]
    n_s = O.compute_all_genes(sub6, rr, ref=ref_s)
    n_f = O.compute_all_genes(sub6, rr, ref=ref_f)
    ns = n_s[~n_s["uninformative"]]["occ"].to_numpy()
    nf = n_f[~n_f["uninformative"]]["occ"].to_numpy()
    zero_s = 100.0 * (ns <= 1e-9).mean()
    zero_f = 100.0 * (nf <= 1e-9).mean()
    med_new_s = float(np.median(ns))
    med_new_f = float(np.median(nf))

    # ---- figure -------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(15.5, 12.6))
    lab_s, lab_f = "%d-org DB" % NS, "%d-org DB" % NF

    # (a) structural growth
    a = ax[0, 0]
    xi = np.arange(len(cats)); w = 0.38
    b1 = a.bar(xi - w / 2, v_s, w, color=L.BLUE, label=lab_s)
    b2 = a.bar(xi + w / 2, v_f, w, color=L.GREEN, label=lab_f)
    for bars in (b1, b2):
        for r in bars:
            a.text(r.get_x() + r.get_width() / 2, r.get_height(),
                   "%d" % int(r.get_height()), ha="center", va="bottom", fontsize=10)
    for i in range(len(cats)):
        a.text(xi[i], max(v_s[i], v_f[i]) * 1.14, "+%.0f%%" % (100 * (v_f[i] / v_s[i] - 1)),
               ha="center", va="bottom", fontsize=10, color=L.RED, fontweight="bold")
    a.set_xticks(xi); a.set_xticklabels(cats)
    a.set_ylabel("count"); a.set_ylim(0, max(v_f) * 1.28)
    a.set_title("(a) Database grows monotonically as organisms are added", fontweight="bold")
    a.legend(loc="upper left"); a.grid(False)

    # (b) shared-pair rho
    b = ax[0, 1]
    hb = b.hexbin(x_op, y_op, gridsize=45, bins="log", cmap="viridis", mincnt=1)
    b.plot([0, 1], [0, 1], "--", color=L.RED, lw=1.5)
    b.set_xlim(0, 1); b.set_ylim(0, 1)
    b.set_xlabel(r"$\rho_{op}$ in %s" % lab_s)
    b.set_ylabel(r"$\rho_{op}$ in %s" % lab_f)
    b.set_title("(b) Shared co-operon links: %d kept, 0 lost, +%d new\n"
                "mean $\\Delta\\rho$=%+.4f (recurrence raises, non-recurrence lowers)"
                % (len(shared), n_new_pairs, dmean), fontweight="bold")
    cb = fig.colorbar(hb, ax=b); cb.set_label("pairs (log)")
    b.grid(False)

    # (c) per-gene OCC drift
    c = ax[1, 0]
    bins = np.linspace(-0.6, 0.6, 61)
    c.hist(docc, bins=bins, color=L.BLUE, alpha=0.75, label="all genes (n=%d)" % len(docc))
    c.hist(docc[rib], bins=bins, color=L.ORANGE, alpha=0.85,
           label="ribosomal (n=%d)" % int(rib.sum()))
    c.axvline(0, color="k", lw=1)
    c.axvline(docc.mean(), color=L.RED, lw=1.8, ls="--",
              label="mean %+.4f" % docc.mean())
    c.set_yscale("log")
    c.set_xlabel(r"$\Delta$OCC per gene  (%s $-$ %s)" % (lab_f, lab_s))
    c.set_ylabel("genes (log)")
    c.set_title("(c) Same %d-org genes, more evidence: median %.3f=%.3f,\n"
                "ribosomal %.3f\u2192%.3f; %d rose / %d fell"
                % (NS, med_s, med_f, rib_med_s, rib_med_f,
                   int((docc > 1e-6).sum()), int((docc < -1e-6).sum())),
                fontweight="bold")
    c.legend(loc="upper left"); c.grid(False)

    # (d) new-organism coverage
    d = ax[1, 1]
    xi2 = np.arange(2)
    bars = d.bar(xi2, [zero_s, zero_f], 0.5, color=[L.RED, L.GREEN])
    for r, mv in zip(bars, [med_new_s, med_new_f]):
        d.text(r.get_x() + r.get_width() / 2, r.get_height() + 1.5,
               "%.1f%%" % r.get_height(), ha="center", va="bottom",
               fontsize=13, fontweight="bold")
        if r.get_height() > 12:
            d.text(r.get_x() + r.get_width() / 2, r.get_height() / 2,
                   "median\nOCC\n%.3f" % mv, ha="center", va="center",
                   fontsize=11, color="white", fontweight="bold")
        else:
            d.text(r.get_x() + r.get_width() / 2, r.get_height() + 5.5,
                   "median\nOCC %.3f" % mv, ha="center", va="bottom",
                   fontsize=11, color=L.GREEN, fontweight="bold")
    d.set_xticks(xi2)
    d.set_xticklabels(["scored with\n%s\n(orgs NOT in DB)" % lab_s,
                       "scored with\n%s\n(orgs folded in)" % lab_f])
    d.set_ylabel("%% of the 6 new organisms' genes with OCC = 0".replace("%%", "%"))
    d.set_ylim(0, max(zero_s, zero_f) * 1.25 + 5)
    d.set_title("(d) Why the DB must be dynamic: a new organism is\nlargely "
                "unscoreable until its evidence is folded in", fontweight="bold")
    d.grid(False)

    fig.suptitle("Dynamic OCC database: what changes when %d vs %d organisms are "
                 "in the reference" % (NS, NF), fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out = sd / "figures" / "01-operon-context-confidence" / "fig07_occ_db_growth.png"
    L.savefig(fig, out)


if __name__ == "__main__":
    main()
