"""c3_occ.py - Operon Context Confidence (OCC) factor.

Scores each gene in [0, 1] from pan-genome operon co-occurrence of its operon
partners, pooled by clean_descriptor and restricted to operons with an
informative majority. The reference holds additive counts, so organisms fold in
incrementally; finalize_reference() derives the per-pair reliabilities.
"""
import csv
import re
from collections import defaultdict

import numpy as np
from scipy.stats import beta as _beta

csv.field_size_limit(10_000_000)

DEFAULTS = dict(
    delta=0.05,          # Jeffreys lower-bound quantile
    lam=1.0,             # enrichment saturation constant
    use_enrichment=True,  # apply the by-chance safeguard w
    link_floor=0.0,      # ignore links with rho below this (0 = keep all)
)


def _norm(d):
    """Returns the lower-cased, stripped descriptor ("" for None)."""
    return (d or "").strip().lower()


def build_contig_map(run_root, organisms):
    """Returns {(organism, feature_id): contig id}, parsed from gene_id in labeled-genes.tsv."""
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


def _key(a, b):
    """Returns the unordered pair as a sorted tuple."""
    return (a, b) if a < b else (b, a)


REF_VERSION = 2   # bump when the counts schema changes


def new_reference(params=None):
    """Returns an empty OCC reference holding additive pan-genome counts.

    rho_adj / rho_op are filled later by finalize_reference()."""
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    return dict(
        version=REF_VERSION,
        # --- additive counts (the persistent database) --------------------
        present=defaultdict(set),      # descriptor -> {organisms present}
        adj_org=defaultdict(set),      # pair -> {organisms adjacent}
        coop_org=defaultdict(set),     # pair -> {organisms same-operon}
        adj_inst=defaultdict(int),     # pair -> # adjacency instances
        coop_inst=defaultdict(int),    # pair -> # same-operon instances
        deg_adj=defaultdict(int),      # descriptor -> adjacency incidence
        deg_op=defaultdict(int),       # descriptor -> same-operon incidence
        M_adj=0, M_op=0,               # totals for the enrichment null
        n_qualifying_operons=0,
        organisms_added=set(),         # guards against double-counting
        # --- derived (filled by finalize_reference) -----------------------
        rho_adj={}, rho_op={},
        params=p, finalized=False,
    )


def _accumulate_organism(ref, org, org_genes, contig):
    """Adds one organism's qualifying operons to ref's additive counts."""
    uninf = {f: bool(u) for f, u in
             zip(org_genes["feature_id"], org_genes["uninformative"])}
    cln = {f: _norm(c) for f, c in
           zip(org_genes["feature_id"], org_genes["clean_descriptor"])}

    for oid, sub in org_genes.groupby("operon_id"):
        if not str(oid).startswith("operon_"):
            continue
        recs = sub.sort_values("start").to_dict("records")
        # Keeps only operons with a strict informative majority.
        inf_recs = [r for r in recs if not uninf[r["feature_id"]]
                    and cln[r["feature_id"]]]
        n_inf = len(inf_recs)
        n_unf = len(recs) - n_inf
        if n_inf <= n_unf or n_inf == 0:
            continue
        ref["n_qualifying_operons"] += 1

        # Records which organisms carry each informative descriptor.
        inf_descs = {cln[r["feature_id"]] for r in inf_recs}
        for d in inf_descs:
            ref["present"][d].add(org)

        # Co-operon channel: every unique informative descriptor pair.
        dl = sorted(inf_descs)
        for i in range(len(dl)):
            for j in range(i + 1, len(dl)):
                k = (dl[i], dl[j])
                ref["coop_org"][k].add(org)
                ref["coop_inst"][k] += 1
                ref["deg_op"][dl[i]] += 1
                ref["deg_op"][dl[j]] += 1
                ref["M_op"] += 1

        # Adjacency channel: consecutive informative genes on the same contig.
        for a, b in zip(recs[:-1], recs[1:]):
            fa, fb = a["feature_id"], b["feature_id"]
            if contig.get((org, fa)) != contig.get((org, fb)):
                continue
            if uninf[fa] or uninf[fb]:
                continue
            da, db = cln[fa], cln[fb]
            if not da or not db or da == db:
                continue
            k = _key(da, db)
            ref["adj_org"][k].add(org)
            ref["adj_inst"][k] += 1
            ref["deg_adj"][da] += 1
            ref["deg_adj"][db] += 1
            ref["M_adj"] += 1


def update_reference(ref, genes, run_root, organisms=None, skip_existing=True):
    """Folds organisms into the reference in place and returns it un-finalized.

    organisms defaults to every organism in genes not yet in ref; organisms
    already present are skipped, or raise when skip_existing is False."""
    all_orgs = sorted(genes["organism"].unique())
    if organisms is not None:
        all_orgs = [o for o in all_orgs if o in set(organisms)]
    to_add = []
    for org in all_orgs:
        if org in ref["organisms_added"]:
            if skip_existing:
                continue
            raise ValueError("organism already in OCC reference: %s" % org)
        to_add.append(org)
    if not to_add:
        return ref

    contig = build_contig_map(run_root, to_add)
    gg = genes[genes["organism"].isin(to_add)]
    for org in to_add:
        _accumulate_organism(ref, org, gg[gg["organism"] == org], contig)
        ref["organisms_added"].add(org)
    ref["finalized"] = False
    return ref


def finalize_reference(ref):
    """Derives rho_adj / rho_op from the counts and returns ref.

    rho = Jeffreys lower bound BetaInv(delta; k+1/2, copres-k+1/2) times the
    enrichment weight lift/(lift+lambda), with lift = observed / (deg_a*deg_b/2M)."""
    delta = ref["params"]["delta"]
    lam = ref["params"]["lam"]
    use_w = ref["params"]["use_enrichment"]
    present = ref["present"]

    _lb_cache = {}

    def _lb(k, n):
        """Returns the memoised Jeffreys lower bound for k of n (scipy.stats.beta)."""
        if n <= 0:
            return 0.0
        k = min(k, n)
        ck = (k, n)
        v = _lb_cache.get(ck)
        if v is None:
            v = float(_beta.ppf(delta, k + 0.5, n - k + 0.5))
            _lb_cache[ck] = v
        return v

    def _w(inst, da, db, deg, M):
        """Returns the not-by-chance weight from the configuration-model lift."""
        if not use_w or M <= 0:
            return 1.0
        E = deg[da] * deg[db] / (2.0 * M)
        if E <= 0:
            return 1.0
        lift = inst / E
        return lift / (lift + lam)

    rho_adj = {}
    for k, orgs in ref["adj_org"].items():
        a, b = k
        copres = len(present[a] & present[b])
        rho_adj[k] = _lb(len(orgs), copres) * _w(ref["adj_inst"][k], a, b,
                                                 ref["deg_adj"], ref["M_adj"])
    rho_op = {}
    for k, orgs in ref["coop_org"].items():
        a, b = k
        copres = len(present[a] & present[b])
        rho_op[k] = _lb(len(orgs), copres) * _w(ref["coop_inst"][k], a, b,
                                                ref["deg_op"], ref["M_op"])
    ref["rho_adj"] = rho_adj
    ref["rho_op"] = rho_op
    ref["finalized"] = True
    return ref


def save_reference(ref, path):
    """Pickles the reference (counts and derived tables) to path."""
    import pickle
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(ref, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def load_reference(path):
    """Loads a pickled reference and checks its schema version."""
    import pickle
    with open(path, "rb") as fh:
        ref = pickle.load(fh)
    if ref.get("version") != REF_VERSION:
        raise ValueError("OCC reference version %r != expected %r; rebuild it"
                         % (ref.get("version"), REF_VERSION))
    return ref


def build_reference(genes, run_root, params=None):
    """Builds and finalizes a reference from every organism in genes in one go."""
    ref = new_reference(params)
    update_reference(ref, genes, run_root)
    finalize_reference(ref)
    return ref


def rho_adj(a, b, ref):
    """Returns the adjacency-channel reliability of a descriptor pair (0 if unseen)."""
    if not a or not b or a == b:
        return 0.0
    return ref["rho_adj"].get(_key(_norm(a), _norm(b)), 0.0)


def rho_op(a, b, ref):
    """Returns the co-operon-channel reliability of a descriptor pair (0 if unseen)."""
    if not a or not b or a == b:
        return 0.0
    return ref["rho_op"].get(_key(_norm(a), _norm(b)), 0.0)


def occ_for_gene(descriptor, adjacent_descs, cooperon_descs, ref, detail=False):
    """Returns the noisy-OR OCC, 1 - prod(1 - rho), over a gene's operon partners.

    adjacent_descs are immediate neighbours, cooperon_descs the other co-members;
    detail=True returns a dict with partner count and best partner/rho/channel.
    """
    d = _norm(descriptor)
    if not d:
        return (dict(occ=float("nan"), n_partners=0, best_partner="",
                     best_rho=0.0, best_channel="") if detail else float("nan"))
    floor = ref["params"]["link_floor"]
    prod = 1.0
    used = 0
    best_rho = 0.0
    best_partner = ""
    best_channel = ""
    for n in adjacent_descs:
        nn = _norm(n)
        if not nn or nn == d:
            continue
        r = rho_adj(d, nn, ref)
        if r >= floor and r > 0.0:
            prod *= (1.0 - r)
            used += 1
            if r > best_rho:
                best_rho, best_partner, best_channel = r, nn, "adj"
    for c in cooperon_descs:
        cc = _norm(c)
        if not cc or cc == d:
            continue
        r = rho_op(d, cc, ref)
        if r >= floor and r > 0.0:
            prod *= (1.0 - r)
            used += 1
            if r > best_rho:
                best_rho, best_partner, best_channel = r, cc, "op"
    occ = 0.0 if used == 0 else 1.0 - prod
    if detail:
        return dict(occ=occ, n_partners=used, best_partner=best_partner,
                    best_rho=best_rho, best_channel=best_channel)
    return occ


def compute_all_genes(genes, run_root, params=None, ref=None):
    """Returns a per-gene OCC DataFrame for every gene in a qualifying operon.

    The median informative OCC is stored in .attrs["occ0"]; ref is built if not given.
    """
    import pandas as pd
    if ref is None:
        ref = build_reference(genes, run_root, params)
    elif not ref.get("finalized"):
        finalize_reference(ref)

    g = genes.copy()
    organisms = sorted(g["organism"].unique())
    contig = build_contig_map(run_root, organisms)
    g["contig"] = [contig.get((o, f)) for o, f in zip(g["organism"], g["feature_id"])]
    uninf = {(o, f): bool(u) for o, f, u in
             zip(g["organism"], g["feature_id"], g["uninformative"])}
    cln = {(o, f): _norm(c) for o, f, c in
           zip(g["organism"], g["feature_id"], g["clean_descriptor"])}

    rows = []
    for (org, oid), sub in g.groupby(["organism", "operon_id"]):
        if not str(oid).startswith("operon_"):
            continue
        recs = sub.sort_values("start").to_dict("records")
        inf_recs = [r for r in recs if not uninf[(org, r["feature_id"])]
                    and cln[(org, r["feature_id"])]]
        n_inf = len(inf_recs)
        n_unf = len(recs) - n_inf
        if n_inf <= n_unf or n_inf == 0:
            continue  # operon does not qualify -> OCC undefined

        member_occ = {}
        for idx, r in enumerate(recs):
            fid = r["feature_id"]
            d = cln[(org, fid)]
            if uninf[(org, fid)] or not d:
                continue
            # Immediate informative neighbours on the same contig.
            adj = []
            for nb in (recs[idx - 1] if idx > 0 else None,
                       recs[idx + 1] if idx + 1 < len(recs) else None):
                if nb is None:
                    continue
                if nb["contig"] != r["contig"]:
                    continue
                nd = cln[(org, nb["feature_id"])]
                if not uninf[(org, nb["feature_id"])] and nd and nd != d:
                    adj.append(nd)
            adj_set = set(adj)
            coop = [cln[(org, rr["feature_id"])] for k, rr in enumerate(recs)
                    if k != idx and not uninf[(org, rr["feature_id"])]
                    and cln[(org, rr["feature_id"])]
                    and cln[(org, rr["feature_id"])] != d
                    and cln[(org, rr["feature_id"])] not in adj_set]
            det = occ_for_gene(d, adj, coop, ref, detail=True)
            member_occ[fid] = det["occ"]
            rows.append(dict(organism=org, feature_id=fid, clean_descriptor=d,
                             operon_id=oid, uninformative=False,
                             n_inf_context=n_inf - 1, n_partners=det["n_partners"],
                             best_partner=det["best_partner"], best_rho=det["best_rho"],
                             best_channel=det["best_channel"], occ=det["occ"]))

        # Uninformative members inherit the mean OCC of the informative members.
        coherence = float(np.mean(list(member_occ.values()))) if member_occ else 0.0
        for r in recs:
            fid = r["feature_id"]
            if uninf[(org, fid)]:
                rows.append(dict(organism=org, feature_id=fid,
                                 clean_descriptor=cln[(org, fid)],
                                 operon_id=oid, uninformative=True,
                                 n_inf_context=n_inf, n_partners=0,
                                 best_partner="", best_rho=0.0,
                                 best_channel="inherited", occ=coherence))

    df = pd.DataFrame(rows)
    if len(df):
        df.attrs["occ0"] = float(df.loc[~df["uninformative"], "occ"].median())
    df.attrs["ref"] = ref
    return df


def apply_occ(base_prob, occ, occ0, beta=1.0):
    """Returns base_prob shifted in log-odds by beta * (occ - occ0)."""
    p = min(max(float(base_prob), 1e-6), 1 - 1e-6)
    logit = np.log(p / (1 - p)) + beta * (float(occ) - float(occ0))
    return 1.0 / (1.0 + np.exp(-logit))
