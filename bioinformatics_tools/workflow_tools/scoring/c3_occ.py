"""c3_occ.py - Operon Context Confidence (OCC) engine.

A dynamic, per-candidate-gene factor that reads a gene's operon neighbourhood and
returns an INDEPENDENT reliability score in [0, 1], derived purely from
pan-genome operon co-occurrence.  It does NOT use C1/C2/C4.

------------------------------------------------------------------------------
DERIVATION (granular, biology-faithful)
------------------------------------------------------------------------------
Gene identity = clean_descriptor (functional name), lower-cased.  Everything is
organism-agnostic: statistics are pooled across all genomes by descriptor.

GATE - which operons provide valid context.
  For an operon with n_inf informative and n_unf uninformative members, the
  operon is CONSIDERED iff  n_inf > n_unf  (strict informative majority).
  * all-hypothetical operons are dropped (they would confuse the statistics);
  * informative-minority / tie operons are dropped;
  * uninformative members inside a qualifying operon are allowed but contribute
    NO positive evidence (they carry no descriptor statistics).
  The gate is applied BOTH when building the pan-genome statistics and when
  scoring a candidate.

PAN-GENOME PRIMITIVES (built only from qualifying operons), per descriptor pair
(a, b), a != b:
  present(d)   = # genomes where d is an informative member of a qualifying operon
  copres(a,b)  = |present(a) & present(b)|          (both functions available)
  adj(a,b)     = # genomes where a,b are IMMEDIATE operon neighbours
  coop(a,b)    = # genomes where a,b share an operon (any distance)
  => always  adj <= coop <= copres.

STEP 1 - recurrence-aware conditional co-occurrence (the "probability").
  Raw adj/copres can be a deceptive 1.0 at tiny support; recurrence is what
  earns trust.  So use the Jeffreys lower confidence bound:
      pi_adj(a,b) = BetaInv(delta; adj + 1/2, copres - adj + 1/2)
      pi_op (a,b) = BetaInv(delta; coop + 1/2, copres - coop + 1/2)
  delta = 0.05 (a 95% lower bound).  Behaviour: 20/20 -> 0.91, 2/2 -> 0.43,
  1/1 -> 0.23, 1/15 -> 0.01.  Deep, conserved partnerships score high; thin or
  coincidental ones collapse to ~0.

STEP 2 - enrichment safeguard (not-by-chance), from a configuration-null model.
      E(a,b)   = deg(a) * deg(b) / (2 M)     (expected co-occurrences by chance)
      lift     = observed(a,b) / E(a,b)
      w(a,b)   = lift / (lift + lambda)      (lambda = 1)
  For genuine partnerships lift is enormous so w ~ 1; only near-chance links are
  damped.  Computed separately for the adjacency and co-operon channels.

STEP 3 - link reliability (per pair, two channels):
      rho_adj(a,b) = pi_adj(a,b) * w_adj(a,b)
      rho_op (a,b) = pi_op (a,b) * w_op (a,b)

STEP 4 - per-gene aggregation (noisy-OR over the neighbourhood).
  For candidate g with informative immediate neighbours N (adjacency channel) and
  non-adjacent informative co-members C (co-operon channel):
      OCC(g) = 1 - PROD_{n in N} (1 - rho_adj(g, n))
                 * PROD_{c in C} (1 - rho_op (g, c))
  Semantics: the confidence that g's operon placement is corroborated by at least
  one genuinely conserved partnership, compounded over all partners.  A single
  rock-solid partner drives OCC -> 1; coincidental links (rho ~ 0) leave it
  untouched; an isolated / all-hypothetical-neighbour gene gets OCC = 0.

  Uninformative candidate embedded in a qualifying operon: it has no descriptor,
  so it INHERITS the operon's coherence = mean OCC of the informative members
  (guilt-by-association with a demonstrably real operon).
------------------------------------------------------------------------------
DYNAMIC REFERENCE DATABASE (incremental, per newly-labeled organism)
------------------------------------------------------------------------------
The reference stores only ADDITIVE pan-genome counts (present-sets, adj/coop
organism-sets and instance counts, degrees, M totals, #qualifying operons).
None of these depend on any cross-organism ordering, so a freshly labeled
organism folds in WITHOUT reprocessing the ones already in the database:

    ref = new_reference()                 # or load_reference(path)
    update_reference(ref, genes, run_root)          # add all new organisms
    #   ... later, when organism #22 is labeled ...
    update_reference(ref, genes22, run_root, organisms=["org22"])
    finalize_reference(ref)               # derive rho_adj / rho_op ONCE
    save_reference(ref, path)

  * update_reference is CHEAP (counts only, ~1.2 s/organism) and idempotent
    (organisms_added guards against double-counting).
  * finalize_reference derives the per-pair reliabilities from the counts; it
    is a PURE function of the counts and is memoised over the few distinct
    (k, n) integer pairs, so it is ~0.1 s regardless of database size.
  * PATTERN: accumulate organisms with update_reference, then finalize ONCE
    just before scoring - do not finalize after every organism.
  * EXACTNESS: folding organisms in one-at-a-time yields bit-for-bit identical
    reliabilities to a from-scratch build_reference() over all of them (adding
    an organism changes M and the copres of touched descriptors, so rho shifts
    globally - but finalize recomputes every pair, keeping it exact).
  * build_reference() is the batch convenience wrapper (new + update-all +
    finalize).
------------------------------------------------------------------------------
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
    return (d or "").strip().lower()


def build_contig_map(run_root, organisms):
    """(organism, feature_id) -> contig id, parsed from gene_id."""
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
    return (a, b) if a < b else (b, a)


REF_VERSION = 2   # bump when the counts schema changes


def new_reference(params=None):
    """Create an empty, dynamic OCC reference (a growable "database").

    Holds only the ADDITIVE pan-genome counts; the per-pair reliabilities
    (rho_adj / rho_op) are derived from them by finalize_reference().  Because
    every field is additive across organisms, a newly labeled organism can be
    folded in with update_reference() WITHOUT reprocessing the existing ones."""
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
    """Fold ONE organism's qualifying operons into ref's additive counts."""
    uninf = {f: bool(u) for f, u in
             zip(org_genes["feature_id"], org_genes["uninformative"])}
    cln = {f: _norm(c) for f, c in
           zip(org_genes["feature_id"], org_genes["clean_descriptor"])}

    for oid, sub in org_genes.groupby("operon_id"):
        if not str(oid).startswith("operon_"):
            continue
        recs = sub.sort_values("start").to_dict("records")
        # informative-majority gate
        inf_recs = [r for r in recs if not uninf[r["feature_id"]]
                    and cln[r["feature_id"]]]
        n_inf = len(inf_recs)
        n_unf = len(recs) - n_inf
        if n_inf <= n_unf or n_inf == 0:
            continue
        ref["n_qualifying_operons"] += 1

        # present (informative members of a qualifying operon)
        inf_descs = {cln[r["feature_id"]] for r in inf_recs}
        for d in inf_descs:
            ref["present"][d].add(org)

        # --- co-operon channel: unique informative descriptor pairs ---------
        dl = sorted(inf_descs)
        for i in range(len(dl)):
            for j in range(i + 1, len(dl)):
                k = (dl[i], dl[j])
                ref["coop_org"][k].add(org)
                ref["coop_inst"][k] += 1
                ref["deg_op"][dl[i]] += 1
                ref["deg_op"][dl[j]] += 1
                ref["M_op"] += 1

        # --- adjacency channel: consecutive informative genes (same contig) -
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
    """Fold one or more organisms into the reference IN PLACE (the dynamic path).

    genes       : DataFrame with at least the organisms to add.
    organisms   : which organisms to add; default = every organism in ``genes``
                  not already present in the reference.
    skip_existing : silently skip organisms already added (idempotent); set
                  False to raise instead.
    Leaves the reference un-finalized (rho tables stale) -> call
    finalize_reference() before scoring.  Returns ref."""
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


def _copy_counts(ref):
    """Deep copy of the additive count structures only (derived rho tables are
    left empty for finalize_reference to recompute).  Used by leave-one-out so
    subtraction never mutates the shared/on-disk reference."""
    return dict(
        version=ref["version"],
        present=defaultdict(set, {d: set(s) for d, s in ref["present"].items()}),
        adj_org=defaultdict(set, {k: set(s) for k, s in ref["adj_org"].items()}),
        coop_org=defaultdict(set, {k: set(s) for k, s in ref["coop_org"].items()}),
        adj_inst=defaultdict(int, dict(ref["adj_inst"])),
        coop_inst=defaultdict(int, dict(ref["coop_inst"])),
        deg_adj=defaultdict(int, dict(ref["deg_adj"])),
        deg_op=defaultdict(int, dict(ref["deg_op"])),
        M_adj=ref["M_adj"], M_op=ref["M_op"],
        n_qualifying_operons=ref["n_qualifying_operons"],
        organisms_added=set(ref["organisms_added"]),
        rho_adj={}, rho_op={},
        params=dict(ref["params"]), finalized=False,
    )


def subtract_organism(ref, genes, run_root, organism, token=None,
                      require_present=True):
    """Return a COPY of ``ref`` with one genome's contribution removed - the
    exact inverse of update_reference() for one organism (leave-one-out).

    Every field the reference stores is ADDITIVE across organisms, so removing a
    genome is exact subtraction: its qualifying-operon contribution is recomputed
    in isolation (identical code path to how it was added, reading the genome's
    own files under ``organism`` in ``run_root``) and taken back out of the
    present-/adjacency-/co-operon-sets, the instance/degree counts and the M
    totals.  Keys that drop to an empty set or zero count are pruned, so the
    result is bit-for-bit identical to a from-scratch build over the OTHER
    genomes (verified).

    ``organism`` names the genome to read/recompute; ``token`` is the identifier
    it is stored under in the reference's membership sets (default: ``organism``).
    They differ only when a genome is matched to the reference by CONTENT
    (fingerprint) but was stored under a different name - then read it under its
    current name yet subtract the stored token.

    Leaves the copy un-finalized; call finalize_reference() before scoring.
    Raises if ``token`` is not in the reference (unless require_present)."""
    tok = token if token is not None else organism
    if tok not in ref["organisms_added"]:
        if require_present:
            raise ValueError("organism not in OCC reference: %s" % tok)
        return _copy_counts(ref)

    # recompute this genome's own contribution in isolation (read by its name)
    g = new_reference(ref["params"])
    contig = build_contig_map(run_root, [organism])
    _accumulate_organism(g, organism, genes[genes["organism"] == organism], contig)

    out = _copy_counts(ref)
    # membership sets: g holds exactly one token (`organism`); remove the STORED
    # token `tok` from every set the genome touched.
    for d in g["present"]:
        if d in out["present"]:
            out["present"][d].discard(tok)
            if not out["present"][d]:
                del out["present"][d]
    for bucket in ("adj_org", "coop_org"):
        for k in g[bucket]:
            if k in out[bucket]:
                out[bucket][k].discard(tok)
                if not out[bucket][k]:
                    del out[bucket][k]
    # instance / degree / M counts: token-independent magnitudes, subtract them
    for bucket in ("adj_inst", "coop_inst", "deg_adj", "deg_op"):
        for k, c in g[bucket].items():
            out[bucket][k] -= c
            if out[bucket][k] <= 0:
                out[bucket].pop(k, None)
    out["M_adj"] -= g["M_adj"]
    out["M_op"] -= g["M_op"]
    out["n_qualifying_operons"] -= g["n_qualifying_operons"]
    out["organisms_added"].discard(tok)
    out["finalized"] = False
    return out


def leave_one_out(ref, genes, run_root, organism, token=None,
                  require_present=True):
    """Convenience: subtract_organism() + finalize_reference().  Returns a
    finalized copy of ``ref`` as if the genome had never been added - the
    reference to score it against for a genuine out-of-sample score."""
    loo = subtract_organism(ref, genes, run_root, organism, token=token,
                            require_present=require_present)
    finalize_reference(loo)
    return loo


def finalize_reference(ref):
    """Derive the per-pair reliabilities rho_adj / rho_op from the current
    counts.  Pure function of the counts, so it can be re-run cheaply after
    every update_reference().  Returns ref."""
    delta = ref["params"]["delta"]
    lam = ref["params"]["lam"]
    use_w = ref["params"]["use_enrichment"]
    present = ref["present"]

    _lb_cache = {}

    def _lb(k, n):
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
    """Persist the reference database (counts + derived tables) to disk."""
    import pickle
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(ref, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def load_reference(path):
    """Load a reference database written by save_reference()."""
    import pickle
    with open(path, "rb") as fh:
        ref = pickle.load(fh)
    if ref.get("version") != REF_VERSION:
        raise ValueError("OCC reference version %r != expected %r; rebuild it"
                         % (ref.get("version"), REF_VERSION))
    return ref


# ---------------------------------------------------------------------------
# Members sidecar: genome content-fingerprint -> organism token.
#
# The reference keys organisms by an opaque token (their name today).  To do
# leave-one-out by CONTENT identity - and to dedupe a genome resubmitted under a
# different name - we keep a small sidecar next to the pickle mapping each
# added genome's fingerprint (c3_lib.genome_fingerprint) to the token it was
# stored under.  It is written alongside (never inside) the pickle, so the
# on-disk reference schema is untouched.
# ---------------------------------------------------------------------------
def members_sidecar_path(reference_path):
    from pathlib import Path
    return Path(str(reference_path) + ".members.tsv")


def load_members(reference_path):
    """Return {genome_fingerprint: organism_token} recorded next to the ref
    (empty dict if the sidecar does not exist)."""
    p = members_sidecar_path(reference_path)
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2 and row[0]:
                    out[row[0]] = row[1]
    return out


def record_member(reference_path, fingerprint, organism):
    """Idempotently upsert fingerprint->organism into the sidecar.  No-op for an
    empty fingerprint.  Caller should hold the reference lock (the OCC updater
    already does)."""
    if not fingerprint:
        return
    members = load_members(reference_path)
    if members.get(fingerprint) == organism:
        return
    members[fingerprint] = organism
    p = members_sidecar_path(reference_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for fp, org in sorted(members.items(), key=lambda x: (x[1], x[0])):
            w.writerow([fp, org])


# ---------------------------------------------------------------------------
# Pool-stats sidecar: per-genome descriptive counts (total genes, operonic vs
# singleton genes, informative vs uninformative operons) for figure provenance.
# The OCC itself keeps only qualifying-operon descriptor statistics, so these
# whole-genome tallies live alongside it. Keyed by organism token; carries the
# content fingerprint too. Being per-genome, figures can aggregate over any
# subset -- e.g. the pool MINUS the organism being reported (leave-one-out).
# ---------------------------------------------------------------------------
POOL_STAT_COLS = ("total_genes", "operonic_genes", "singleton_genes",
                  "n_operons", "n_informative_operons", "n_uninformative_operons")


def pool_stats_sidecar_path(reference_path):
    from pathlib import Path
    return Path(str(reference_path) + ".genome_stats.tsv")


def load_pool_stats(reference_path):
    """Return {organism_token: {stat: int, ...}} from the pool-stats sidecar
    (empty dict if absent). Each value has the POOL_STAT_COLS keys."""
    p = pool_stats_sidecar_path(reference_path)
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                org = (row.get("organism") or "").strip()
                if not org:
                    continue
                out[org] = {c: int(float(row.get(c, 0) or 0)) for c in POOL_STAT_COLS}
    return out


def record_pool_stats(reference_path, fingerprint, organism, stats):
    """Idempotently upsert one genome's pool stats into the sidecar. ``stats`` is
    a dict over POOL_STAT_COLS (e.g. c3_lib.genome_pool_stats())."""
    if not organism:
        return
    p = pool_stats_sidecar_path(reference_path)
    existing = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                o = (row.get("organism") or "").strip()
                if o:
                    existing[o] = row
    existing[organism] = dict(fingerprint=fingerprint or "", organism=organism,
                              **{c: int(stats.get(c, 0)) for c in POOL_STAT_COLS})
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(("fingerprint", "organism") + POOL_STAT_COLS)
        for o in sorted(existing):
            r = existing[o]
            w.writerow([r.get("fingerprint", ""), o]
                       + [int(float(r.get(c, 0) or 0)) for c in POOL_STAT_COLS])


def build_reference(genes, run_root, params=None):
    """Convenience batch build = new + update(all organisms) + finalize.

    Produces exactly the same reliabilities as folding the organisms in one at a
    time via update_reference(); use the incremental path for the dynamic
    database, this one for a from-scratch build."""
    ref = new_reference(params)
    update_reference(ref, genes, run_root)
    finalize_reference(ref)
    return ref


def rho_adj(a, b, ref):
    if not a or not b or a == b:
        return 0.0
    return ref["rho_adj"].get(_key(_norm(a), _norm(b)), 0.0)


def rho_op(a, b, ref):
    if not a or not b or a == b:
        return 0.0
    return ref["rho_op"].get(_key(_norm(a), _norm(b)), 0.0)


def occ_for_gene(descriptor, adjacent_descs, cooperon_descs, ref, detail=False):
    """Dynamic OCC for one informative candidate gene.

    descriptor      : candidate clean_descriptor (informative).
    adjacent_descs  : informative descriptors of the candidate's IMMEDIATE
                      operon neighbours (in the genome being scored).
    cooperon_descs  : informative descriptors of the candidate's NON-adjacent
                      operon co-members.
    Returns OCC in [0, 1] (or, if detail=True, a dict with the OCC plus the
    partner count, best partner and best link reliability).
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
    """Compute OCC for every gene that sits in a qualifying operon.

    Returns a DataFrame: organism, feature_id, clean_descriptor, operon_id,
    uninformative, n_inf_context, n_partners, best_partner, best_rho,
    best_channel, occ, plus a fitted neutral pivot in .attrs["occ0"].
    A prebuilt ``ref`` (from build_reference) may be supplied to avoid rebuilding.
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

        # position index for adjacency lookup
        member_occ = {}
        for idx, r in enumerate(recs):
            fid = r["feature_id"]
            d = cln[(org, fid)]
            if uninf[(org, fid)] or not d:
                continue
            # immediate informative neighbours
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

        # uninformative members inherit operon coherence = mean informative OCC
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
    """Recommended log-odds application: shift base_prob by OCC around pivot occ0."""
    p = min(max(float(base_prob), 1e-6), 1 - 1e-6)
    logit = np.log(p / (1 - p)) + beta * (float(occ) - float(occ0))
    return 1.0 / (1.0 + np.exp(-logit))
