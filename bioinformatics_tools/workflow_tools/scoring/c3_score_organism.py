#!/usr/bin/env python3
"""
C3 organism scorer -- Operon Context Confidence (OCC).

Scores every gene of one organism by how well its operon neighbourhood recurs
across the pan-genome reference (leave-one-out by default).
  adjacency mode: C3 = geometric mean of rho_adj over the operon's UniOP-linked
                  adjacencies with OCC evidence; rho_adj is the Jeffreys posterior
                  mean (k+0.5)/(n+1) times the enrichment weight.
  hybrid mode:    per-gene max of best adjacent rho_adj and best co-member rho_op.
Pairs with a hypothetical member or no OCC evidence are excluded; an operon with
none left scores 0. Singletons, non-coding genes and all-hypothetical operons
score 0.5 (--allhyp-value). Both variants are written side by side.

Usage:
    python3 c3_score_organism.py \
        --operon-info    /path/to/labeled-genes-operon-info.tsv \
        --genes-file     /path/to/labeled-genes.tsv \
        --operon-results /path/to/operon/operon_results.tsv \
        --reference      /path/to/reference_data/occ_reference.pkl \
        --output         /path/to/scored-...-c3-operonic-context-confidence.tsv \
        [--allhyp-value 0.5]
"""
import argparse
import csv
import math
import sys
from pathlib import Path

# c3_occ and c3_lib sit next to this script and match how the reference was built.
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import c3_occ            # noqa: E402
import c3_lib as L       # noqa: E402

csv.field_size_limit(10_000_000)

_NEUTRAL = 0.5


# ---------------------------------------------------------------------------
# Adjacency channel: Jeffreys posterior mean here (c3_score); c3_occ.rho_adj
# gives the 95% lower bound (c3_lowerbound).
# ---------------------------------------------------------------------------
def _enrich_w(key, na, nb, ref):
    """Returns the adjacency enrichment weight lift/(lift+lambda), as in c3_occ.finalize_reference."""
    p = ref["params"]
    if not p.get("use_enrichment", True):
        return 1.0
    M = ref["M_adj"]
    if M <= 0:
        return 1.0
    E = ref["deg_adj"][na] * ref["deg_adj"][nb] / (2.0 * M)
    if E <= 0:
        return 1.0
    lift = ref["adj_inst"][key] / E
    lam = p.get("lam", 1.0)
    return lift / (lift + lam)


def rho_adj_mean(a, b, ref):
    """Returns the Jeffreys posterior mean (k+0.5)/(n+1) of P(adjacent | both present) times the enrichment weight.

    0.0 for empty, identical or never-adjacent descriptors.
    """
    na, nb = c3_occ._norm(a), c3_occ._norm(b)
    if not na or not nb or na == nb:
        return 0.0
    key = c3_occ._key(na, nb)
    orgs = ref["adj_org"].get(key)
    if not orgs:
        return 0.0
    if na not in ref["present"] or nb not in ref["present"]:
        return 0.0
    n = len(ref["present"][na] & ref["present"][nb])
    if n <= 0:
        return 0.0
    k = min(len(orgs), n)
    mean = (k + 0.5) / (n + 1.0)
    return mean * _enrich_w(key, na, nb, ref)


def enrich_weight(a, b, ref):
    """Returns the adjacency enrichment weight alone, 0.0 for never-adjacent pairs."""
    na, nb = c3_occ._norm(a), c3_occ._norm(b)
    if not na or not nb or na == nb:
        return 0.0
    key = c3_occ._key(na, nb)
    if not ref["adj_org"].get(key):
        return 0.0
    return _enrich_w(key, na, nb, ref)


# ---------------------------------------------------------------------------
# Membership channel (co-operon, order-free), used by the hybrid C3.
# ---------------------------------------------------------------------------
def _enrich_w_op(key, na, nb, ref):
    """Returns the co-operon enrichment weight lift/(lift+lambda)."""
    p = ref["params"]
    if not p.get("use_enrichment", True):
        return 1.0
    M = ref.get("M_op", 0)
    if M <= 0:
        return 1.0
    E = ref["deg_op"][na] * ref["deg_op"][nb] / (2.0 * M)
    if E <= 0:
        return 1.0
    lift = ref["coop_inst"][key] / E
    return lift / (lift + p.get("lam", 1.0))


def rho_op_mean(a, b, ref):
    """Returns the Jeffreys posterior mean (k+0.5)/(n+1) of P(share an operon | both present) times the enrichment weight."""
    na, nb = c3_occ._norm(a), c3_occ._norm(b)
    if not na or not nb or na == nb:
        return 0.0
    key = c3_occ._key(na, nb)
    orgs = ref.get("coop_org", {}).get(key)
    if not orgs:
        return 0.0
    if na not in ref["present"] or nb not in ref["present"]:
        return 0.0
    n = len(ref["present"][na] & ref["present"][nb])
    if n <= 0:
        return 0.0
    k = min(len(orgs), n)
    return (k + 0.5) / (n + 1.0) * _enrich_w_op(key, na, nb, ref)


def rho_op_lb(a, b, ref):
    """Returns the finalized Jeffreys 95% lower bound for the membership channel."""
    na, nb = c3_occ._norm(a), c3_occ._norm(b)
    if not na or not nb or na == nb:
        return 0.0
    return ref.get("rho_op", {}).get(c3_occ._key(na, nb), 0.0)


# ---------------------------------------------------------------------------
def geomean(vals):
    """Returns the geometric mean via math.log/exp; 0.0 if any term <= 0, None if empty."""
    vals = list(vals)
    if not vals:
        return None
    if any(v <= 0.0 for v in vals):
        return 0.0
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def is_informative(uninformative, clean_descriptor):
    """Returns True when the gene is not uninformative and has a descriptor."""
    if bool(uninformative):
        return False
    return bool(str(clean_descriptor).strip())


# ---------------------------------------------------------------------------
# Conflict signal: a mate's consensus partner (>= _CONFLICT_MIN_SUPPORT genomes)
# differs from this descriptor, which has <= _CONFLICT_MINORITY of that support.
# Reported as c3_descriptor_conflict; not applied to the score.
_CONFLICT_MIN_SUPPORT = 3
_CONFLICT_MINORITY = 0.34


def build_partner_consensus(ref):
    """Returns {descriptor: Counter{co-operon partner: #genomes}} from the order-free coop_org counts (collections.Counter)."""
    from collections import Counter, defaultdict
    partners = defaultdict(Counter)
    for key, orgs in ref.get("coop_org", ref.get("adj_org", {})).items():
        try:
            a, b = key
        except (TypeError, ValueError):
            continue
        n = len(orgs) if hasattr(orgs, "__len__") else int(orgs)
        partners[a][b] += n
        partners[b][a] += n
    return partners


def descriptor_conflict(descriptor, neigh_descriptors, partners):
    """Returns the strongest contradiction (0-1) of a descriptor by its mates' consensus partners."""
    worst = 0.0
    for a in neigh_descriptors:
        pc = partners.get(a)
        if not pc:
            continue
        cons_descriptor, cons_n = pc.most_common(1)[0]
        our_n = pc.get(descriptor, 0)
        if (cons_n >= _CONFLICT_MIN_SUPPORT and cons_descriptor != descriptor
                and our_n <= _CONFLICT_MINORITY * cons_n):
            worst = max(worst, (cons_n - our_n) / cons_n)
    return worst


def load_pairwise_probs(operon_results_path):
    """Reads operon_results.tsv with csv into {frozenset{feature_id, neighbour_id}: UniOP pair probability}."""
    pmap = {}
    p = Path(operon_results_path)
    if not p.is_file():
        print(f"[c3-scorer] WARNING: operon_results not found: {p} "
              f"(all pairs will be treated as unlinked)", file=sys.stderr)
        return pmap
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            for gcol, pcol in (
                ("OPERON_downstream_gene_id", "OPERON_downstream_pairwise_probability"),
                ("OPERON_upstream_gene_id", "OPERON_upstream_pairwise_probability"),
            ):
                other = row.get(gcol, "")
                praw = row.get(pcol, "")
                if fid and other and praw not in ("", "NA", "None", "nan"):
                    try:
                        pmap[frozenset((fid, other))] = float(praw)
                    except ValueError:
                        pass
    return pmap


# ---------------------------------------------------------------------------
def compute(genes, ref, pmap, allhyp_value=_NEUTRAL, c3_mode="adjacency"):
    """Scores every gene with pandas groupby per operon and returns {feature_id: {c3, c3_lowerbound, ...}}.

    Computes both the adjacency and hybrid C3; c3_mode picks which is the primary c3.
    In hybrid mode uninformative members inherit their operon's mean.
    """
    per_gene = {}

    g = genes.copy()
    g["operon_id"] = g["operon_id"].astype(str)
    g["inf_flag"] = [
        is_informative(u, d)
        for u, d in zip(g["uninformative"], g["clean_descriptor"])
    ]
    partners = build_partner_consensus(ref)   # cross-genome label consensus (for conflict)

    # ---- singletons (no operon_ id: NOT_IN_AN_OPERON / non-coding) ----------
    sing = g[~g["operon_id"].str.startswith("operon_")]
    for r in sing.itertuples(index=False):
        per_gene[r.feature_id] = dict(
            c3=_NEUTRAL, c3_lowerbound=_NEUTRAL, gene_class="singleton",
            c3_adjacency=_NEUTRAL, c3_lb_adjacency=_NEUTRAL,
            c3_hybrid=_NEUTRAL, c3_lb_hybrid=_NEUTRAL,
            operon_id="", n_pairs=0, n_pairs_used=0, n_supported=0,
            n_neutral=0, n_nolink=0, descriptor_conflict=0.0, operon_significance=0.0,
            operon_ambiguity=0.0,
        )

    # ---- operons ------------------------------------------------------------
    opg = g[g["operon_id"].str.startswith("operon_")]
    for oid, members in opg.groupby("operon_id"):
        members = members.sort_values("start")
        recs = list(members.itertuples(index=False))
        n_inf = sum(1 for x in recs if x.inf_flag)
        all_hyp = (n_inf == 0)

        # C3 is conservation only (UniOP probability enters later as C2); a
        # never-seen adjacency is excluded rather than scored.
        terms_mean, terms_lb, w_terms = [], [], []
        n_zero = n_neutral = n_supported = n_nolink = 0
        for a, b in zip(recs[:-1], recs[1:]):
            P = pmap.get(frozenset((a.feature_id, b.feature_id)))
            if P is None:                      # not UniOP-linked (contig break)
                n_nolink += 1
                continue
            if not a.inf_flag or not b.inf_flag:   # hypothetical -> cannot assess -> exclude
                n_neutral += 1
                continue
            rho_l = c3_occ.rho_adj(a.clean_descriptor, b.clean_descriptor, ref)
            rho_m = rho_adj_mean(a.clean_descriptor, b.clean_descriptor, ref)
            if rho_m > 0.0 or rho_l > 0.0:     # OCC has conservation evidence
                terms_mean.append(rho_m)       # PURE conservation (no P_ab)
                terms_lb.append(rho_l)
                w_terms.append(enrich_weight(a.clean_descriptor, b.clean_descriptor, ref))
                n_supported += 1
            else:                              # never-seen -> no evidence -> EXCLUDE (novel != penalty)
                n_zero += 1

        n_used = len(terms_mean)
        # Operon significance: geometric mean enrichment over supported adjacencies
        operon_sig = geomean(w_terms) if w_terms else 0.0
        if all_hyp:
            c3 = c3_lb = allhyp_value
            gclass = "all_hypothetical_operon"
        elif n_used == 0:                      # no OCC conservation evidence -> novel -> 0 (no boost)
            c3 = c3_lb = 0.0
            gclass = "operon_member"
        else:
            c3 = geomean(terms_mean)
            c3_lb = geomean(terms_lb)
            gclass = "operon_member"

        # Operon ambiguity: fraction of adjacent pairs with a hypothetical member (review note only)
        n_pairs_total = len(recs) - 1
        operon_ambiguity = (n_neutral / n_pairs_total) if n_pairs_total > 0 else 0.0

        # ---- Hybrid variant: per-gene max(best adjacency, best co-member) ------
        hyb_per, hyb_mean, hyb_mean_lb = {}, 0.0, 0.0
        if not all_hyp:
            inf_pos = [i for i, r in enumerate(recs)
                       if r.inf_flag and str(r.clean_descriptor).strip()]
            for i in inf_pos:
                d = recs[i].clean_descriptor
                bm = bl = 0.0
                for j in (i - 1, i + 1):                 # immediate informative neighbours
                    if 0 <= j < len(recs) and recs[j].inf_flag:
                        m = rho_adj_mean(d, recs[j].clean_descriptor, ref)
                        if m > bm:
                            bm, bl = m, c3_occ.rho_adj(d, recs[j].clean_descriptor, ref)
                for j in inf_pos:                        # co-members, any distance
                    if j != i:
                        m = rho_op_mean(d, recs[j].clean_descriptor, ref)
                        if m > bm:
                            bm, bl = m, rho_op_lb(d, recs[j].clean_descriptor, ref)
                hyb_per[i] = (bm, bl)
            if hyb_per:
                hyb_mean = sum(v[0] for v in hyb_per.values()) / len(hyb_per)
                hyb_mean_lb = sum(v[1] for v in hyb_per.values()) / len(hyb_per)

        for idx, x in enumerate(recs):
            # Conflict is checked against all operon mates (order-free).
            mates = [recs[j].clean_descriptor for j in range(len(recs)) if j != idx]
            conflict = (descriptor_conflict(x.clean_descriptor, mates, partners)
                        if x.inf_flag else 0.0)
            c3_adj, c3_adj_lb = c3, c3_lb                # adjacency: shared operon-level
            if all_hyp:
                c3_hyb, c3_hyb_lb = c3, c3_lb            # = allhyp_value (both equal)
            elif idx in hyb_per:
                c3_hyb, c3_hyb_lb = hyb_per[idx]         # informative: own best partner
            elif x.inf_flag:
                c3_hyb, c3_hyb_lb = 0.0, 0.0             # informative, no partner -> novel
            else:
                c3_hyb, c3_hyb_lb = hyb_mean, hyb_mean_lb  # uninformative inherits coherence
            g_c3, g_lb = ((c3_hyb, c3_hyb_lb) if c3_mode == "hybrid"
                          else (c3_adj, c3_adj_lb))       # primary variant
            per_gene[x.feature_id] = dict(
                c3=g_c3, c3_lowerbound=g_lb,
                c3_adjacency=c3_adj, c3_lb_adjacency=c3_adj_lb,
                c3_hybrid=c3_hyb, c3_lb_hybrid=c3_hyb_lb,
                gene_class=gclass, operon_id=oid,
                n_pairs=len(recs) - 1, n_pairs_used=n_used,
                n_supported=n_supported, n_neutral=(n_neutral + n_zero),
                n_nolink=n_nolink, descriptor_conflict=conflict,
                operon_significance=operon_sig, operon_ambiguity=operon_ambiguity,
                c3_mode=c3_mode,
            )
    return per_gene


# ---------------------------------------------------------------------------
def _breakdown_and_formula(info, raw_operon_id):
    """Returns the (reason, formula) text for one gene's C3."""
    cls = info["gene_class"]
    c3 = info["c3"]
    if cls == "singleton":
        if raw_operon_id == "NOT_APPLICABLE_NON_CODING":
            reason = "non-coding (neutral)"
        else:
            reason = "singleton, not in an operon (neutral)"
        return reason, f"neutral = {c3:.4f}"
    if cls == "all_hypothetical_operon":
        return ("all-hypothetical operon (neutral, equivalent to singleton)",
                f"neutral = {c3:.4f}")
    # operon_member
    if info.get("c3_mode") == "hybrid":
        if c3 <= 0.0:
            return ("operon member, no conserved adjacency or co-member partnership "
                    "(novel context, neutral)", f"max(rho_adj, rho_op) = {c3:.4f}")
        return ("operon context (hybrid): strongest of this gene's adjacency and "
                "co-membership partnerships across genomes",
                f"max(best rho_adj, best rho_op) = {c3:.4f}")
    if info["n_pairs_used"] == 0:
        return ("operon member, no UniOP-linked informative adjacency (neutral)",
                f"neutral = {c3:.4f}")
    reason = (f"operon context: {info['n_supported']}/{info['n_pairs_used']} "
              f"OCC-supported adjacencies, {info['n_neutral']} neutral"
              + (f", {info['n_nolink']} unlinked" if info["n_nolink"] else ""))
    formula = (f"geomean(rho_adj) over {info['n_pairs_used']} "
               f"adjacencies = {c3:.4f}")
    return reason, formula


def main():
    """Loads the OCC reference (optionally leave-one-out), scores the organism and writes the C3 table with csv."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--operon-info", required=True,
                    help="labeled-genes-operon-info.tsv")
    ap.add_argument("--genes-file", required=True, help="labeled-genes.tsv")
    ap.add_argument("--operon-results", required=True,
                    help="operon/operon_results.tsv (UniOP per-pair probabilities)")
    ap.add_argument("--reference", required=True,
                    help="OCC reference pickle (occ_reference.pkl)")
    ap.add_argument("--output", required=True, help="output TSV path")
    ap.add_argument("--allhyp-value", type=float, default=_NEUTRAL,
                    help="C3 for an all-hypothetical operon (default 0.5)")
    ap.add_argument("--c3-mode", choices=("adjacency", "hybrid"), default="adjacency",
                    help="'adjacency' (default, shipped): operon-level geomean(rho_adj). "
                         "'hybrid': per-gene max(best adjacency, best co-member) -- "
                         "order-free, robust to a hypothetical breaking the chain.")
    ap.add_argument("--leave-one-out", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="If the genome being scored is already in the OCC "
                         "reference (matched by content fingerprint, else name), "
                         "remove its own contribution before scoring so C3 is a "
                         "genuine out-of-sample estimate (default: on).")
    args = ap.parse_args()

    for p in (args.operon_info, args.genes_file, args.reference):
        if not Path(p).is_file():
            print(f"[c3-scorer] ERROR: not found: {p}", file=sys.stderr)
            sys.exit(1)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Organism name from the first labeled row
    organism = ""
    with open(args.genes_file, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            organism = (row.get("organism_name") or "").strip()
            break

    print(f"[c3-scorer] organism: {organism or '(unknown)'}")
    print("[c3-scorer] loading OCC reference ...")
    ref = c3_occ.load_reference(args.reference)
    if not ref.get("finalized"):
        c3_occ.finalize_reference(ref)
    print(f"  reference: {len(ref.get('organisms_added', []))} organisms, "
          f"{len(ref['rho_adj'])} rho_adj entries")

    print("[c3-scorer] joining labeled-genes + operon-info ...")
    # compute_hash=True adds aa_hash for the leave-one-out content fingerprint.
    genes = L.load_organism(organism, Path(args.genes_file),
                            Path(args.operon_info), compute_hash=True)
    print(f"  genes: {len(genes)}")

    # --- Leave-one-out: removes this genome from the reference when present,
    # matched by content fingerprint, else by name.
    if args.leave_one_out:
        run_root = Path(args.genes_file).resolve().parents[2]
        fp = L.genome_fingerprint(genes["aa_hash"])
        members = c3_occ.load_members(args.reference)
        token = members.get(fp) if fp else None
        if token is None and organism in ref.get("organisms_added", set()):
            token = organism  # fallback: matched by name
        if token is not None:
            n_before = len(ref["organisms_added"])
            ref = c3_occ.leave_one_out(ref, genes, run_root, organism, token=token)
            how = "fingerprint" if (fp and members.get(fp) == token) else "name"
            print(f"[c3-scorer] leave-one-out: removed this genome "
                  f"(matched by {how} as '{token}') from the reference "
                  f"-> {n_before} - 1 = {len(ref['organisms_added'])} organisms, "
                  f"{len(ref['rho_adj'])} rho_adj entries")
        else:
            print("[c3-scorer] leave-one-out: this genome is not in the reference "
                  "-> scoring against the full reference (already out-of-sample)")

    pmap = load_pairwise_probs(args.operon_results)
    print(f"  UniOP per-pair probabilities: {len(pmap)}")

    print(f"[c3-scorer] C3 mode: {args.c3_mode}")
    per_gene = compute(genes, ref, pmap, allhyp_value=args.allhyp_value,
                       c3_mode=args.c3_mode)

    # ---- one output row per operon-info row, in order -----------------------
    output_columns = [
        "feature_id", "organism_name", "best_consensus_product_descriptor",
        "product_descriptor_source", "product_descriptor_source_id", "operon_id",
        "c3_score", "c3_descriptor_conflict", "c3_operon_significance",
        "c3_operon_ambiguity",
        "c3_signal_breakdown", "c3_formula",
        "c3_gene_class", "c3_operon_id", "c3_n_pairs_used", "c3_n_supported",
        "c3_n_neutral", "c3_n_unlinked", "c3_lowerbound",
        # both variants; c3_score is the --c3-mode selection
        "c3_score_adjacency", "c3_lowerbound_adjacency",
        "c3_score_hybrid", "c3_lowerbound_hybrid",
    ]
    n = 0
    score_sum = 0.0
    cls_counts = {}
    with open(args.operon_info, newline="") as fh_in, \
            open(args.output, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")
        writer = csv.DictWriter(fh_out, fieldnames=output_columns,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            raw_oid = (row.get("operon_id") or "").strip()
            info = per_gene.get(fid, dict(
                c3=_NEUTRAL, c3_lowerbound=_NEUTRAL, gene_class="singleton",
                c3_adjacency=_NEUTRAL, c3_lb_adjacency=_NEUTRAL,
                c3_hybrid=_NEUTRAL, c3_lb_hybrid=_NEUTRAL,
                operon_id="", n_pairs=0, n_pairs_used=0, n_supported=0,
                n_neutral=0, n_nolink=0, descriptor_conflict=0.0,
                operon_significance=0.0, operon_ambiguity=0.0))
            reason, formula = _breakdown_and_formula(info, raw_oid)

            out_row = {col: row.get(col, "") for col in output_columns}
            out_row["c3_score"] = f"{info['c3']:.4f}"
            out_row["c3_descriptor_conflict"] = f"{info['descriptor_conflict']:.4f}"
            out_row["c3_operon_significance"] = f"{info['operon_significance']:.4f}"
            out_row["c3_operon_ambiguity"] = f"{info['operon_ambiguity']:.4f}"
            out_row["c3_signal_breakdown"] = reason
            out_row["c3_formula"] = formula
            out_row["c3_gene_class"] = info["gene_class"]
            out_row["c3_operon_id"] = info["operon_id"]
            out_row["c3_n_pairs_used"] = info["n_pairs_used"]
            out_row["c3_n_supported"] = info["n_supported"]
            out_row["c3_n_neutral"] = info["n_neutral"]
            out_row["c3_n_unlinked"] = info["n_nolink"]
            out_row["c3_lowerbound"] = f"{info['c3_lowerbound']:.4f}"
            out_row["c3_score_adjacency"] = f"{info['c3_adjacency']:.4f}"
            out_row["c3_lowerbound_adjacency"] = f"{info['c3_lb_adjacency']:.4f}"
            out_row["c3_score_hybrid"] = f"{info['c3_hybrid']:.4f}"
            out_row["c3_lowerbound_hybrid"] = f"{info['c3_lb_hybrid']:.4f}"
            writer.writerow(out_row)

            n += 1
            score_sum += info["c3"]
            cls_counts[info["gene_class"]] = cls_counts.get(info["gene_class"], 0) + 1

    print(f"\n[c3-scorer] wrote {n} genes -> {args.output}")
    if n:
        print(f"  mean C3: {score_sum / n:.4f}")
    print(f"  classes: {cls_counts}")


if __name__ == "__main__":
    main()
