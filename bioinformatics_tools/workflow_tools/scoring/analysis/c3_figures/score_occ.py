"""score_occ.py - CLI for the dynamic Operon Context Confidence (OCC) database.

Operationalises c3_occ.py's incremental reference: the OCC database stores only
ADDITIVE pan-genome co-occurrence counts, so every newly-labeled organism folds
in WITHOUT reprocessing the organisms already in it.  This CLI exposes that
lifecycle - build once, then update per new organism, then score.

Gene identity = clean_descriptor (functional name).  OCC is derived purely from
pan-genome operon co-occurrence (figs 01-05), INDEPENDENT of C1/C2/C4.  Read-only
prototype; the production scorer is untouched.

CONVENTIONS (match the c3fig_* scripts)
  --stats-dir  = .../scoring/c3-genes-comprehensive-stats  (holds _cache/genes.pkl)
  run_root     = stats_dir.parents[1] (the run dir with <org>/labeling/…)  - the
                 adjacency channel reads contigs from there.
  --db         = the persistent OCC database pickle
                 (default: <stats_dir>/_cache/occ_reference.pkl)

SUBCOMMANDS
  build   Build the database from scratch from every organism in the cache.
            score_occ.py build  --stats-dir DIR [--db PATH]
  update  Fold organism(s) into the database (create it if missing), then save.
          Default = every organism in the cache not already in the database - so
          re-running after a new organism is labeled just adds that organism.
            score_occ.py update --stats-dir DIR [--organisms A B …] [--db PATH]
  score   Load the database and write per-gene OCC for every gene in the cache.
            score_occ.py score  --stats-dir DIR [--db PATH] -o OUT.tsv
  status  Print a summary of the database (organisms, #pairs, params).
            score_occ.py status --db PATH

TYPICAL DYNAMIC WORKFLOW
  # first build (all organisms currently labeled)
  score_occ.py build  --stats-dir .../c3-genes-comprehensive-stats
  # later, when organism "Klebsiella_pneumoniae" is newly labeled into a run:
  score_occ.py update --stats-dir <that run>/scoring/c3-genes-comprehensive-stats \
                      --organisms Klebsiella_pneumoniae
  # score whenever needed
  score_occ.py score  --stats-dir .../c3-genes-comprehensive-stats -o per_gene_occ.tsv
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import c3_occ as O


def _resolve(stats_dir):
    stats_dir = Path(stats_dir).resolve()
    cache = stats_dir / "_cache" / "genes.pkl"
    if not cache.is_file():
        sys.exit("[score_occ] no cache at %s (run c3fig_00_build_cache.py first)" % cache)
    run_root = stats_dir.parents[1]
    return stats_dir, cache, run_root


def _default_db(stats_dir):
    return Path(stats_dir) / "_cache" / "occ_reference.pkl"


def cmd_build(args):
    stats_dir, cache, run_root = _resolve(args.stats_dir)
    db_path = Path(args.db) if args.db else _default_db(stats_dir)
    genes = L.load_cache(cache)
    orgs = sorted(genes["organism"].unique())
    print("[score_occ] building OCC database from %d organisms: %s"
          % (len(orgs), ", ".join(orgs)))
    ref = O.build_reference(genes, run_root)
    O.save_reference(ref, db_path)
    print("[score_occ] saved -> %s  (%d qualifying operons, %d adj / %d op pairs)"
          % (db_path, ref["n_qualifying_operons"], len(ref["rho_adj"]),
             len(ref["rho_op"])))


def cmd_update(args):
    stats_dir, cache, run_root = _resolve(args.stats_dir)
    db_path = Path(args.db) if args.db else _default_db(stats_dir)
    genes = L.load_cache(cache)
    if db_path.is_file():
        ref = O.load_reference(db_path)
        print("[score_occ] loaded database %s (%d organisms already in it)"
              % (db_path, len(ref["organisms_added"])))
    else:
        ref = O.new_reference()
        print("[score_occ] no database at %s - creating a new one" % db_path)

    before = set(ref["organisms_added"])
    O.update_reference(ref, genes, run_root, organisms=args.organisms)
    added = sorted(set(ref["organisms_added"]) - before)
    if not added:
        print("[score_occ] nothing to add (all requested organisms already present)")
    else:
        print("[score_occ] added %d organism(s): %s" % (len(added), ", ".join(added)))
    O.finalize_reference(ref)  # cheap (memoised); leaves the DB ready to score
    O.save_reference(ref, db_path)
    print("[score_occ] saved -> %s  (now %d organisms, %d adj / %d op pairs)"
          % (db_path, len(ref["organisms_added"]), len(ref["rho_adj"]),
             len(ref["rho_op"])))


def cmd_score(args):
    stats_dir, cache, run_root = _resolve(args.stats_dir)
    db_path = Path(args.db) if args.db else _default_db(stats_dir)
    genes = L.load_cache(cache)
    if db_path.is_file():
        ref = O.load_reference(db_path)
        print("[score_occ] scoring with database %s (%d organisms)"
              % (db_path, len(ref["organisms_added"])))
    else:
        print("[score_occ] no database at %s - building one on the fly" % db_path)
        ref = O.build_reference(genes, run_root)
    df = O.compute_all_genes(genes, run_root, ref=ref)
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.drop(columns=[c for c in df.columns if c not in (
        "organism", "feature_id", "clean_descriptor", "operon_id", "uninformative",
        "n_inf_context", "n_partners", "best_partner", "best_rho", "best_channel",
        "occ")]).to_csv(out, sep="\t", index=False)
    occ0 = df.attrs.get("occ0", float("nan"))
    n_inf = int((~df["uninformative"]).sum())
    print("[score_occ] wrote %d rows (%d informative) -> %s"
          % (len(df), n_inf, out))
    print("[score_occ] neutral pivot occ0 (informative median) = %.4f" % occ0)


def cmd_status(args):
    if not args.db:
        sys.exit("[score_occ] status needs --db PATH")
    db_path = Path(args.db)
    if not db_path.is_file():
        sys.exit("[score_occ] no database at %s" % db_path)
    ref = O.load_reference(db_path)
    orgs = sorted(ref["organisms_added"])
    print("OCC database: %s" % db_path)
    print("  version               : %s" % ref.get("version"))
    print("  finalized             : %s" % ref.get("finalized"))
    print("  organisms (%d)         : %s" % (len(orgs), ", ".join(orgs)))
    print("  qualifying operons    : %d" % ref["n_qualifying_operons"])
    print("  adjacency pairs (rho) : %d" % len(ref["rho_adj"]))
    print("  co-operon pairs (rho) : %d" % len(ref["rho_op"]))
    print("  descriptors present   : %d" % len(ref["present"]))
    print("  params                : %s" % ref["params"])


def main():
    ap = argparse.ArgumentParser(
        description="Dynamic Operon Context Confidence (OCC) database CLI.")
    sub = ap.add_subparsers(dest="cmd")
    sub.required = True

    b = sub.add_parser("build", help="build the database from scratch")
    b.add_argument("--stats-dir", required=True)
    b.add_argument("--db", default=None)
    b.set_defaults(func=cmd_build)

    u = sub.add_parser("update", help="fold organism(s) into the database")
    u.add_argument("--stats-dir", required=True)
    u.add_argument("--organisms", nargs="+", default=None,
                   help="organisms to add (default: all not yet in the database)")
    u.add_argument("--db", default=None)
    u.set_defaults(func=cmd_update)

    s = sub.add_parser("score", help="write per-gene OCC to a TSV")
    s.add_argument("--stats-dir", required=True)
    s.add_argument("--db", default=None)
    s.add_argument("-o", "--out", required=True)
    s.set_defaults(func=cmd_score)

    st = sub.add_parser("status", help="print a summary of the database")
    st.add_argument("--db", required=True)
    st.set_defaults(func=cmd_status)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
