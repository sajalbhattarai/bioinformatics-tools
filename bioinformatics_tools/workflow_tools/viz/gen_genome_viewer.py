#!/usr/bin/env python3
"""gen_genome_viewer.py -- builds a self-contained interactive HTML genome viewer
from FINAL_ANNOTATION_WITH_CONFIDENCE.tsv.

Gene mode colours genes by confidence tier with review flags; operon mode colours
each operon from a categorical cycle and highlights its members. The per-gene
data (plus evidence trail and fingerprints) is embedded as JSON and gzip/base64
in TEMPLATE, so the page works offline.

Usage: gen_genome_viewer.py <FINAL.tsv> <out.html>
"""
import base64
import csv
import gzip
import json
import re
import sys
from pathlib import Path

csv.field_size_limit(10 ** 8)

TIERS = ["highest", "high", "medium", "fair", "low"]          # index 0..4; -1 = non-coding
TIER_IDX = {t: i for i, t in enumerate(TIERS)}

# Evidence trail entries (display name, group, candidate columns) read from the
# consolidated per-tool matrix; decision databases first, then domain
# signatures, then specialised callers.
EVIDENCE = [
    ("RAST", "decision", ["RAST_description"]),
    ("COG", "decision", ["COG_description"]),
    ("Pfam", "decision", ["PFAM_description"]),
    ("KEGG", "decision", ["KEGG_description"]),
    ("eggNOG", "decision", ["EGGNOG_description"]),
    ("UniProt", "decision", ["UNIPROT_description"]),
    ("PGAP", "decision", ["PGAP_description"]),
    ("TIGRFAM", "decision", ["TIGRFAM_description"]),
    ("NCBIfam", "decision", ["INTERPRO_NCBIFAM_description"]),
    ("InterPro", "domain", ["INTERPRO_description"]),
    ("Gene3D", "domain", ["INTERPRO_GENE3D_description"]),
    ("SUPERFAMILY", "domain", ["INTERPRO_SUPERFAMILY_description"]),
    ("PANTHER", "domain", ["INTERPRO_PANTHER_description"]),
    ("CDD", "domain", ["INTERPRO_CDD_description"]),
    ("SMART", "domain", ["INTERPRO_SMART_description"]),
    ("PRINTS", "domain", ["INTERPRO_PRINTS_description"]),
    ("PROSITE", "domain", ["INTERPRO_PROSITE_PROFILES_description", "INTERPRO_PROSITE_PATTERNS_description"]),
    ("HAMAP", "domain", ["INTERPRO_HAMAP_description"]),
    ("FunFam", "domain", ["INTERPRO_FUNFAM_description"]),
    ("PIRSF", "domain", ["INTERPRO_PIRSF_description"]),
    ("SFLD", "domain", ["INTERPRO_SFLD_description"]),
    ("TCDB", "special", ["TCDB_family_description"]),
    ("MEROPS", "special", ["MEROPS_description"]),
    ("dbCAN", "special", ["DBCAN_description"]),
    ("RAST subsystem", "special", ["RASTTK_subsystem_description"]),
]
EV_NAMES = [e[0] for e in EVIDENCE]
EV_GROUPS = [e[1] for e in EVIDENCE]
_UNINF = re.compile(r"^(hypothetical|uncharacter|unknown|putative uncharacter|domain of unknown"
                    r"|conserved (hypothetical|protein)|duf\d)", re.I)
# alignment-based tools carry a % identity worth showing on the evidence trail
IDENTITY_COL = {"UniProt": "UNIPROT_percent_identity", "COG": "COG_identity",
                "MEROPS": "MEROPS_percent_identity", "TCDB": "TCDB_percent_identity"}


def fmt_pid(v):
    """Formats a fraction or percentage as "NN% id", empty when not numeric."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return ""
    if x <= 1.0:
        x *= 100
    return f"{round(x)}% id"


# Gene-call columns carry the caller's prefix (RAST_start, PRODIGAL_start), so
# they are matched by the part after the prefix.
GENE_CALLERS = ("rast", "rasttk", "prodigal")


def col(row, name):
    """Returns a row value by column name, ignoring "Column-X: " prefixes (regex) and the gene-caller prefix."""
    want = name.strip().lower()
    head, _, rest = want.partition("_")
    alt = rest if head in GENE_CALLERS and rest else ""
    for k in row:
        bare = re.sub(r"^Column-[A-Z]+:\s*", "", k or "").strip().lower()
        if bare == want:
            return row[k] or ""
        if alt:
            k_head, _, k_rest = bare.partition("_")
            if k_rest == alt and k_head in GENE_CALLERS:
                return row[k] or ""
    return ""


def clean_ev(d):
    """Strips source and gnl| prefixes from an evidence description by regex and caps it at 120 characters."""
    d = re.sub(r"^[A-Za-z][\w /()]*?:\s*", "", str(d)).strip()      # drop "JCVI:"/"KEGG:" prefixes
    d = re.sub(r"^gnl\|[^|]*\|[^|]*\|\S*\s*", "", d).strip()         # drop gnl|DB|acc| prefixes
    return d[:120]                                                   # keep [EC:...] tags — the EC evidence


def contig_of(gid):
    """Returns the contig part of a "<contig>_<start><strand><len>" gene_id (regex)."""
    return re.sub(r"_[0-9]+[+-][0-9]+$", "", gid)


def num(v, nd=2):
    """Returns v rounded to nd decimals, or None when not numeric."""
    try:
        return round(float(v), nd)
    except (TypeError, ValueError):
        return None


# Locations of the consolidated matrix relative to the FINAL table's directory,
# before and after reorganize_outputs.py; --consolidated overrides them.
_CONS_NAME = "consolidated-merged-all-columns.tsv"
_CONS_CANDIDATES = (
    Path("per-tool-phased-output") / "consolidation" / _CONS_NAME,  # reorganized
    Path("consolidation") / _CONS_NAME,                             # mid-run, FINAL at top
    Path("..") / "consolidation" / _CONS_NAME,                      # mid-run, FINAL in scoring/
)


def find_consolidated(final_path, explicit=None):
    """Returns the consolidated matrix path from the explicit path or known locations, or None.

    A missing explicit path exits with an error.
    """
    if explicit:
        p = Path(explicit)
        if not p.is_file():
            raise SystemExit(f"--consolidated not found: {p}")
        return p
    base = Path(final_path).resolve().parent
    for rel in _CONS_CANDIDATES:
        p = (base / rel).resolve()
        if p.is_file():
            return p
    return None


def load_fingerprints(final_path):
    """Reads labeled-genes-fingerprint-full.tsv with csv into {feature_id: {h, l, f, raw}}, or {} when absent.

    Parses "pattern hash: X || label: Y || fingerprint: F1: v | F2: v ..."; empty fields are dropped.
    """
    base = Path(final_path).resolve().parent
    for rel in ("../fingerprint/labeled-genes-fingerprint-full.tsv", "fingerprint/labeled-genes-fingerprint-full.tsv"):
        p = (base / rel).resolve()
        if p.is_file():
            break
    else:
        return {}
    out = {}
    for r in csv.DictReader(open(p, newline=""), delimiter="\t"):
        cell = r.get("fingerprint", "") or ""
        parts = {k.strip(): v.strip() for k, _, v in (seg.partition(":") for seg in cell.split(" || "))}
        fields = []
        for seg in parts.get("fingerprint", "").split(" | "):
            k, _, v = seg.partition(":")
            if k.strip() and v.strip():
                fields.append([k.strip(), v.strip()])
        out[r.get("feature_id", "")] = {"h": parts.get("pattern hash", ""), "l": parts.get("label", ""), "f": fields,
                                        "raw": parts.get("fingerprint", "")}
    print(f"fingerprints from {p}", file=sys.stderr)
    return out


def main():
    """Reads the FINAL table and evidence with csv, packs the per-gene data and writes the HTML viewer."""
    argv = [a for a in sys.argv[1:] if a != "--artifact"]
    explicit_cons = None
    if "--consolidated" in argv:
        i = argv.index("--consolidated")
        explicit_cons = argv[i + 1]
        del argv[i:i + 2]
    if len(argv) < 2:
        raise SystemExit(
            "usage: gen_genome_viewer.py <FINAL.tsv> <out.html> "
            "[--consolidated <consolidated-merged-all-columns.tsv>] [--artifact]"
        )
    final, out = argv[0], argv[1]
    rows = list(csv.DictReader(open(final, newline=""), delimiter="\t"))
    if not rows:
        raise SystemExit(f"no rows in {final}")
    organism = col(rows[0], "organism_name")

    # consolidated per-tool matrix (for the evidence trail); join by feature_id
    cons_path = find_consolidated(final, explicit_cons)
    cons = {}
    if cons_path is not None:
        print(f"evidence trail from {cons_path}", file=sys.stderr)
        for r in csv.DictReader(open(cons_path, newline=""), delimiter="\t"):
            cons[r.get("feature_id", "")] = r
    else:
        # Not fatal: only the evidence trail is empty.
        print(f"WARNING: no {_CONS_NAME} near {final}; "
              "evidence trail will be empty", file=sys.stderr)

    def evidence_of(fid):
        """Returns the evidence trail [[tool index, description, identity, informative]] for a gene."""
        cr = cons.get(fid)
        if not cr:
            return []
        # Every tool with a result is listed, informative or not.
        trail, seen = [], set()
        for ti, (name, grp, colnames) in enumerate(EVIDENCE):
            val = ""
            for c in colnames:
                v = clean_ev(cr.get(c, "") or "")
                if v and len(v) > 2:
                    val = v
                    break
            if not val:
                continue
            inf = 0 if _UNINF.match(val) else 1
            key = re.sub(r"[^a-z0-9]", "", val.lower())
            if grp == "domain":                    # domain sigs: skip uninformative + duplicates (noise)
                if not inf or key in seen:
                    continue
            seen.add(key)
            mc = IDENTITY_COL.get(name)
            m = fmt_pid(cr.get(mc, "")) if mc else ""
            trail.append([ti, val, m, inf])
        return trail

    def uniprot_hit(fid):
        """Returns UniProt's best hit (identity, entry, description, informative) for a gene, or None."""
        cr = cons.get(fid)
        if not cr:
            return None
        pid = cr.get("UNIPROT_percent_identity", "").strip()
        desc = (cr.get("UNIPROT_description", "") or "").strip()
        en = (cr.get("UNIPROT_entry_name", "") or "").strip()
        if not pid and not desc:
            return None
        try:
            pidv = round(float(pid))
        except (TypeError, ValueError):
            pidv = None
        return {"pid": pidv, "en": en, "desc": desc[:80],
                "inf": bool(desc) and not _UNINF.match(clean_ev(desc))}

    # contig lengths (max end), largest first
    clen = {}
    for r in rows:
        gid = col(r, "gene_id")
        if not gid:
            continue
        try:
            e = int(col(r, "RAST_end"))
        except ValueError:
            continue
        c = contig_of(gid)
        clen[c] = max(clen.get(c, 0), e)
    order = sorted(clen, key=lambda c: -clen[c])
    cidx = {c: i for i, c in enumerate(order)}
    contigs = [{"name": c, "len": clen[c]} for c in order]

    fps = load_fingerprints(final)
    genes = []
    # Longer per-gene report text in XF order, gzipped into the page and unpacked
    # only when a report is downloaded.
    XF = ["gene_id", "FEATURE_TYPE", "ENVELOPE", "ENVELOPE_reason", "best_consensus_product_descriptor_source_hierarchy_order",
          "best_consensus_product_descriptor_source_audit_trail", "specialized_database_hits", "localization_and_topology_hits",
          "C1_reasoning", "C2_score_formula", "C2_score_reasoning", "C3_score_formula", "C3_score_reasoning",
          "C4_score_formula", "C4_score_reasoning", "PRELIMINARY_confidence_C1_C4_formula", "PRELIMINARY_confidence_C1_C4_reasoning",
          "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_formula", "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_reasoning",
          "DOES_OPERON_CONTEXT_IMPROVE_CONFIDENCE?", "NEEDS_REVIEW_REASON"]
    extras = []
    for r in rows:
        gid = col(r, "gene_id")
        if not gid:
            continue
        try:
            s, e = int(col(r, "RAST_start")), int(col(r, "RAST_end"))
        except ValueError:
            continue
        tier = col(r, "CONFIDENCE_TIER").strip()
        ti = TIER_IDX.get(tier, -1)
        opid = col(r, "UniOP_OPERON_id").strip()
        in_op = col(r, "IS_IN_OPERON?").strip().lower() == "yes" and opid.startswith("operon_")
        rv = col(r, "NEEDS_REVIEW?").strip().lower() == "yes"
        c4v = num(col(r, "C4_score_EC_conflict"))
        genes.append({
            "s": s, "e": e,
            "st": -1 if col(r, "RAST_strand").strip() == "-" else 1,
            "ci": cidx[contig_of(gid)],
            "ti": ti,
            "nm": col(r, "best_consensus_product_descriptor") or col(r, "BEST_PRODUCT_DESCRIPTOR(copied_here_for_convenience)"),
            "op": opid if in_op else None,
            "pr": num(col(r, "UniOP_operon_probability")),
            "c1": num(col(r, "C1_score_database_coverage")),
            "c2": num(col(r, "C2_score_pairwise_genes_UniOP_probability")),
            "c3": num(col(r, "C3_score_operon_context")),
            "c3h": num(col(r, "C3_score_operon_context_hybrid")),
            "c4": c4v,
            "pre": num(col(r, "PRELIMINARY_confidence_C1_C4")),
            "fin": num(col(r, "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT")),
            "finh": num(col(r, "ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_hybrid")),
            "tih": TIER_IDX.get(col(r, "CONFIDENCE_TIER_hybrid").strip(), -1),
            "imp": 1 if col(r, "DOES_OPERON_CONTEXT_IMPROVE_CONFIDENCE?").strip().lower() == "yes" else 0,
            "ecs": col(r, "EC_EVIDENCE_STATUS"),
            "c4r": col(r, "C4_score_reasoning")[:200] if (c4v is not None and c4v < 1) else "",
            "rv": 1 if rv else 0,
            "rr": col(r, "NEEDS_REVIEW_REASON") if rv else "",
            "src": col(r, "best_consensus_product_descriptor_source"),
            "up": uniprot_hit(col(r, "RAST_feature_id")),
            "ev": evidence_of(col(r, "RAST_feature_id")),
            "fid": col(r, "RAST_feature_id"),
            "fp": ({k: v for k, v in fps[col(r, "RAST_feature_id")].items() if k != "raw"} if col(r, "RAST_feature_id") in fps else None),
        })
        extras.append([col(r, k) for k in XF] + [(fps.get(col(r, "RAST_feature_id")) or {}).get("raw", "")])

    n_op = len({g["op"] for g in genes if g["op"]})
    n_flag = sum(g["rv"] for g in genes)
    # Display label is the genome identifier as given.
    short = organism

    data = {
        "organism": organism, "short": short,
        "contigs": contigs, "genes": genes,
        "nOperons": n_op, "nFlag": n_flag,
        "totLen": sum(clen.values()),
        "evNames": EV_NAMES, "evGroups": EV_GROUPS,
    }
    packed = base64.b64encode(gzip.compress(json.dumps({"fields": XF + ["fingerprint_line"], "rows": extras},
                                                        separators=(",", ":")).encode("utf-8"), 9)).decode("ascii")
    html = TEMPLATE.replace("/*__DATA__*/", json.dumps(data, separators=(",", ":"))).replace("__XTRA__", packed)

    if "--artifact" in sys.argv:
        # Artifact mode emits only title, style and body content; the host supplies the document shell.
        style = re.search(r"<style>.*?</style>", html, re.S).group(0)
        inner = re.search(r"<body>(.*)</body>", html, re.S).group(1)
        html = (f"<title>{short} genome — MARGIE confidence viewer</title>\n"
                f"{style}\n{inner}")

    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    print(f"wrote {out}  ({short}: {len(genes)} genes, {len(order)} replicons, "
          f"{n_op} operons, {n_flag} flagged)  {out.stat().st_size//1024} KB")


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MARGIE interactive genome viewer</title>
<style>
  :root{
    /* pure white page, black text, no shading (single theme, by request) */
    --surface:#ffffff; --ink:#000000; --muted:#000000; --line:#cfcfcf; --panel:#ffffff;
    --btn:#ffffff; --btn-ink:#000000; --shadow:none;
    /* confidence tiers, the same five colours the app's tables use */
    --t0:#1F77FF; --t1:#00B84D; --t2:#FFCC00; --t3:#FF8C00; --t4:#EE2233; --tn:#bdbdbd;
    --operon:#1667e0; --nonop:#c79a5c; --flag:#666666;
  }
  *{box-sizing:border-box}
  [hidden]{display:none !important}
  html,body{margin:0;background:#ffffff;color:#000000;
    /* Serif throughout, per request. Times New Roman first, with a
       platform-serif fallback chain so Linux and older browsers do not
       silently drop to a sans default. */
    font-family:"Times New Roman",Times,Georgia,"Liberation Serif","DejaVu Serif",serif;
    -webkit-font-smoothing:antialiased}
  .wrap{display:flex;flex-wrap:wrap;gap:16px;padding:18px;max-width:1280px;margin:0 auto}
  .left{flex:1 1 560px;min-width:340px}
  .right{flex:1 1 360px;min-width:300px;max-width:440px}
  h1{font-size:18px;margin:0 0 2px;overflow-wrap:anywhere;word-break:break-word}
  h1 em{font-style:italic}
  .sub{color:var(--muted);font-size:13px;margin-bottom:12px}
  .modebar{display:flex;gap:0;margin:6px 0 4px;border:1px solid var(--line);border-radius:3px;overflow:hidden;width:max-content;box-shadow:var(--shadow)}
  .modebar button{font-family:inherit;font-size:14px;padding:7px 20px;border:0;background:var(--btn);color:var(--btn-ink);cursor:pointer;transition:background .12s}
  .modebar button.on{background:var(--ink);color:var(--surface)}
  .opts{font-size:12px;color:var(--muted);margin:6px 0 2px;display:flex;gap:14px;align-items:center}
  .opts label{cursor:pointer;user-select:none}
  .plate{background:#ffffff;border:1px solid var(--line);border-radius:4px;padding:6px;box-shadow:var(--shadow)}
  /* The circle reads better small: it leaves room for the scorecard beside it. */
  #circPlate{max-width:560px}
  .viewbar{display:flex;flex-wrap:wrap;align-items:center;gap:10px 16px;margin:8px 0 6px}
  .modebar.small button{font-size:13px;padding:5px 14px}
  .linctl{display:flex;flex-wrap:wrap;align-items:center;gap:10px 14px;font-size:12.5px;color:var(--muted)}
  .linctl select{font-family:inherit;font-size:12.5px;padding:2px 4px;border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--ink)}
  .linpos{font-variant-numeric:tabular-nums}
  .linplate{max-width:none}
  .linscroll{position:relative;overflow-x:auto;overflow-y:hidden}
  .linspacer{height:1px}
  #lin{position:sticky;left:0;width:100%;height:210px;display:block}
  .strandkey{font-size:12px;margin-top:6px}
  .zsel{fill:rgba(31,119,255,.14);stroke:#1F77FF;stroke-width:1.5;stroke-dasharray:4 3}
  .zbtn.on{background:var(--ink);color:#fff;border-color:var(--ink)}
  svg{width:100%;height:auto;display:block;touch-action:none}
  .gene{cursor:pointer}
  .sel{stroke:#111;stroke-width:2;paint-order:stroke}
  .dim{opacity:.16}
  .hi{stroke:#111;stroke-width:1.1}
  .plate{position:relative}
  .zoombar{position:absolute;top:6px;right:6px;z-index:2;display:flex;gap:4px;align-items:center}
  .zbtn{font-family:inherit;font-size:14px;line-height:1;width:26px;height:26px;
        border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--ink);
        cursor:pointer;padding:0;display:flex;align-items:center;justify-content:center}
  .zbtn:hover{background:#f4f4f4}
  .zreset{font-family:inherit;font-size:11.5px;padding:4px 8px;border:1px solid var(--line);
          border-radius:6px;background:#fff;color:var(--ink);cursor:pointer}
  .zreset:disabled{opacity:.35;cursor:default}
  .zhint{font-size:11.5px;color:var(--muted);text-align:right;margin-top:-2px}
  .legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:10px;font-size:12.5px}
  .legend span{display:inline-flex;align-items:center;gap:6px;color:var(--ink)}
  .sw{width:14px;height:14px;border-radius:3px;display:inline-block;border:1px solid rgba(0,0,0,.12)}
  .panel{background:var(--panel);border:1px solid var(--line);border-radius:4px;padding:14px 16px;min-height:280px}
  .panel .empty{color:var(--muted);font-size:13px;line-height:1.5}
  .ptitle{font-size:17px;margin:0 0 2px;line-height:1.25}
  .ptag{display:inline-block;font-size:12px;padding:1px 7px;border-radius:2px;color:#fff;margin-bottom:8px}
  .cardhead{display:flex;justify-content:space-between;align-items:flex-start;gap:12px}
  .cardhead-l{min-width:0}
  .dlbtn{flex:0 0 auto;font-family:inherit;font-size:11.5px;padding:5px 11px;border:1px solid var(--line);
    background:var(--btn);color:var(--btn-ink);border-radius:3px;cursor:pointer;white-space:nowrap;transition:background .12s}
  .dlbtn:hover{background:var(--ink);color:var(--surface);border-color:var(--ink)}
  /* A boxed table, not free text: every figure sits in its own cell. */
  .kv{display:grid;grid-template-columns:auto 1fr;font-size:13px;margin:10px 0;
      border:1px solid var(--line);border-radius:3px;overflow:hidden}
  .kv b,.kv span{padding:5px 10px;border-top:1px solid var(--line)}
  .kv b:first-child,.kv b:first-child + span{border-top:0}
  .kv b{color:var(--muted);font-weight:normal;border-right:1px solid var(--line);background:#fafafa}
  .kv span{font-variant-numeric:tabular-nums}
  .two{display:inline-grid;grid-auto-flow:column;gap:0 10px}
  .two i{font-style:normal;color:var(--muted)}
  .bar{height:9px;border-radius:0;background:#ffffff;border:1px solid #cfcfcf;overflow:hidden;margin:2px 0}
  .bar>i{display:block;height:100%}
  .members{margin-top:10px;max-height:260px;overflow:auto;border-top:1px solid var(--line);padding-top:8px}
  .mrow{display:flex;align-items:center;gap:8px;font-size:12px;padding:3px 6px;border-radius:6px;cursor:pointer}
  .mrow:hover{background:rgba(127,127,127,.16)}
  .dot{width:10px;height:10px;border-radius:50%;flex:0 0 auto}
  .mrow .mnm{flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .mrow .mfin{color:var(--muted);font-variant-numeric:tabular-nums}
  .mrow .mnum{width:16px;flex:0 0 auto;text-align:right;color:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}
  .omap{width:100%;height:auto;display:block;margin:8px 0 2px}
  .omap polygon{cursor:pointer}
  .flagtag{color:var(--flag);font-size:11.5px;margin-top:6px}
  .trail{margin-top:12px;border-top:1px solid var(--line);padding-top:9px;max-height:340px;overflow:auto}
  .fphash{margin-top:4px;font-size:11px}
  /* A tall window (full screen, or a big monitor): the map grows with the height
     and the page uses the whole width instead of leaving margins. */
  @media (min-height:820px){
    .wrap{max-width:none}
    #circPlate{max-width:min(100%,calc(100vh - 180px))}
    .right{max-width:520px}
  }
  /* contigs: the list view, and the backbone's alternating shades */
  .ctgplate{padding:0}
  .ctghead{display:flex;justify-content:space-between;gap:12px;padding:8px 12px;border-bottom:1px solid var(--line);font-size:12.5px}
  .ctglist{max-height:min(640px,calc(100vh - 260px));overflow:auto}
  .crow{display:grid;grid-template-columns:minmax(0,1fr) 90px 60px 140px;gap:10px;align-items:center;padding:5px 12px;border-bottom:1px solid #eee;font-size:12.5px;cursor:pointer}
  .crow:hover,.crow.on{background:#f3f6ff}
  .crow .cnm{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .crow .cnum{text-align:right;font-variant-numeric:tabular-nums}
  .tbar,.kv .tbar{display:flex;height:8px;border:1px solid #cfcfcf;overflow:hidden;padding:0;width:100%}
  .tbar i{display:block;height:100%}
  .bb{cursor:pointer}
  .mrow .mnum.kb{width:auto;white-space:nowrap}
  .trail h4{font-size:12px;margin:0 0 2px;color:var(--muted);font-weight:normal;letter-spacing:.5px;text-transform:uppercase}
  .trail .cnt{color:var(--muted);font-size:12px;margin-bottom:4px}
  .grp{font-size:11.5px;color:var(--muted);margin:9px 0 2px;font-style:italic}
  .erow{display:grid;grid-template-columns:88px 1fr;gap:8px;font-size:11.5px;padding:1.5px 0;line-height:1.3}
  .erow .etool{color:var(--muted);white-space:nowrap}
  .erow.win .etool{color:var(--operon);font-weight:bold}
  .erow.win .edesc{font-weight:bold}
  .chosen{display:inline-block;font-size:9.5px;color:#fff;background:var(--operon);border-radius:2px;padding:1px 7px;margin-left:7px;vertical-align:1px}
  .emetric{display:inline-block;font-size:10.5px;color:var(--operon);font-variant-numeric:tabular-nums}
  .erow.uninf{opacity:.5}
  .ecflag{font-size:11.5px;color:#c0143c;margin:6px 0 2px;line-height:1.35}
  .uphit{margin-top:8px;font-size:11.5px;color:var(--muted);line-height:1.4}
  .uphit b{color:var(--ink);font-weight:normal}
  .uphit .uninf{color:#c0143c}
  .oplink{color:var(--operon);cursor:pointer;text-decoration:underline;text-underline-offset:2px}
  #tip{position:fixed;pointer-events:none;background:#ffffff;color:#000000;font-size:11.5px;
    padding:6px 9px;border:1px solid #000000;max-width:280px;opacity:0;transition:opacity .08s;z-index:9;line-height:1.35}
  #tip b{color:#000000;font-weight:bold}
  .foot{font-size:11.5px;color:#000000;margin-top:12px}
</style>
</head>
<body>
<div class="wrap">
  <div class="left">
    <h1><span id="org"></span></h1>
    <div class="sub" id="sub"></div>
    <div class="modebar">
      <button id="mGene" class="on">Gene mode</button>
      <button id="mOperon">Operon mode</button>
      <button id="mReview">Review flags</button>
    </div>
    <div class="opts">
      <label><input type="checkbox" id="showFlags"> show review flags (grey)</label>
      <span id="hint"></span>
    </div>
    <div class="viewbar">
      <div class="modebar small">
        <button id="vCirc" class="on" type="button">Circular</button>
        <button id="vLin" type="button">Linear</button>
        <button id="vCtg" type="button" hidden>Contigs</button>
      </div>
      <span id="linctl" class="linctl" hidden>
        <label><span id="linWhat">replicon</span> <select id="linContig"></select></label>
        <label>genes shown <select id="linN"><option>10</option><option>20</option><option>40</option></select></label>
        <label>strands <select id="linStr"><option value="one">one line</option><option value="two">two strands</option></select></label>
        <span id="linPos" class="linpos"></span>
      </span>
    </div>
    <div class="plate" id="circPlate">
      <div class="zoombar">
        <button id="zarea" class="zbtn" type="button" title="Drag a box to zoom into it (or hold Shift and drag)" aria-label="Zoom into an area"><svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true"><rect x="2.5" y="2.5" width="11" height="11" fill="none" stroke="currentColor" stroke-width="1.4" stroke-dasharray="3 2"/></svg></button>
        <button id="zout" class="zbtn" type="button" title="Zoom out" aria-label="Zoom out">−</button>
        <button id="zin" class="zbtn" type="button" title="Zoom in" aria-label="Zoom in">+</button>
        <button id="zreset" class="zreset" type="button" title="Reset zoom (or double-click the map)" disabled>Reset</button>
      </div>
      <svg id="map" viewBox="0 0 780 780" aria-label="circular genome map: scroll to zoom, drag to pan, drag a box to zoom into it"></svg>
      <div class="zhint">+ / − or scroll to zoom | drag to pan | Shift-drag (or the dotted button) to box a region | double-click to reset</div>
    </div>
    <div class="plate linplate" id="linPlate" hidden>
      <div class="linscroll" id="linScroll"><div class="linspacer" id="linSpacer"></div><svg id="lin" aria-label="linear gene map"></svg></div>
      <div class="zhint">scroll sideways through the replicon | click a gene for its scorecard</div>
    </div>
    <div class="plate ctgplate" id="ctgPlate" hidden>
      <div class="ctghead"><span id="ctgCount"></span><span>longest first | click a contig for its genes</span></div>
      <div class="ctglist" id="ctgList"></div>
    </div>
    <div class="legend" id="legend"></div>
    <div class="foot">Every value is read verbatim from FINAL_ANNOTATION_WITH_CONFIDENCE.tsv. Click a gene or operon for details.</div>
  </div>
  <div class="right">
    <div class="panel" id="panel"><div class="empty">Hover to preview, click to pin.<br><br>In <b>Gene mode</b>, each arc is a gene coloured by its confidence tier. In <b>Operon mode</b>, each operon takes its own colour and grey arcs are non-operonic; selecting an operon shows its member genes and its context scores. <b>Linear</b> lays the same genes out along the replicon, ten at a time.</div></div>
  </div>
</div>
<div id="tip"></div>
<script type="application/octet-stream" id="xtra">__XTRA__</script>
<script>
const D = /*__DATA__*/;
const TIER_NAMES=["highest","high","medium","fair","low"];
// Ordered tiers, blue (highest) through green and amber to red (low): the
// same five colours the app's tables and rings use, so one gene reads the
// same everywhere. Bright and far apart at a glance, which a single-hue ramp
// was not on a dense map. Keep in sync with TIER in make_circular_genome.py.
const TIER_COL=["#1F77FF","#00B84D","#FFCC00","#FF8C00","#EE2233"];
// FLAG is a RESERVED status colour -- never reused as a tier step, so an
// alarm can never be confused with a ranking. Operon/non-operon are a separate
// categorical pair used only in operon mode, where the tier ramp is not shown.
// Operon identity is CATEGORICAL -- unlike tiers, which are ordered -- so
// distinct hues are the right rule here, exactly where a ramp would be wrong.
// Cycled by a hash of the operon id, so an operon keeps its colour across
// renders and modes instead of depending on iteration order.
//
// Bright by request, and still spread along the blue-yellow axis that survives
// colour-vision deficiency, with the hues separated in lightness too. Six is
// the ceiling; the cycle repeats past that, which is fine because the colour
// delineates neighbouring operons rather than naming them (hover names the
// operon, click opens its card).
const OPERON_CYCLE=["#1F77FF","#FF8C00","#00B84D","#B65CFF","#00C2D1","#EE2233"];
// Retained for the legend swatch and the operon card tag.
const OPERON=OPERON_CYCLE[0];
const NONCODE="#d5d5d5", NONOP="#6b6b6b", FLAG="#b32b1e";
// Stable per-operon colour: hash the id so a given operon keeps its colour
// across renders and modes, instead of depending on iteration order.
function opColor(id){
  let h=0; for(let i=0;i<id.length;i++) h=(h*31+id.charCodeAt(i))|0;
  return OPERON_CYCLE[Math.abs(h)%OPERON_CYCLE.length];
}
const CX=390, CY=396, START=90;
// A genome in more than one piece shows no contig names, ticks or positions:
// with a draft's hundreds of contigs they piled up on each other, and even a
// few crowd the ring. The pieces sit end to end with a small gap between.
// Past MANY contigs the gaps share at most 60 degrees and each contig gets its
// length's share of the rest (a fixed 4-degree gap and floor each wrapped the
// ring round itself).
const MANY=24;
const DRAFT=D.contigs.length>MANY, MULTI=D.contigs.length>1;
const GAP=DRAFT?Math.min(4,60/D.contigs.length):4, MINSPAN=DRAFT?0:4;
const R={bbO:352,bbI:343, fO:335,fI:306, rO:301,rI:272, tick:262};

// ---- layout: each contig an arc, width ~ length, floor for tiny ones ----
const order=D.contigs.map((c,i)=>i);
const total=D.contigs.reduce((a,c)=>a+c.len,0);
const spanTotal=360-order.length*GAP;
let floored=D.contigs.filter(c=>spanTotal*c.len/total<MINSPAN);
let fixed=floored.length*MINSPAN;
let bigLen=D.contigs.filter(c=>spanTotal*c.len/total>=MINSPAN).reduce((a,c)=>a+c.len,0)||1;
let rest=Math.max(spanTotal-fixed,MINSPAN);
const lay=[]; let cur=START;
D.contigs.forEach(c=>{
  cur-=GAP;
  const span=(spanTotal*c.len/total<MINSPAN)?MINSPAN:rest*c.len/bigLen;
  lay.push({start:cur,span}); cur-=span;
});
const ang=(ci,pos)=>lay[ci].start - lay[ci].span*(pos/D.contigs[ci].len);
const pol=(r,deg)=>{const a=deg*Math.PI/180;return [CX+r*Math.cos(a), CY - r*Math.sin(a)];};
const P=(r,d)=>{const p=pol(r,d);return p[0].toFixed(1)+","+p[1].toFixed(1);};

// ---- build operon index ----
const operons={};
D.genes.forEach((g,i)=>{ if(g.op){ (operons[g.op]=operons[g.op]||[]).push(i); } });

// ---- SVG scaffold: backbone + ticks ----
const svg=document.getElementById("map");
const NS="http://www.w3.org/2000/svg";
function el(tag,attrs){const e=document.createElementNS(NS,tag);for(const k in attrs)e.setAttribute(k,attrs[k]);return e;}
function annulusStrip(ci,ri,ro,fill){
  // backbone drawn as a fan of small quads across the contig arc (no arc-flag math)
  const s=lay[ci].start, sp=lay[ci].span, seg=Math.max(2,Math.ceil(sp/2));
  for(let k=0;k<seg;k++){
    const a0=s-sp*k/seg, a1=s-sp*(k+1)/seg;
    const q=el("polygon",{points:`${P(ro,a0)} ${P(ro,a1)} ${P(ri,a1)} ${P(ri,a0)}`,fill,stroke:"none",class:"bb"});
    q.dataset.ci=ci; svg.appendChild(q);
  }
}
// A genome in pieces: the backbone alternates black and grey from one contig
// to the next, so where each begins and ends shows without labels.
const CONTIG_SHADE=["#000000","#8c8c8c"];
order.forEach((ci,k)=>{
  annulusStrip(ci,R.bbI,R.bbO,MULTI?CONTIG_SHADE[k%2]:"#000000");
  const L=D.contigs[ci].len, step=L>3e6?1e6:5e5;
  // Positions only along a single replicon (see MULTI above).
  if(!MULTI) for(let p=0;p<=L;p+=step){
    const a=ang(ci,p),[x0,y0]=pol(R.bbO,a),[x1,y1]=pol(R.bbO+8,a);
    svg.appendChild(el("line",{x1:x0,y1:y0,x2:x1,y2:y1,stroke:"#000000","stroke-width":.8}));
    if(L>=25e4 && p%step===0){const[tx,ty]=pol(R.bbO+20,a);
      const t=el("text",{x:tx,y:ty,"font-size":10,fill:"#000000","text-anchor":"middle","dominant-baseline":"middle"});
      t.textContent=(p/1e6).toFixed(1); svg.appendChild(t);}
  }
});

// ---- gene layer (event-delegated) ----
const geneLayer=el("g",{}); svg.appendChild(geneLayer);
const flagLayer=el("g",{}); svg.appendChild(flagLayer);
const nodes=[];
D.genes.forEach((g,i)=>{
  let aHi=ang(g.ci,g.s), aLo=ang(g.ci,g.e); if(aHi<aLo){const t=aHi;aHi=aLo;aLo=t;}
  if(aHi-aLo<0.05) aHi=aLo+0.05;
  const ro=g.st>0?R.fO:R.rO, ri=g.st>0?R.fI:R.rI;
  const poly=el("polygon",{points:`${P(ro,aHi)} ${P(ro,aLo)} ${P(ri,aLo)} ${P(ri,aHi)}`,class:"gene"});
  poly.dataset.i=i; if(g.op)poly.dataset.op=g.op;
  geneLayer.appendChild(poly); nodes.push(poly);
  if(g.rv){const mid=(aHi+aLo)/2;const[x0,y0]=pol(R.tick-6,mid),[x1,y1]=pol(R.tick+6,mid);
    const ln=el("line",{x1:x0,y1:y0,x2:x1,y2:y1,stroke:FLAG,"stroke-width":.9}); ln.dataset.flag=i; flagLayer.appendChild(ln);}
});
// centre label
const cInfo=el("text",{x:CX,y:CY+2,"font-size":14,fill:"#000000","text-anchor":"middle"});
cInfo.textContent=(D.totLen/1e6).toFixed(2)+" Mb  |  "+D.genes.length.toLocaleString()+" genes"; svg.appendChild(cInfo);
const cInfo2=el("text",{x:CX,y:CY+22,"font-size":12.5,fill:"#000000","text-anchor":"middle"});
svg.appendChild(cInfo2);

// ---- colouring ----
let mode="gene", selKind=null, selVal=null;
// ---- which genes to show ----
// The page that embeds the map chooses, by its address: #show=all (the
// default), a tier (highest, high, medium, fair, low), noncoding, flagged,
// operonic or nonoperonic. The rest stay on the ring, greyed like review
// mode's unflagged genes, so where the chosen ones sit still reads.
const HIDDEN="#ececea";
let SHOW="all";
function readShow(){const m=/(?:^#|&)show=([a-z]+)/.exec(location.hash);SHOW=m?m[1]:"all";}
function shown(g){
  switch(SHOW){
    case "all": return true;
    case "noncoding": return g.ti<0;
    case "flagged": return !!g.rv;
    case "operonic": return !!g.op;
    case "nonoperonic": return !g.op && g.ti>=0;
    default: { const k=TIER_NAMES.indexOf(SHOW); return k<0 || g.ti===k; }
  }
}
function geneFill(g){
  if(!shown(g)) return HIDDEN;
  if(mode==="gene") return g.ti<0?NONCODE:TIER_COL[g.ti];
  if(mode==="review") return g.rv?(g.ti<0?NONCODE:TIER_COL[g.ti]):"#ececea";
  return g.op?opColor(g.op):(g.ti<0?NONCODE:NONOP);
}
function paint(){
  D.genes.forEach((g,i)=>nodes[i].setAttribute("fill",geneFill(g)));
  document.getElementById("showFlags").checked ? flagLayer.style.display="" : flagLayer.style.display="none";
  flagLayer.querySelectorAll("line").forEach(l=>{ l.style.display=shown(D.genes[+l.dataset.flag])?"":"none"; });
  const nShown=SHOW==="all"?0:D.genes.filter(shown).length;
  cInfo2.textContent = SHOW!=="all"
    ? "showing "+nShown.toLocaleString()+" "+(SHOW==="noncoding"?"non-coding":SHOW==="operonic"?"in operons":SHOW==="nonoperonic"?"not in operons":SHOW)
    : mode==="operon"
    ? D.nOperons.toLocaleString()+" operons"
    : D.nFlag.toLocaleString()+" flagged for review";
  renderLegend(); applySelection();
  if(typeof LIN!=="undefined" && LIN.visible()) LIN.render();
}
function renderLegend(){
  const L=document.getElementById("legend"); L.innerHTML="";
  const items = mode==="operon"
    ? OPERON_CYCLE.map((c,i)=>[i===0?"operon (colour cycles)":"",c])
        .concat([["non-operonic",NONOP],["non-coding",NONCODE]])
    : mode==="review"
    ? TIER_NAMES.map((n,i)=>["flagged | "+n,TIER_COL[i]]).concat([["not flagged","#ececea"]])
    : TIER_NAMES.map((n,i)=>[n,TIER_COL[i]]).concat([["non-coding",NONCODE]]);
  items.forEach(([n,c])=>{const s=document.createElement("span");
    s.innerHTML=`<i class="sw" style="background:${c}"></i>${n}`; L.appendChild(s);});
  if(document.getElementById("showFlags").checked){const s=document.createElement("span");
    s.innerHTML=`<i class="sw" style="background:${FLAG}"></i>review flag`; L.appendChild(s);}
  {const s=document.createElement("span");
    s.innerHTML=`outer ring: + strand (5′→3′ clockwise) | inner ring: − strand (5′→3′ anticlockwise)`; L.appendChild(s);}
  if(MULTI){const s=document.createElement("span");
    s.innerHTML=`<i class="sw" style="background:linear-gradient(90deg,#000 50%,#8c8c8c 50%)"></i>contigs (backbone alternates; hover or click one)`; L.appendChild(s);}
}

// ---- contigs ----
let VIEWNOW="circ";
const contigGenes=D.contigs.map(()=>[]);
D.genes.forEach((g,i)=>contigGenes[g.ci].push(i));
const contigRank=D.contigs.map((c,i)=>i).sort((a,b)=>D.contigs[b].len-D.contigs[a].len);
function tierBar(list){
  const n=[0,0,0,0,0,0]; list.forEach(i=>{ const t=D.genes[i].ti; n[t<0?5:t]++; });
  const tot=list.length||1, col=[...TIER_COL,NONCODE];
  return `<span class="tbar" title="${TIER_NAMES.map((t,k)=>t+" "+n[k]).join(" | ")} | non-coding ${n[5]}">${n.map((v,k)=>v?`<i style="width:${v/tot*100}%;background:${col[k]}"></i>`:"").join("")}</span>`;
}
function contigCard(ci){
  const c=D.contigs[ci], list=contigGenes[ci];
  const ops=new Set(list.map(i=>D.genes[i].op).filter(Boolean)), fl=list.filter(i=>D.genes[i].rv).length;
  let h=`<div class="ptitle">${esc(c.name)}</div><span class="ptag" style="background:#555">contig ${contigRank.indexOf(ci)+1} of ${D.contigs.length}</span>`;
  h+=`<div class="kv"><b>length</b><span>${c.len.toLocaleString()} bp</span><b>genes</b><span>${list.length.toLocaleString()}</span>`
    +`<b>operons</b><span>${ops.size.toLocaleString()}</span><b>flagged</b><span>${fl.toLocaleString()}</span></div>`;
  h+=`<div class="kv"><b>confidence tiers</b><span>${tierBar(list)}</span></div>`;
  h+=`<div class="cardhead" style="margin:8px 0"><button class="dlbtn" data-ctg-map="${ci}">Show on the map</button><button class="dlbtn" data-ctg-lin="${ci}">Linear view</button></div>`;
  h+=`<div class="trail"><h4>genes, in order</h4>`;
  list.slice().sort((a,b)=>D.genes[a].s-D.genes[b].s).forEach(i=>{ const g=D.genes[i], [tn,tc]=tierBadge(g.ti);
    h+=`<div class="mrow" data-gene-open="${i}"><span class="dot" style="background:${tc}"></span><span class="mnm">${esc(g.nm||"(unnamed)")}</span><span class="mnum kb">${(g.s/1000).toFixed(1)} kb</span></div>`; });
  return h+`</div>`;
}
function showContig(ci){ selKind="contig"; selVal=ci; panel.innerHTML=contigCard(ci); applySelection();
  document.querySelectorAll(".crow").forEach(r=>r.classList.toggle("on",+r.dataset.ci===ci)); }
function renderContigList(){
  const box=document.getElementById("ctgList"); if(!box) return;
  document.getElementById("ctgCount").textContent=`${D.contigs.length.toLocaleString()} contigs | ${D.totLen.toLocaleString()} bp`;
  box.innerHTML=`<div class="crow" style="cursor:default;font-weight:bold"><span>contig</span><span class="cnum">length (bp)</span><span class="cnum">genes</span><span>tiers</span></div>`
    +contigRank.map(ci=>`<div class="crow" data-ci="${ci}"><span class="cnm" title="${esc(D.contigs[ci].name)}">${esc(D.contigs[ci].name)}</span>`
      +`<span class="cnum">${D.contigs[ci].len.toLocaleString()}</span><span class="cnum">${contigGenes[ci].length}</span>${tierBar(contigGenes[ci])}</div>`).join("");
}
document.addEventListener("click",e=>{
  const row=e.target.closest(".crow[data-ci]"); if(row){ showContig(+row.dataset.ci); return; }
  const m=e.target.closest("[data-ctg-map]"); if(m){ VIEW.set("circ"); showContig(+m.dataset.ctgMap); return; }
  const l=e.target.closest("[data-ctg-lin]"); if(l){ VIEW.set("lin"); LIN.showContig(+l.dataset.ctgLin); return; }
  const o=e.target.closest("[data-gene-open]"); if(o){ showGene(+o.dataset.geneOpen); if(LIN.visible()) LIN.goto(+o.dataset.geneOpen); }
});
svg.addEventListener("mousemove",e=>{ const t=e.target.closest(".bb"); if(!t||!MULTI) return;
  const ci=+t.dataset.ci, c=D.contigs[ci];
  tipShow(`<b>${esc(c.name)}</b><br>contig ${contigRank.indexOf(ci)+1} of ${D.contigs.length} | ${c.len.toLocaleString()} bp | ${contigGenes[ci].length} genes`,e.clientX,e.clientY); });
svg.addEventListener("click",e=>{ const t=e.target.closest(".bb"); if(t&&MULTI) showContig(+t.dataset.ci); });

// ---- selection / highlight ----
function clearFX(){nodes.forEach(n=>n.classList.remove("sel","dim","hi"));}
/*
 * Tell MARGIE (when this page is shown inside it) what is open: the view, the
 * colouring, and the gene, operon or contig selected. Chat with the genome
 * takes it as the subject of the next question.
 */
function tell(){
  if(window.parent===window) return;
  let sel=null;
  if(selKind==="gene"&&selVal!=null){ const g=D.genes[selVal]; sel={kind:"gene",id:g.fid||"",label:g.nm||"",contig:D.contigs[g.ci].name,operon:g.op||""}; }
  else if(selKind==="operon"&&operons[selVal]) sel={kind:"operon",id:selVal,label:operons[selVal].length+" genes"};
  else if(selKind==="contig"&&selVal!=null){ const c=D.contigs[selVal]; sel={kind:"contig",id:c.name,label:c.len.toLocaleString()+" bp, "+contigGenes[selVal].length+" genes"}; }
  window.parent.postMessage({margie:"selection",genome:D.short,sel,view:VIEWNOW,mode},"*");
}
function applySelection(){
  tell();
  clearFX();
  if(selKind==="operon" && operons[selVal]){
    const set=new Set(operons[selVal]);
    nodes.forEach((n,i)=>{ if(set.has(i))n.classList.add("hi"); else n.classList.add("dim"); });
  } else if(selKind==="contig" && selVal!=null){
    const set=new Set(contigGenes[selVal]);
    nodes.forEach((n,i)=>{ if(set.has(i))n.classList.add("hi"); else n.classList.add("dim"); });
  } else if(selKind==="gene" && selVal!=null){
    nodes[selVal].classList.add("sel");
  }
}
function pct(v){return v==null?"—":(+v).toFixed(2);}   // model scores are 0–1; show as decimals
function tierBadge(ti){return ti<0?["non-coding","#999"]:[TIER_NAMES[ti],TIER_COL[ti]];}
function bar(v,c){return `<div class="bar"><i style="width:${Math.round((v||0)*100)}%;background:${c}"></i></div>`;}

function geneCard(i){
  const g=D.genes[i], [tn,tc]=tierBadge(g.ti);
  const strand=g.st>0?"+":"−", cn=D.contigs[g.ci].name;
  let h=`<div class="cardhead"><div class="cardhead-l"><div class="ptitle">${esc(g.nm||"(unnamed)")}</div>`
    +`<span class="ptag" style="background:${tc}">${tn}</span></div>`
    +`<button class="dlbtn" data-gene-dl="${i}" title="Everything about this gene: identity, every score with its formula and reasoning, the fingerprint, the evidence trail, its operon's genes, and any AI interpretation">Download gene report</button></div>`;
  h+=`<div class="kv">`;
  h+=`<b>location</b><span>${cn}:${g.s.toLocaleString()}–${g.e.toLocaleString()} (${strand})</span>`;
  h+=`<b>operon</b><span>${g.op?`<a class="oplink" data-op-link="${g.op}">${g.op}</a> | P = `+pct(g.pr):"none (singleton)"}</span>`;
  h+=`</div>`;
  h+=`<div class="kv">`;
  h+=`<b>C1 database coverage</b><span>${pct(g.c1)}</span>`;
  h+=`<b>C2 operon membership</b><span>${pct(g.c2)} (${g.op?"in an operon":"singleton"})</span>`;
  h+=`<b>C3 operon context</b><span>${pct(g.c3)} (adj) | ${pct(g.c3h)} (hyb)</span>`;
  h+=`<b>C4 EC agreement</b><span>${pct(g.c4)}${g.ecs?" | "+esc(g.ecs):""}</span>`;
  h+=`</div>`;
  h+=`<div class="kv"><b>preliminary (C1 × C4)</b><span>${pct(g.pre)}</span></div>${bar(g.pre,"#888")}`;
  h+=`<div class="kv"><b>final confidence</b><span>${pct(g.fin)} (adj) | ${pct(g.finh)} (hyb)</span></div>${bar(g.fin,tc)}`;
  h+=`<div class="kv"><b>confidence tier</b><span>${tierBadge(g.ti)[0]} (adj) | ${tierBadge(g.tih)[0]} (hyb)</span></div>`;
  if(g.up)h+=`<div class="uphit"><b>UniProt best hit</b> ${g.up.pid!=null?g.up.pid+"% id":"—"}`
    +` | ${esc(g.up.en||"")} | ${esc(g.up.desc||"")}${g.up.inf?"":` <span class="uninf">(uninformative, not used)</span>`}</div>`;
  if(g.rv)h+=`<div class="flagtag">⚑ flagged for review: ${esc(g.rr||"")}</div>`;
  h+=evidenceTrail(g);
  h+=fingerprintBox(g);
  return h;
}
// The gene's fingerprint: its label, and the fields that are filled in (the hash kept small).
function fingerprintBox(g){
  const fp=g.fp;
  if(!fp) return "";
  let h=`<div class="trail"><h4>fingerprint</h4>`;
  if(fp.l) h+=`<div class="cnt">label: <b>${esc(fp.l)}</b></div>`;
  (fp.f||[]).forEach(([k,v])=>{ h+=`<div class="erow"><span class="etool">${esc(k.replace(/_/g," "))}</span><span class="edesc">${esc(v)}</span></div>`; });
  if(fp.h) h+=`<div class="cnt fphash">pattern hash ${esc(fp.h)}</div>`;
  return h+`</div>`;
}
// ---- the gene report ----------------------------------------------------
// A plain-text report of one gene: identity, confidence (every score with its
// formula and reasoning), fingerprint (one line, then database by database),
// the evidence trail, the genes of its operon (the same sections, indented),
// and, when the page is shown inside MARGIE, the AI conversation about it.
const SEP="=====x====x====";
let XTRA=null;
async function extras(){
  if(XTRA) return XTRA;
  try{
    const b64=document.getElementById("xtra").textContent.trim();
    const bin=Uint8Array.from(atob(b64),c=>c.charCodeAt(0));
    const d=await new Response(new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"))).json();
    XTRA=d.rows.map(r=>Object.fromEntries(d.fields.map((f,k)=>[f,r[k]||""])));
  }catch(err){ XTRA=D.genes.map(()=>({})); }
  return XTRA;
}
/** Fingerprint fields grouped by database: RAST_description -> RAST: description. */
function byDatabase(fields){
  const out=[], at={};
  (fields||[]).forEach(([k,v])=>{ const j=k.indexOf("_"), db=j>0?k.slice(0,j):k, f=j>0?k.slice(j+1).replace(/_/g," "):"value";
    if(!(db in at)){ at[db]=out.length; out.push([db,[]]); } out[at[db]][1].push([f,v]); });
  return out;
}
function geneSections(i,X,pad){
  const g=D.genes[i], x=X[i]||{}, cn=D.contigs[g.ci].name, P=(n)=>pad+"  ".repeat(n);
  const L=[], kv=(n,k,v)=>{ if(v!==undefined&&v!==null&&String(v).trim()!=="") L.push(`${P(n)}${k}: ${v}`); };
  const sub=(n,k,v)=>kv(n+1,k,v);
  L.push(`${pad}GENE: ${g.nm||"(unnamed)"}`);
  kv(1,"Feature id",g.fid); kv(1,"Gene id",x.gene_id); kv(1,"Type",x.FEATURE_TYPE);
  kv(1,"Location",`${cn}:${g.s.toLocaleString()}-${g.e.toLocaleString()} (${g.st>0?"+":"-"})`);
  kv(1,"Operon",g.op?`${g.op} (P = ${pct(g.pr)})`:"none (singleton)");
  kv(1,"Name chosen from",g.src); kv(1,"Source rank",x.best_consensus_product_descriptor_source_hierarchy_order);
  kv(1,"Envelope",x.ENVELOPE); kv(1,"Envelope reason",x.ENVELOPE_reason);
  L.push(pad+SEP);
  L.push(`${pad}CONFIDENCE`);
  kv(1,"Tier",`${tierBadge(g.ti)[0]} (adj) | ${tierBadge(g.tih)[0]} (hyb)`);
  kv(1,"Final confidence",`${pct(g.fin)} (adj) | ${pct(g.finh)} (hyb)`);
  sub(1,"formula",x.ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_formula); sub(1,"reasoning",x.ADJUSTED_CONFIDENCE_WITH_OPERON_CONTEXT_reasoning);
  kv(1,"Operon context changes it",x["DOES_OPERON_CONTEXT_IMPROVE_CONFIDENCE?"]);
  kv(1,"Preliminary (C1 x C4)",pct(g.pre));
  sub(1,"formula",x.PRELIMINARY_confidence_C1_C4_formula); sub(1,"reasoning",x.PRELIMINARY_confidence_C1_C4_reasoning);
  kv(1,"C1 database coverage",pct(g.c1)); sub(1,"reasoning",x.C1_reasoning);
  kv(1,"C2 operon membership",`${pct(g.c2)} (${g.op?"in an operon":"singleton"})`);
  sub(1,"formula",x.C2_score_formula); sub(1,"reasoning",x.C2_score_reasoning);
  kv(1,"C3 operon context",`${pct(g.c3)} (adj) | ${pct(g.c3h)} (hyb)`);
  sub(1,"formula",x.C3_score_formula); sub(1,"reasoning",x.C3_score_reasoning);
  kv(1,"C4 EC agreement",`${pct(g.c4)}${g.ecs?" | "+g.ecs:""}`);
  sub(1,"formula",x.C4_score_formula); sub(1,"reasoning",x.C4_score_reasoning);
  kv(1,"Review",g.rv?`flagged: ${x.NEEDS_REVIEW_REASON||g.rr||""}`:"not flagged");
  L.push(pad+SEP);
  L.push(`${pad}FINGERPRINT`);
  if(g.fp){
    kv(1,"Label",g.fp.l); kv(1,"Pattern hash",g.fp.h);
    kv(1,"Fingerprint (one line)",x.fingerprint_line||(g.fp.f||[]).map(([k,v])=>`${k}: ${v}`).join(" | "));
    L.push(`${P(1)}By database:`);
    byDatabase(g.fp.f).forEach(([db,items])=>{ L.push(`${P(2)}${db}`); items.forEach(([f,v])=>L.push(`${P(3)}${f}: ${v}`)); });
  } else L.push(`${P(1)}no fingerprint for this gene`);
  L.push(pad+SEP);
  L.push(`${pad}EVIDENCE TRAIL (every database's call)`);
  const ev=g.ev||[]; let last=null;
  if(!ev.length) L.push(`${P(1)}no per-tool record for this feature`);
  ev.forEach(r=>{ const grp=D.evGroups[r[0]], name=D.evNames[r[0]];
    if(grp!==last){ L.push(`${P(1)}${GRP_LABEL[grp]||grp}`); last=grp; }
    L.push(`${P(2)}${name}${r[2]?" ("+r[2]+")":""}: ${r[1]}${isWin(name,g.src)?"  [chosen]":""}${r[3]?"":"  [uninformative]"}`); });
  if(g.up) kv(1,"UniProt best hit",`${g.up.en||""} | ${g.up.pid!=null?g.up.pid+"% id | ":""}${g.up.desc||""}${g.up.inf?"":" (uninformative, not used)"}`);
  kv(1,"Name audit trail",x.best_consensus_product_descriptor_source_audit_trail);
  kv(1,"Specialised database hits",x.specialized_database_hits);
  kv(1,"Localisation and topology",x.localization_and_topology_hits);
  return L;
}
/** The AI conversation about this gene, asked of MARGIE when the page is shown inside it. */
function askAI(i){
  return new Promise(res=>{
    if(window.parent===window){ res(""); return; }
    const g=D.genes[i], x=(XTRA||[])[i]||{}, id="r"+Math.random().toString(36).slice(2);
    const done=(t)=>{ removeEventListener("message",on); res(t); };
    const on=(e)=>{ const d=e.data||{}; if(d.margie==="ai-text"&&d.id===id) done(String(d.text||"")); };
    addEventListener("message",on);
    window.parent.postMessage({margie:"ai-for-gene",id,genome:D.short,fid:g.fid||"",gene:x.gene_id||"",product:g.nm||""},"*");
    setTimeout(()=>done(""),1500);
  });
}
async function geneReport(i){
  const X=await extras(), g=D.genes[i];
  const L=[`MARGIE gene report`,`Genome: ${D.short}`,`Written: ${new Date().toLocaleString()}`,SEP];
  L.push(...geneSections(i,X,""));
  L.push(SEP);
  if(g.op&&operons[g.op]){
    const mem=operons[g.op];
    L.push(`OPERON ${g.op} (P = ${pct(g.pr)}): ${mem.length} genes, in order`);
    mem.forEach((j,k)=>{ L.push(""); L.push(`    [${k+1}/${mem.length}]${j===i?" (this gene)":""}`); L.push(...geneSections(j,X,"    ")); });
  } else L.push(`OPERON: none (singleton)`);
  L.push(SEP);
  const ai=await askAI(i);
  if(ai){ L.push(`AI INTERPRETATION (from Chat with the genome; check it against the evidence above)`); L.push(ai); L.push(SEP); }
  return L.join("\n")+"\n";
}
async function downloadGene(i){
  const g=D.genes[i], name=(g.fid||("gene_"+i)).replace(/[^A-Za-z0-9._-]+/g,"_");
  const a=document.createElement("a");
  a.href=URL.createObjectURL(new Blob([await geneReport(i)],{type:"text/plain"}));
  a.download=`${D.short}_${name}_gene_report.txt`; a.click();
  setTimeout(()=>URL.revokeObjectURL(a.href),1500);
}
document.addEventListener("click",e=>{ const b=e.target.closest("[data-gene-dl]"); if(b) downloadGene(+b.dataset.geneDl); });
const GRP_LABEL={decision:"decision databases (set C1)",domain:"domain / family signatures",special:"specialised callers"};
function isWin(name,src){ if(!src)return false; const a=name.toUpperCase().replace(/[^A-Z]/g,""),b=src.toUpperCase().replace(/[^A-Z]/g,"");
  return a&&b&&(b.indexOf(a)>=0||a.indexOf(b)>=0); }
function evidenceTrail(g){
  const ev=g.ev||[];
  if(!ev.length) return `<div class="trail"><h4>evidence trail</h4><div class="cnt">no per-tool record for this feature.</div></div>`;
  const infN=ev.filter(r=>r[3]).length;
  let h=`<div class="trail"><h4>evidence trail</h4>`
       +`<div class="cnt">${infN} of ${ev.length} databases returned an informative name`
       +`${g.src?` | chosen: <b>${esc(g.src)}</b>`:""}</div>`;
  if(g.c4!=null && g.c4<1)
    h+=`<div class="ecflag">⚑ EC conflict: C4 = ${pct(g.c4)}${g.ecs?" ("+esc(g.ecs)+")":""}. `
      +`${esc(g.c4r||"tools disagree on the EC number; compare the [EC:…] tags below.")}</div>`;
  let last=null;
  ev.forEach(row=>{
    const ti=row[0], desc=row[1], metric=row[2]||"", inf=row[3];
    const grp=D.evGroups[ti], name=D.evNames[ti];
    if(grp!==last){ h+=`<div class="grp">${GRP_LABEL[grp]||grp}</div>`; last=grp; }
    const win=isWin(name,g.src);
    h+=`<div class="erow${win?" win":""}${inf?"":" uninf"}"><span class="etool">${esc(name)}`
      +`${metric?` <span class="emetric">${esc(metric)}</span>`:""}</span>`
      +`<span class="edesc">${esc(desc)}${win?`<span class="chosen">chosen</span>`:""}</span></div>`;
  });
  return h+`</div>`;
}
function operonCard(id){
  const idx=operons[id];
  const sorted=idx.map(i=>[i,D.genes[i]]).sort((a,b)=>a[1].s-b[1].s);
  const span=Math.max(...sorted.map(x=>x[1].e))-Math.min(...sorted.map(x=>x[1].s));
  const raised=sorted.filter(([,g])=>g.imp).length;             // FINAL's own flag, not recomputed
  const flagged=sorted.filter(([,g])=>g.rv).length;
  let h=`<div class="cardhead"><div class="cardhead-l">`
       +`<div class="ptitle">${id}</div>`
       +`<span class="ptag" style="background:${opColor(id)}">${idx.length} genes | operon</span></div>`
       +`<button class="dlbtn" data-dl="${id}" title="Download this operon map as an image">Download map (PNG)</button></div>`;
  h+=`<div class="kv">`;                                        // structural facts only (no recomputed scores)
  h+=`<b>span</b><span>${(span/1000).toFixed(1)} kb</span>`;
  h+=`<b>raised by operon context</b><span>${raised} / ${idx.length} gene${idx.length>1?"s":""}</span>`;
  h+=`<b>flagged for review</b><span>${flagged} / ${idx.length}</span>`;
  h+=`</div>`;
  h+=operonArrowMap(sorted);
  h+=`<div class="cnt">arrows point 5′→3′, coloured by tier | per-gene C1–C4 and final (adj) / (hyb) in the downloaded map | click an arrow to open a gene</div>`;
  h+=`<div class="members">`;
  sorted.forEach(([i,g],k)=>{
    const c=g.ti<0?NONCODE:TIER_COL[g.ti];
    h+=`<div class="mrow" data-goto="${i}"><span class="mnum">${k+1}</span><span class="dot" style="background:${c}"></span>`
      +`<span class="mnm">${esc(g.nm||"(unnamed)")}</span>`
      +`<span class="mfin">${pct(g.fin)}${g.rv?" ⚑":""}</span></div>`;
  });
  h+=`</div>`;
  return h;
}
function operonArrowMap(sorted){
  const lo=Math.min(...sorted.map(x=>x[1].s)), hi=Math.max(...sorted.map(x=>x[1].e));
  const W=380,M=8,H=46,y=15,h=20,rng=Math.max(1,hi-lo);
  const sx=v=>M+(W-2*M)*(v-lo)/rng;
  let s=`<svg viewBox="0 0 ${W} ${H}" class="omap" preserveAspectRatio="xMidYMid meet">`;
  s+=`<line x1="${M}" y1="${y+h/2}" x2="${W-M}" y2="${y+h/2}" stroke="#c2c2bd" stroke-width="1.6"/>`;
  sorted.forEach(([i,g],k)=>{
    let x0=sx(g.s), x1=sx(g.e); if(x1-x0<4){const m=(x0+x1)/2;x0=m-2;x1=m+2;}
    const hd=Math.min(9,(x1-x0)*0.55), c=g.ti<0?NONCODE:TIER_COL[g.ti], yt=y,yb=y+h,ym=y+h/2;
    const pts=g.st>0
      ? `${x0.toFixed(1)},${yt} ${(x1-hd).toFixed(1)},${yt} ${x1.toFixed(1)},${ym} ${(x1-hd).toFixed(1)},${yb} ${x0.toFixed(1)},${yb}`
      : `${x1.toFixed(1)},${yt} ${(x0+hd).toFixed(1)},${yt} ${x0.toFixed(1)},${ym} ${(x0+hd).toFixed(1)},${yb} ${x1.toFixed(1)},${yb}`;
    s+=`<polygon points="${pts}" fill="${c}" stroke="rgba(0,0,0,.2)" stroke-width="0.6" data-goto="${i}">`
      +`<title>${esc((k+1)+". "+(g.nm||"(unnamed)"))} | ${g.st>0?"+":"−"} | final ${pct(g.fin)}${g.rv?" | ⚑":""}</title></polygon>`;
    s+=`<text x="${((x0+x1)/2).toFixed(1)}" y="${y-4}" font-size="9" fill="#8a8a86" text-anchor="middle">${k+1}</text>`;
  });
  return s+`</svg>`;
}
function esc(s){return String(s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));}

const panel=document.getElementById("panel");
function showGene(i){selKind="gene";selVal=i;panel.innerHTML=geneCard(i);applySelection();}
function showOperon(id){selKind="operon";selVal=id;panel.innerHTML=operonCard(id);applySelection();}

// ---- interaction ----
const tip=document.getElementById("tip");
function tipShow(html,x,y){tip.innerHTML=html;tip.style.opacity=1;
  tip.style.left=Math.min(x+14,innerWidth-tip.offsetWidth-8)+"px";
  tip.style.top=(y+14)+"px";}
function tipHide(){tip.style.opacity=0;}

geneLayer.addEventListener("mousemove",e=>{
  const t=e.target.closest(".gene"); if(!t){tipHide();return;}
  const i=+t.dataset.i, g=D.genes[i];
  if(mode==="operon"&&g.op){
    tipShow(`<b>${g.op}</b> | ${operons[g.op].length} genes<br>${esc((g.nm||"").slice(0,60))}`,e.clientX,e.clientY);
  } else {
    const[tn]=tierBadge(g.ti);
    tipShow(`${esc((g.nm||"(unnamed)").slice(0,64))}<br><b>${tn}</b> | final ${pct(g.fin)}${g.rv?" | ⚑":""}`,e.clientX,e.clientY);
  }
});
geneLayer.addEventListener("mouseleave",tipHide);
geneLayer.addEventListener("click",e=>{
  const t=e.target.closest(".gene"); if(!t)return;
  const i=+t.dataset.i, g=D.genes[i];
  if(mode==="operon"&&g.op) showOperon(g.op); else showGene(i);
});
panel.addEventListener("click",e=>{
  const dl=e.target.closest("[data-dl]"); if(dl){ downloadOperon(dl.dataset.dl); return; }
  const ol=e.target.closest("[data-op-link]"); if(ol){ setMode("operon"); showOperon(ol.dataset.opLink); return; }
  const m=e.target.closest("[data-goto]"); if(m) showGene(+m.dataset.goto);
});

// ---- downloadable standalone operon figure (white bg, theme-independent) ----
function revShort(rr){                          // full review sentence -> short trigger tag(s)
  if(!rr) return "yes";
  const s=rr.toLowerCase(), t=[];
  if(s.includes("ec conflict")) t.push("EC conflict");
  if(s.includes("low confidence")) t.push("low conf.");
  if(s.includes("operon inference ambig")) t.push("operon ambig.");
  return t.length?t.join("; "):(rr.length>20?rr.slice(0,18)+"…":rr);
}
function operonFigureSVG(id){
  const sorted=operons[id].map(i=>D.genes[i]).sort((a,b)=>a.s-b.s), n=sorted.length;
  const lo=Math.min(...sorted.map(g=>g.s)), hi=Math.max(...sorted.map(g=>g.e)), rng=Math.max(1,hi-lo);
  const W=1480, MX=32, aw=W-2*MX, FF="'Times New Roman', Times, Georgia, serif", INK="#000000";
  const arrY=104, arrH=30;                                 // arrow band; numbers above, intergenic below
  const headY=arrY+arrH+48, rowH=21;                       // table header baseline
  const legY=headY+24+rowH*n+16, H=legY+26;
  const sx=v=>MX+aw*(v-lo)/rng;
  const XE=s=>String(s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
  const dec=v=>v==null?"—":(+v).toFixed(2);
  const igd=g=>g<=0?"‹1 bp":g<1000?g+" bp":(g/1000).toFixed(1)+" kb";
  const T=(x,y,s,sz,w,anc)=>`<text x="${x}" y="${y}" font-family="${FF}" `
    +`font-size="${sz}" fill="${INK}"${w?` font-weight="${w}"`:''}${anc?` text-anchor="${anc}"`:''}>${XE(s)}</text>`;
  let s=`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`;
  s+=`<rect width="${W}" height="${H}" fill="#ffffff"/>`;
  s+=T(MX,30,id,20,700);
  s+=T(MX,52,D.organism,13);                               // genome identifier / filename, verbatim
  s+=T(MX,72,`${n} genes  |  ${((hi-lo)/1000).toFixed(1)} kb region  |  arrow length ∝ gene length  |  intergenic distances shown below arrows`,12.5);
  s+=`<line x1="${MX}" y1="${arrY+arrH/2}" x2="${W-MX}" y2="${arrY+arrH/2}" stroke="#cccccc" stroke-width="1.2"/>`;
  const ax=[]; let prevEnd=null;
  sorted.forEach((g,k)=>{
    let x0=sx(g.s), x1=sx(g.e); if(x1-x0<7) x1=x0+7;
    if(prevEnd!=null && x0<prevEnd+2){ const w=x1-x0; x0=prevEnd+2; x1=x0+w; }
    prevEnd=x1; ax.push([x0,x1]);
    // Arrows are coloured by POSITION IN THE OPERON, not by tier. Tier is a
    // sequential ramp, so every member of one operon fell in the same hue
    // family and the figure read as monochrome. Identity is unambiguous here
    // because each arrow is numbered and the table below repeats the number,
    // so colour is free to do a different job. Tier is still shown -- as the
    // swatch and named column in that table.
    const hd=Math.min(12,(x1-x0)*0.5),
          c=g.ti<0?NONCODE:OPERON_CYCLE[k%OPERON_CYCLE.length],
          yt=arrY,yb=arrY+arrH,ym=arrY+arrH/2;
    const pts=g.st>0?`${x0.toFixed(1)},${yt} ${(x1-hd).toFixed(1)},${yt} ${x1.toFixed(1)},${ym} ${(x1-hd).toFixed(1)},${yb} ${x0.toFixed(1)},${yb}`
                    :`${x1.toFixed(1)},${yt} ${(x0+hd).toFixed(1)},${yt} ${x0.toFixed(1)},${ym} ${(x0+hd).toFixed(1)},${yb} ${x1.toFixed(1)},${yb}`;
    s+=`<polygon points="${pts}" fill="${c}" stroke="${g.rv?'#c0143c':'#000000'}" stroke-width="${g.rv?1.8:0.7}"/>`;
    s+=T((x0+x1)/2,arrY-10,String(k+1),11,null,'middle');
  });
  for(let k=0;k<n-1;k++)
    s+=T((ax[k][1]+ax[k+1][0])/2,arrY+arrH+15,igd(sorted[k+1].s-sorted[k].e-1),10,null,'middle');
  // table: # | swatch | gene product | location | C1 | C2 | C3 (adj)/(hyb) | C4 | final (adj)/(hyb) | tier | review
  const CX={loc:432, c1:710,c2:770,c3:903,c4:965,fin:1090,tier:1180,rev:W-MX};
  s+=`<line x1="${MX}" y1="${headY+7}" x2="${W-MX}" y2="${headY+7}" stroke="#000000" stroke-width="0.8"/>`;
  s+=T(MX,headY,'#',11.5,700);
  // Two swatch columns precede the product name: arrow colour, then tier.
  s+=T(MX+21.5,headY,'map',8.5,700,'middle');
  s+=T(MX+41.5,headY,'tier',8.5,700,'middle');
  s+=T(MX+58,headY,'gene product',11.5,700);
  s+=T(CX.loc,headY,'location (bp)',11.5,700);
  s+=T(CX.c1,headY,'C1',11.5,700,'end');
  s+=T(CX.c2,headY,'C2',11.5,700,'end');
  s+=T(CX.c3,headY,'C3 adj/hyb',11.5,700,'end');
  s+=T(CX.c4,headY,'C4',11.5,700,'end');
  s+=T(CX.fin,headY,'final adj/hyb',11.5,700,'end');
  s+=T(CX.tier,headY,'tier',11.5,700,'end');
  s+=T(CX.rev,headY,'review',11.5,700,'end');
  sorted.forEach((g,k)=>{
    // Two swatches per row: the arrow's own colour (identity, matches the map
    // above) and the tier colour (confidence). Showing only one of them would
    // make the table disagree with the figure.
    const y=headY+24+rowH*k,
          c=g.ti<0?NONCODE:OPERON_CYCLE[k%OPERON_CYCLE.length],
          ct=g.ti<0?NONCODE:TIER_COL[g.ti],
          tn=g.ti<0?'non-coding':TIER_NAMES[g.ti];
    s+=T(MX,y,String(k+1),11);
    s+=`<rect x="${MX+16}" y="${y-9}" width="11" height="11" fill="${c}" stroke="#000000" stroke-width="0.5"/>`;
    s+=`<rect x="${MX+36}" y="${y-9}" width="11" height="11" fill="${ct}" stroke="#000000" stroke-width="0.5"/>`;
    s+=T(MX+58,y,(g.nm||'(unnamed)').slice(0,54),12);
    s+=T(CX.loc,y,`${g.s.toLocaleString()}–${g.e.toLocaleString()} ${g.st>0?'+':'−'}`,11);
    s+=T(CX.c1,y,dec(g.c1),12,null,'end');
    s+=T(CX.c2,y,dec(g.c2),12,null,'end');
    s+=T(CX.c3,y,`${dec(g.c3)}/${dec(g.c3h)}`,12,null,'end');
    s+=T(CX.c4,y,dec(g.c4),12,null,'end');
    s+=T(CX.fin,y,`${dec(g.fin)}/${dec(g.finh)}`,12,700,'end');
    s+=T(CX.tier,y,tn,11,null,'end');
    s+=T(CX.rev,y,g.rv?revShort(g.rr):'',11,null,'end');
  });
  let lx=MX;
  TIER_NAMES.concat(['non-coding']).forEach((nm,i)=>{
    s+=`<rect x="${lx}" y="${legY}" width="12" height="12" fill="${i<5?TIER_COL[i]:NONCODE}" stroke="#000000" stroke-width="0.5"/>`;
    s+=T(lx+16,legY+10,nm,11); lx+=16+nm.length*6.3+20;
  });
  s+=T(lx+4,legY+10,'red outline = flagged for review',11);
  return s+`</svg>`;
}
function dlSVG(svg,id){
  const a=document.createElement('a');
  a.href='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(svg);
  a.download=id+'_operon_map.svg'; a.click();
}
function downloadOperon(id){
  const svg=operonFigureSVG(id);
  const img=new Image();
  img.onload=()=>{ try{
      const sc=3, cv=document.createElement('canvas');   // 3× → high-resolution PNG (~3360 px wide)
      cv.width=img.width*sc; cv.height=img.height*sc;
      const ctx=cv.getContext('2d'); ctx.setTransform(sc,0,0,sc,0,0); ctx.drawImage(img,0,0);
      cv.toBlob(b=>{ if(!b){dlSVG(svg,id);return;} const a=document.createElement('a');
        a.href=URL.createObjectURL(b); a.download=id+'_operon_map.png'; a.click();
        setTimeout(()=>URL.revokeObjectURL(a.href),1500); },'image/png');
    }catch(err){ dlSVG(svg,id); } };
  img.onerror=()=>dlSVG(svg,id);
  img.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(svg);
}

// ---- linear view: the same arrows as an operon map, along the replicon ----
// The whole replicon is scrollable, but only the genes in view are drawn, so
// 4,000 genes cost no more than 12 polygons at a time.
const LIN=(function(){
  const plateC=document.getElementById("circPlate"), plateL=document.getElementById("linPlate");
  const scroll=document.getElementById("linScroll"), spacer=document.getElementById("linSpacer");
  const svgL=document.getElementById("lin"), ctl=document.getElementById("linctl");
  const selC=document.getElementById("linContig"), selN=document.getElementById("linN"), posEl=document.getElementById("linPos");
  const bC=document.getElementById("vCirc"), bL=document.getElementById("vLin");
  if(!plateL||!scroll||!svgL) return {render(){},visible(){return false;}};

  const byContig=D.contigs.map(()=>[]);
  D.genes.forEach((g,i)=>byContig[g.ci].push(i));
  byContig.forEach(a=>a.sort((x,y)=>D.genes[x].s-D.genes[y].s));
  D.contigs.forEach((c,i)=>{ const o=document.createElement("option");
    o.value=i; o.textContent=c.name+"  |  "+(c.len/1e6).toFixed(2)+" Mb ("+byContig[i].length.toLocaleString()+" genes)";
    selC.appendChild(o); });
  if(D.contigs.length<2) selC.parentElement.hidden=true;

  let ci=0, per=10, on=false, frame=0;
  const TOP=58, AH=52;                         // arrow top, arrow height
  // One line: every gene on one track, pointing its way. Two strands: + genes
  // (read 5'->3' left to right) on the upper track, - genes (5'->3' right to
  // left) on the lower one, each labelled on its outer side.
  let two=false;
  const LOW=TOP+AH+44;                          // the lower track's top, in two-strand view
  const heightNow=()=>two?LOW+AH+50:210;
  const selS=document.getElementById("linStr");
  const slot=()=>Math.max(70, scroll.clientWidth/per);

  function layout(){ spacer.style.width=(byContig[ci].length*slot())+"px"; render(); }

  function render(){
    if(!on) return;
    const list=byContig[ci], sw=slot(), W=scroll.clientWidth||600, H=heightNow();
    svgL.setAttribute("viewBox",`0 0 ${W} ${H}`); svgL.style.height=H+"px";
    const left=scroll.scrollLeft;
    const first=Math.max(0,Math.floor(left/sw)-1), last=Math.min(list.length-1,first+per+2);
    let out=two
      ? `<line x1="0" y1="${TOP+AH/2}" x2="${W}" y2="${TOP+AH/2}" stroke="#d7d7d2" stroke-width="2"/>`
        +`<line x1="0" y1="${LOW+AH/2}" x2="${W}" y2="${LOW+AH/2}" stroke="#d7d7d2" stroke-width="2"/>`
        +`<text x="6" y="${TOP+AH+16}" font-size="11" fill="#555">+ strand (5′→3′, left to right)</text>`
        +`<text x="6" y="${LOW-6}" font-size="11" fill="#555">− strand (5′→3′, right to left)</text>`
      : `<line x1="0" y1="${TOP+AH/2}" x2="${W}" y2="${TOP+AH/2}" stroke="#d7d7d2" stroke-width="2"/>`;
    const short=(nm)=>esc(nm.length>Math.floor(sw/7)?nm.slice(0,Math.max(6,Math.floor(sw/7)-1))+"…":nm);
    for(let k=first;k<=last;k++){
      const i=list[k], g=D.genes[i];
      const x=k*sw-left, w=sw*0.82, x0=x+sw*0.09, x1=x0+w, cx=x+sw/2;
      const top=two&&g.st<0?LOW:TOP;
      const hd=Math.min(16,w*0.3), yt=top, yb=top+AH, ym=top+AH/2;
      const pts=g.st>0
        ? `${x0},${yt} ${x1-hd},${yt} ${x1},${ym} ${x1-hd},${yb} ${x0},${yb}`
        : `${x1},${yt} ${x0+hd},${yt} ${x0},${ym} ${x0+hd},${yb} ${x1},${yb}`;
      const c=geneFill(g), nm=(g.nm||"(unnamed)");
      // Narrow slots get the position only; wider ones the score too.
      const kb=`${(g.s/1000).toFixed(1)}–${(g.e/1000).toFixed(1)} kb`+(sw>=150?` | ${pct(g.fin)}${g.rv?" | ⚑":""}`:"");
      out+=`<polygon class="gene" data-i="${i}" points="${pts}" fill="${c}" stroke="rgba(0,0,0,.35)" stroke-width="0.8"></polygon>`;
      if(!two){
        out+=`<text x="${cx}" y="${TOP-26}" font-size="11" fill="#000" text-anchor="middle">${k+1}</text>`
          +`<text x="${cx}" y="${TOP-10}" font-size="11.5" fill="#000" text-anchor="middle">${short(nm)}</text>`
          +`<text x="${cx}" y="${yb+16}" font-size="10.5" fill="#444" text-anchor="middle">${(g.s/1000).toFixed(1)}–${(g.e/1000).toFixed(1)} kb</text>`
          +`<text x="${cx}" y="${yb+31}" font-size="10.5" fill="#444" text-anchor="middle">${g.st>0?"+":"−"} | ${pct(g.fin)}${g.rv?" | ⚑":""}</text>`;
      } else if(g.st>0){                       // upper track: labels above
        out+=`<text x="${cx}" y="${yt-26}" font-size="11" fill="#000" text-anchor="middle">${k+1} | ${short(nm)}</text>`
          +`<text x="${cx}" y="${yt-10}" font-size="10.5" fill="#444" text-anchor="middle">${kb}</text>`;
      } else {                                  // lower track: labels below
        out+=`<text x="${cx}" y="${yb+16}" font-size="10.5" fill="#444" text-anchor="middle">${kb}</text>`
          +`<text x="${cx}" y="${yb+32}" font-size="11" fill="#000" text-anchor="middle">${k+1} | ${short(nm)}</text>`;
      }
    }
    svgL.innerHTML=out;
    const lo=list[Math.min(first+1,list.length-1)], hi=list[Math.min(last,list.length-1)];
    posEl.textContent=`genes ${Math.min(first+2,list.length).toLocaleString()}–${(last+1).toLocaleString()} of ${list.length.toLocaleString()}`
      +`  |  ${(D.genes[lo].s/1e6).toFixed(3)}–${(D.genes[hi].e/1e6).toFixed(3)} Mb`;
  }

  scroll.addEventListener("scroll",()=>{ if(frame) return; frame=requestAnimationFrame(()=>{frame=0;render();}); });
  addEventListener("resize",()=>{ if(on) layout(); });
  selC.addEventListener("change",()=>{ ci=+selC.value; scroll.scrollLeft=0; layout(); });
  selN.addEventListener("change",()=>{ per=+selN.value||10; layout(); });
  if(selS) selS.addEventListener("change",()=>{ two=selS.value==="two"; layout(); });
  svgL.addEventListener("click",e=>{ const t=e.target.closest(".gene"); if(!t) return;
    const i=+t.dataset.i, g=D.genes[i];
    if(mode==="operon"&&g.op) showOperon(g.op); else showGene(i); });
  svgL.addEventListener("mousemove",e=>{ const t=e.target.closest(".gene"); if(!t){tipHide();return;}
    const g=D.genes[+t.dataset.i], [tn]=tierBadge(g.ti);
    tipShow(`${esc((g.nm||"(unnamed)").slice(0,64))}<br><b>${tn}</b> | final ${pct(g.fin)}${g.op?" | "+esc(g.op):""}`,e.clientX,e.clientY); });
  svgL.addEventListener("mouseleave",tipHide);

  function show(linear){
    on=linear;
    plateC.hidden=linear; plateL.hidden=!linear; ctl.hidden=!linear;
    bC.classList.toggle("on",!linear); bL.classList.toggle("on",linear);
    if(linear) layout();
  }
  /** The linear view at the start of one contig. */
  function showContig(c){ ci=c; selC.value=String(c); scroll.scrollLeft=0; if(on) layout(); }
  /** Bring gene *i* into view, switching replicon if it sits on another. */
  function goto_(i){
    const g=D.genes[i]; if(!on) return;
    if(g.ci!==ci){ ci=g.ci; selC.value=String(ci); layout(); }
    const k=byContig[ci].indexOf(i); if(k<0) return;
    scroll.scrollLeft=Math.max(0,(k-Math.floor(per/2))*slot());
    render();
  }
  return {render, visible:()=>on, goto:goto_, show, showContig};
})();

// ---- the three views: circular, linear, contigs ----
const VIEW=(function(){
  const bC=document.getElementById("vCirc"), bL=document.getElementById("vLin"), bK=document.getElementById("vCtg");
  const plateK=document.getElementById("ctgPlate");
  if(MULTI&&bK) bK.hidden=false;
  const what=document.getElementById("linWhat"); if(what&&MULTI) what.textContent="contig";
  function set(v){
    VIEWNOW=v; tell();
    if(v==="ctg"){ LIN.show(false); document.getElementById("circPlate").hidden=true; plateK.hidden=false; renderContigList(); }
    else { plateK.hidden=true; LIN.show(v==="lin"); }
    bC.classList.toggle("on",v==="circ"); bL.classList.toggle("on",v==="lin"); if(bK) bK.classList.toggle("on",v==="ctg");
  }
  bC.addEventListener("click",()=>set("circ"));
  bL.addEventListener("click",()=>set("lin"));
  if(bK) bK.addEventListener("click",()=>set("ctg"));
  return {set};
})();

// ---- mode toggle ----
function setMode(m){
  mode=m;
  document.getElementById("mGene").classList.toggle("on",m==="gene");
  document.getElementById("mOperon").classList.toggle("on",m==="operon");
  document.getElementById("mReview").classList.toggle("on",m==="review");
  document.getElementById("hint").textContent = m==="gene"
    ? "arcs coloured by confidence tier"
    : m==="review"
    ? "only review-flagged genes are coloured (by tier); the rest are greyed"
    : "each operon takes its own colour | grey = non-operonic | click an operon";
  selKind=null;selVal=null;
  panel.innerHTML=`<div class="empty">${m==="gene"
    ? "Click any gene arc for its full confidence scorecard (C1–C4, preliminary, final, review status)."
    : m==="review"
    ? "Only the "+D.nFlag.toLocaleString()+" genes flagged for review are shown in colour. Click one for its scores; if it sits in an operon, the operon link to jump to that operon."
    : "Click an operon to see its member genes and how genome context changed their confidence."}</div>`;
  paint();
}
document.getElementById("mGene").onclick=()=>setMode("gene");
document.getElementById("mOperon").onclick=()=>setMode("operon");
document.getElementById("mReview").onclick=()=>setMode("review");
document.getElementById("showFlags").onchange=paint;

// ---- init ----

// ---- pan / zoom on the circular map -------------------------------------
// The map is dense at genome scale: 4,600 arcs in a 780px circle means single
// genes are sub-pixel. Zoom is what makes individual operons readable, so it
// drives the viewBox rather than a CSS transform -- vectors stay crisp at any
// magnification and hit-testing keeps working, which a scaled bitmap would lose.
(function(){
  const svg=document.getElementById("map");
  if(!svg) return;
  const BASE={x:0,y:0,w:780,h:780};
  let vb={...BASE};
  // 40x in, and out to a third of the map's size, so the whole circle can be
  // made smaller than the plate rather than only filling it.
  const MIN_W=780/40, MAX_W=780*3;

  function apply(){ svg.setAttribute("viewBox",`${vb.x} ${vb.y} ${vb.w} ${vb.h}`);
    const btn=document.getElementById("zreset");
    if(btn) btn.disabled = (Math.abs(vb.w-BASE.w)<0.5 && Math.abs(vb.x)<0.5 && Math.abs(vb.y)<0.5);
  }
  // Client px -> current viewBox units.
  function toSvg(ev){ const r=svg.getBoundingClientRect();
    return {x:vb.x+(ev.clientX-r.left)/r.width*vb.w, y:vb.y+(ev.clientY-r.top)/r.height*vb.h}; }

  svg.addEventListener("wheel",(e)=>{
    e.preventDefault();
    const p=toSvg(e);
    const k=Math.exp(e.deltaY*0.0015);                 // smooth, direction-preserving
    let w=Math.min(MAX_W,Math.max(MIN_W,vb.w*k));
    const s=w/vb.w;
    // Keep the point under the cursor fixed -- zooming to centre instead makes
    // it impossible to reach a specific operon.
    vb={x:p.x-(p.x-vb.x)*s, y:p.y-(p.y-vb.y)*s, w:w, h:vb.h*s};
    clamp(); apply();
  },{passive:false});

  // Boxing a region: ⬚ makes a plain drag draw the box, and Shift always does.
  const areaBtn=document.getElementById("zarea");
  let areaMode=false, band=null, rect=null;
  if(areaBtn) areaBtn.addEventListener("click",()=>{
    areaMode=!areaMode; areaBtn.classList.toggle("on",areaMode);
    svg.style.cursor=areaMode?"crosshair":"grab";
  });
  const NSU="http://www.w3.org/2000/svg";
  function bandStart(p){
    band={x0:p.x,y0:p.y,x1:p.x,y1:p.y};
    rect=document.createElementNS(NSU,"rect"); rect.setAttribute("class","zsel");
    svg.appendChild(rect); bandDraw();
  }
  function bandDraw(){
    if(!band||!rect) return;
    rect.setAttribute("x",Math.min(band.x0,band.x1)); rect.setAttribute("y",Math.min(band.y0,band.y1));
    rect.setAttribute("width",Math.abs(band.x1-band.x0)); rect.setAttribute("height",Math.abs(band.y1-band.y0));
  }
  function bandEnd(){
    if(!band) return false;
    const w=Math.abs(band.x1-band.x0), h=Math.abs(band.y1-band.y0);
    const x=Math.min(band.x0,band.x1), y=Math.min(band.y0,band.y1);
    if(rect) rect.remove();
    band=null; rect=null;
    if(w<6||h<6) return false;                       // a stray click, not a box
    // Fit the box, keeping the map square: the longer side decides.
    const side=Math.max(Math.min(Math.max(w,h),MAX_W),MIN_W);
    vb={x:x+w/2-side/2, y:y+h/2-side/2, w:side, h:side};
    clamp(); apply();
    return true;
  }

  let drag=null;
  // Capture is deferred until the pointer has actually MOVED. setPointerCapture
  // retargets every subsequent pointer event -- including the click -- to the
  // element holding capture, so capturing on pointerdown meant clicks never
  // reached the individual arcs and the detail panel stopped opening entirely.
  // A plain click must behave exactly as it did before pan existed.
  const DRAG_PX = 4;
  svg.addEventListener("pointerdown",(e)=>{
    if(e.button!==0) return;
    drag={x:e.clientX,y:e.clientY,vx:vb.x,vy:vb.y,moved:false,id:e.pointerId,cap:false,
          box:(areaMode||e.shiftKey)};
  });
  svg.addEventListener("pointermove",(e)=>{
    if(!drag || e.pointerId!==drag.id) return;
    if(drag.box){
      if(!drag.moved){
        if(Math.abs(e.clientX-drag.x)+Math.abs(e.clientY-drag.y) <= DRAG_PX) return;
        drag.moved=true;
        try{ svg.setPointerCapture(drag.id); drag.cap=true; }catch(_){}
        bandStart(toSvg({clientX:drag.x,clientY:drag.y}));
      }
      const p=toSvg(e); band.x1=p.x; band.y1=p.y; bandDraw();
      return;
    }
    if(!drag.moved){
      if(Math.abs(e.clientX-drag.x)+Math.abs(e.clientY-drag.y) <= DRAG_PX) return;
      drag.moved=true;
      // Only now does this become a pan; take capture so the drag survives the
      // pointer leaving the element.
      try{ svg.setPointerCapture(drag.id); drag.cap=true; }catch(_){}
      svg.style.cursor="grabbing";
    }
    const r=svg.getBoundingClientRect();
    const dx=(e.clientX-drag.x)/r.width*vb.w, dy=(e.clientY-drag.y)/r.height*vb.h;
    vb.x=drag.vx-dx; vb.y=drag.vy-dy; clamp(); apply();
  });
  function endDrag(e){
    if(!drag) return;
    const wasDrag=drag.moved, cap=drag.cap, id=drag.id, wasBox=drag.box;
    drag=null;
    if(wasBox) bandEnd();
    svg.style.cursor=areaMode?"crosshair":"grab";
    if(cap){ try{ svg.releasePointerCapture(id); }catch(_){} }
    // Suppress only the click that closes a real drag, so panning across the
    // map does not also select a gene. A click that never moved is untouched.
    if(wasDrag){
      const eat=(ev)=>{ ev.stopPropagation(); ev.preventDefault(); };
      svg.addEventListener("click", eat, {capture:true, once:true});
      setTimeout(()=>svg.removeEventListener("click", eat, {capture:true}), 350);
    }
  }
  svg.addEventListener("pointerup",endDrag);
  svg.addEventListener("pointercancel",endDrag);
  svg.addEventListener("dblclick",(e)=>{ e.preventDefault(); vb={...BASE}; apply(); });

  // Never let the map be dragged entirely off screen.
  function clamp(){
    const pad=vb.w*0.5;
    vb.x=Math.min(Math.max(vb.x,-pad),BASE.w-vb.w+pad);
    vb.y=Math.min(Math.max(vb.y,-pad),BASE.h-vb.h+pad);
  }

  // Buttons zoom about the centre of what is currently shown, so repeated
  // clicks stay on whatever the user has already panned to.
  function zoomBy(k){
    const cx=vb.x+vb.w/2, cy=vb.y+vb.h/2;
    const w=Math.min(MAX_W,Math.max(MIN_W,vb.w*k)), s=w/vb.w;
    vb={x:cx-vb.w*s/2, y:cy-vb.h*s/2, w:w, h:vb.h*s};
    clamp(); apply();
  }
  const zin=document.getElementById("zin"), zout=document.getElementById("zout");
  if(zin) zin.addEventListener("click",()=>zoomBy(1/1.4));
  if(zout) zout.addEventListener("click",()=>zoomBy(1.4));
  const btn=document.getElementById("zreset");
  if(btn) btn.addEventListener("click",()=>{ vb={...BASE}; apply(); });
  svg.style.cursor="grab";
  apply();
})();

document.getElementById("org").textContent=D.short;   // genome identifier / filename, verbatim
document.getElementById("sub").textContent=
  "confidence genome viewer  |  "+(D.totLen/1e6).toFixed(2)+" Mb  |  "
  +D.genes.length.toLocaleString()+" genes  |  "+D.contigs.length.toLocaleString()
  +(DRAFT?" contigs (draft assembly)":" replicon"+(D.contigs.length>1?"s":""))
  +"  |  "+D.nOperons.toLocaleString()+" operons  |  "+D.nFlag.toLocaleString()+" flagged";
readShow();
addEventListener("hashchange",()=>{ readShow(); paint(); });
setMode("gene");
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
