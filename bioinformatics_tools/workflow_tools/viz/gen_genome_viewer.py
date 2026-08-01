#!/usr/bin/env python3
"""gen_genome_viewer.py -- build a SELF-CONTAINED interactive HTML genome viewer
from a FINAL_ANNOTATION_WITH_CONFIDENCE.tsv.

Two modes, toggled in the page:
  * Gene mode   -- every gene arc is coloured by CONFIDENCE_TIER (5 distinct,
                   colour-blind-safe Okabe-Ito hues); review flags shown grey.
  * Operon mode -- operonic genes are blue, non-operonic genes light brown,
                   non-coding grey; hovering/clicking an operon highlights all
                   its member genes and shows the operon's details.

Everything is read from the FINAL table; the data is embedded directly in the
HTML so the file opens offline with no server and no external assets.

Usage: gen_genome_viewer.py <FINAL.tsv> <out.html>
"""
import csv
import json
import re
import sys
from pathlib import Path

csv.field_size_limit(10 ** 8)

TIERS = ["highest", "high", "medium", "fair", "low"]          # index 0..4; -1 = non-coding
TIER_IDX = {t: i for i, t in enumerate(TIERS)}

# Evidence trail: what each database actually called this gene. Read from the
# consolidated per-tool matrix, joined by feature_id. Grouped so the panel can
# show the seven C1 decision databases first, then domain signatures, then the
# specialised callers. (display name, group, candidate columns).
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
    try:
        x = float(v)
    except (TypeError, ValueError):
        return ""
    if x <= 1.0:
        x *= 100
    return f"{round(x)}% id"


def col(row, name):
    for k in row:
        if re.sub(r"^Column-[A-Z]+:\s*", "", k or "").strip().lower() == name.lower():
            return row[k] or ""
    return ""


def clean_ev(d):
    d = re.sub(r"^[A-Za-z][\w /()]*?:\s*", "", str(d)).strip()      # drop "JCVI:"/"KEGG:" prefixes
    d = re.sub(r"^gnl\|[^|]*\|[^|]*\|\S*\s*", "", d).strip()         # drop gnl|DB|acc| prefixes
    return d[:120]                                                   # keep [EC:...] tags — the EC evidence


def contig_of(gid):
    return re.sub(r"_[0-9]+[+-][0-9]+$", "", gid)


def num(v, nd=2):
    try:
        return round(float(v), nd)
    except (TypeError, ValueError):
        return None


# Where the consolidated per-tool matrix sits, relative to the FINAL table's
# own directory. The layout differs depending on WHEN this runs:
#   during the run (phase11 scoring)  -> <organism>/consolidation/...
#                                        with FINAL at <organism>/scoring/
#   after reorganize_outputs.py       -> <organism>/per-tool-phased-output/
#                                        consolidation/... with FINAL at the top
# Generating per-organism at scoring time means the post-reorganize guess alone
# would silently miss (the lookup is existence-guarded), producing a viewer with
# an empty evidence trail. Search both, and let --consolidated override.
_CONS_NAME = "consolidated-merged-all-columns.tsv"
_CONS_CANDIDATES = (
    Path("per-tool-phased-output") / "consolidation" / _CONS_NAME,  # reorganized
    Path("consolidation") / _CONS_NAME,                             # mid-run, FINAL at top
    Path("..") / "consolidation" / _CONS_NAME,                      # mid-run, FINAL in scoring/
)


def find_consolidated(final_path, explicit=None):
    """Resolve the consolidated matrix, or None. An explicit path that does not
    exist is an error rather than a silent downgrade -- if the caller named it,
    they expect the evidence trail."""
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


def main():
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
        # Not fatal: the map is complete without it, only the per-gene evidence
        # trail is empty. Say so loudly rather than shipping a hollow viewer.
        print(f"WARNING: no {_CONS_NAME} near {final}; "
              "evidence trail will be empty", file=sys.stderr)

    def evidence_of(fid):
        cr = cons.get(fid)
        if not cr:
            return []
        # Show EVERY tool that returned anything, not just the informative ones,
        # so nothing (including EC numbers) is hidden. row = [ti, desc, metric, inf].
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
        """UniProt's best BLAST hit, shown for EVERY gene (even when its call was
        uninformative and therefore not selected) so the % identity is visible."""
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

    genes = []
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
        })

    n_op = len({g["op"] for g in genes if g["op"]})
    n_flag = sum(g["rv"] for g in genes)
    # display label = the genome identifier / filename verbatim (a user genome may be
    # "abc.fasta" with no parseable scientific name — the filename always works).
    short = organism

    data = {
        "organism": organism, "short": short,
        "contigs": contigs, "genes": genes,
        "nOperons": n_op, "nFlag": n_flag,
        "totLen": sum(clen.values()),
        "evNames": EV_NAMES, "evGroups": EV_GROUPS,
    }
    html = TEMPLATE.replace("/*__DATA__*/", json.dumps(data, separators=(",", ":")))

    if "--artifact" in sys.argv:
        # Artifact host supplies its own <!doctype>/<html>/<head>/<body>; emit only
        # the page content (title + style + body inner) so nothing is duplicated.
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
    /* bright confidence-tier palette */
    --t0:#0b2842; --t1:#154064; --t2:#256291; --t3:#4184b5; --t4:#6ba3c8; --tn:#c8c8c8;
    --operon:#1667e0; --nonop:#c79a5c; --flag:#666666;
  }
  *{box-sizing:border-box}
  html,body{margin:0;background:#ffffff;color:#000000;
    font-family:Calibri,"Segoe UI",Arial,"Helvetica Neue",Helvetica,sans-serif;
    -webkit-font-smoothing:antialiased}
  .wrap{display:flex;flex-wrap:wrap;gap:16px;padding:18px;max-width:1280px;margin:0 auto}
  .left{flex:1 1 560px;min-width:340px}
  .right{flex:1 1 360px;min-width:300px;max-width:440px}
  h1{font-size:19px;margin:0 0 2px;overflow-wrap:anywhere;word-break:break-word}
  h1 em{font-style:italic}
  .sub{color:var(--muted);font-size:14px;margin-bottom:12px}
  .modebar{display:flex;gap:0;margin:6px 0 4px;border:1px solid var(--line);border-radius:9px;overflow:hidden;width:max-content;box-shadow:var(--shadow)}
  .modebar button{font-family:inherit;font-size:15px;padding:7px 20px;border:0;background:var(--btn);color:var(--btn-ink);cursor:pointer;transition:background .12s}
  .modebar button.on{background:var(--ink);color:var(--surface)}
  .opts{font-size:13px;color:var(--muted);margin:6px 0 2px;display:flex;gap:14px;align-items:center}
  .opts label{cursor:pointer;user-select:none}
  .plate{background:#ffffff;border:1px solid var(--line);border-radius:14px;padding:6px;box-shadow:var(--shadow)}
  svg{width:100%;height:auto;display:block;touch-action:none}
  .gene{cursor:pointer}
  .sel{stroke:#111;stroke-width:2;paint-order:stroke}
  .dim{opacity:.16}
  .hi{stroke:#111;stroke-width:1.1}
  .legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:10px;font-size:13.5px}
  .legend span{display:inline-flex;align-items:center;gap:6px;color:var(--ink)}
  .sw{width:14px;height:14px;border-radius:3px;display:inline-block;border:1px solid rgba(0,0,0,.12)}
  .panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px;min-height:280px}
  .panel .empty{color:var(--muted);font-size:14px;line-height:1.5}
  .ptitle{font-size:18px;margin:0 0 2px;line-height:1.25}
  .ptag{display:inline-block;font-size:12px;padding:2px 9px;border-radius:20px;color:#fff;margin-bottom:8px}
  .cardhead{display:flex;justify-content:space-between;align-items:flex-start;gap:12px}
  .cardhead-l{min-width:0}
  .dlbtn{flex:0 0 auto;font-family:inherit;font-size:12.5px;padding:5px 11px;border:1px solid var(--line);
    background:var(--btn);color:var(--btn-ink);border-radius:8px;cursor:pointer;white-space:nowrap;transition:background .12s}
  .dlbtn:hover{background:var(--ink);color:var(--surface);border-color:var(--ink)}
  .kv{display:grid;grid-template-columns:auto 1fr;gap:3px 12px;font-size:14px;margin:8px 0}
  .kv b{color:var(--muted);font-weight:normal}
  .bar{height:9px;border-radius:0;background:#ffffff;border:1px solid #cfcfcf;overflow:hidden;margin:2px 0}
  .bar>i{display:block;height:100%}
  .members{margin-top:10px;max-height:260px;overflow:auto;border-top:1px solid var(--line);padding-top:8px}
  .mrow{display:flex;align-items:center;gap:8px;font-size:13px;padding:3px 6px;border-radius:6px;cursor:pointer}
  .mrow:hover{background:rgba(127,127,127,.16)}
  .dot{width:10px;height:10px;border-radius:50%;flex:0 0 auto}
  .mrow .mnm{flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .mrow .mfin{color:var(--muted);font-variant-numeric:tabular-nums}
  .mrow .mnum{width:16px;flex:0 0 auto;text-align:right;color:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}
  .omap{width:100%;height:auto;display:block;margin:8px 0 2px}
  .omap polygon{cursor:pointer}
  .flagtag{color:var(--flag);font-size:12.5px;margin-top:6px}
  .trail{margin-top:12px;border-top:1px solid var(--line);padding-top:9px;max-height:340px;overflow:auto}
  .trail h4{font-size:12px;margin:0 0 2px;color:var(--muted);font-weight:normal;letter-spacing:.5px;text-transform:uppercase}
  .trail .cnt{color:var(--muted);font-size:12px;margin-bottom:4px}
  .grp{font-size:11.5px;color:var(--muted);margin:9px 0 2px;font-style:italic}
  .erow{display:grid;grid-template-columns:88px 1fr;gap:8px;font-size:12.5px;padding:1.5px 0;line-height:1.3}
  .erow .etool{color:var(--muted);white-space:nowrap}
  .erow.win .etool{color:var(--operon);font-weight:bold}
  .erow.win .edesc{font-weight:bold}
  .chosen{display:inline-block;font-size:9.5px;color:#fff;background:var(--operon);border-radius:10px;padding:1px 7px;margin-left:7px;vertical-align:1px}
  .emetric{display:inline-block;font-size:10.5px;color:var(--operon);font-variant-numeric:tabular-nums}
  .erow.uninf{opacity:.5}
  .ecflag{font-size:12.5px;color:#c0143c;margin:6px 0 2px;line-height:1.35}
  .uphit{margin-top:8px;font-size:12.5px;color:var(--muted);line-height:1.4}
  .uphit b{color:var(--ink);font-weight:normal}
  .uphit .uninf{color:#c0143c}
  .oplink{color:var(--operon);cursor:pointer;text-decoration:underline;text-underline-offset:2px}
  #tip{position:fixed;pointer-events:none;background:#ffffff;color:#000000;font-size:12.5px;
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
    <div class="plate"><svg id="map" viewBox="0 0 780 780" aria-label="circular genome map"></svg></div>
    <div class="legend" id="legend"></div>
    <div class="foot">Every value is read verbatim from FINAL_ANNOTATION_WITH_CONFIDENCE.tsv. Click a gene or operon for details.</div>
  </div>
  <div class="right">
    <div class="panel" id="panel"><div class="empty">Hover to preview, click to pin.<br><br>In <b>Gene mode</b>, each arc is a gene coloured by its confidence tier. In <b>Operon mode</b>, blue arcs are operonic genes and light-brown arcs are non-operonic; selecting an operon shows its member genes and its context scores.</div></div>
  </div>
</div>
<div id="tip"></div>
<script>
const D = /*__DATA__*/;
const TIER_NAMES=["highest","high","medium","fair","low"];
// Ordered tiers -> sequential single-hue ramp, darkest = highest confidence.
// Validated ordinal (monotone L, gaps >= 0.06, 10 deg hue spread, light end
// 2.66:1 on white). The old rainbow failed all four ordinal checks and was
// indistinguishable in grayscale. Keep in sync with TIER in make_circular_genome.py.
const TIER_COL=["#0b2842","#154064","#256291","#4184b5","#6ba3c8"];
// FLAG is a RESERVED status colour -- never reused as a tier step, so an
// alarm can never be confused with a ranking. Operon/non-operon are a separate
// categorical pair used only in operon mode, where the tier ramp is not shown.
const NONCODE="#c8c8c8", OPERON="#1667e0", NONOP="#c79a5c", FLAG="#b32b1e";
const CX=390, CY=396, GAP=4, MINSPAN=4, START=90;
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
    svg.appendChild(el("polygon",{points:`${P(ro,a0)} ${P(ro,a1)} ${P(ri,a1)} ${P(ri,a0)}`,fill,stroke:"none"}));
  }
}
order.forEach(ci=>{
  annulusStrip(ci,R.bbI,R.bbO,"#000000");
  const L=D.contigs[ci].len, step=L>3e6?1e6:5e5;
  for(let p=0;p<=L;p+=step){
    const a=ang(ci,p),[x0,y0]=pol(R.bbO,a),[x1,y1]=pol(R.bbO+8,a);
    svg.appendChild(el("line",{x1:x0,y1:y0,x2:x1,y2:y1,stroke:"#000000","stroke-width":.8}));
    if(L>=25e4 && p%step===0){const[tx,ty]=pol(R.bbO+20,a);
      const t=el("text",{x:tx,y:ty,"font-size":10,fill:"#000000","text-anchor":"middle","dominant-baseline":"middle"});
      t.textContent=(p/1e6).toFixed(1); svg.appendChild(t);}
  }
  if(order.length>1){const mid=lay[ci].start-lay[ci].span/2,[lx,ly]=pol(R.bbO+40,mid);
    const t=el("text",{x:lx,y:ly,"font-size":11,fill:"#000000","text-anchor":"middle","dominant-baseline":"middle"});
    t.textContent=D.contigs[ci].name+" · "+(D.contigs[ci].len/1e6).toFixed(2)+" Mb"; svg.appendChild(t);}
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
cInfo.textContent=(D.totLen/1e6).toFixed(2)+" Mb · "+D.genes.length.toLocaleString()+" genes"; svg.appendChild(cInfo);
const cInfo2=el("text",{x:CX,y:CY+22,"font-size":12.5,fill:"#000000","text-anchor":"middle"});
svg.appendChild(cInfo2);

// ---- colouring ----
let mode="gene", selKind=null, selVal=null;
function geneFill(g){
  if(mode==="gene") return g.ti<0?NONCODE:TIER_COL[g.ti];
  if(mode==="review") return g.rv?(g.ti<0?NONCODE:TIER_COL[g.ti]):"#ececea";
  return g.op?OPERON:(g.ti<0?NONCODE:NONOP);
}
function paint(){
  D.genes.forEach((g,i)=>nodes[i].setAttribute("fill",geneFill(g)));
  document.getElementById("showFlags").checked ? flagLayer.style.display="" : flagLayer.style.display="none";
  cInfo2.textContent = mode==="operon"
    ? D.nOperons.toLocaleString()+" operons"
    : D.nFlag.toLocaleString()+" flagged for review";
  renderLegend(); applySelection();
}
function renderLegend(){
  const L=document.getElementById("legend"); L.innerHTML="";
  const items = mode==="operon"
    ? [["operonic gene",OPERON],["non-operonic",NONOP],["non-coding",NONCODE]]
    : mode==="review"
    ? TIER_NAMES.map((n,i)=>["flagged · "+n,TIER_COL[i]]).concat([["not flagged","#ececea"]])
    : TIER_NAMES.map((n,i)=>[n,TIER_COL[i]]).concat([["non-coding",NONCODE]]);
  items.forEach(([n,c])=>{const s=document.createElement("span");
    s.innerHTML=`<i class="sw" style="background:${c}"></i>${n}`; L.appendChild(s);});
  if(document.getElementById("showFlags").checked){const s=document.createElement("span");
    s.innerHTML=`<i class="sw" style="background:${FLAG}"></i>review flag`; L.appendChild(s);}
}

// ---- selection / highlight ----
function clearFX(){nodes.forEach(n=>n.classList.remove("sel","dim","hi"));}
function applySelection(){
  clearFX();
  if(selKind==="operon" && operons[selVal]){
    const set=new Set(operons[selVal]);
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
  let h=`<div class="ptitle">${esc(g.nm||"(unnamed)")}</div>`;
  h+=`<span class="ptag" style="background:${tc}">${tn}</span>`;
  h+=`<div class="kv">`;
  h+=`<b>location</b><span>${cn}:${g.s.toLocaleString()}–${g.e.toLocaleString()} (${strand})</span>`;
  h+=`<b>operon</b><span>${g.op?`<a class="oplink" data-op-link="${g.op}">${g.op}</a> · P=`+pct(g.pr):"none (singleton)"}</span>`;
  h+=`</div>`;
  h+=`<div class="kv">`;
  h+=`<b>C1 database coverage</b><span>${pct(g.c1)}</span>`;
  h+=`<b>C2 operon membership</b><span>${pct(g.c2)}</span>`;
  h+=`<b>C3 operon context</b><span>adj ${pct(g.c3)} · hyb ${pct(g.c3h)}</span>`;
  h+=`<b>C4 EC agreement</b><span>${pct(g.c4)}${g.ecs?" · "+esc(g.ecs):""}</span>`;
  h+=`</div>`;
  h+=`<div class="kv"><b>preliminary (C1 × C4)</b><span>${pct(g.pre)}</span></div>${bar(g.pre,"#888")}`;
  h+=`<div class="kv"><b>final confidence</b><span>adj ${pct(g.fin)} · hyb ${pct(g.finh)}</span></div>${bar(g.fin,tc)}`;
  h+=`<div class="kv"><b>confidence tier</b><span>adj ${tierBadge(g.ti)[0]} · hyb ${tierBadge(g.tih)[0]}</span></div>`;
  if(g.up)h+=`<div class="uphit"><b>UniProt best hit</b> ${g.up.pid!=null?g.up.pid+"% id":"—"}`
    +` · ${esc(g.up.en||"")} · ${esc(g.up.desc||"")}${g.up.inf?"":` <span class="uninf">(uninformative — not used)</span>`}</div>`;
  if(g.rv)h+=`<div class="flagtag">⚑ flagged for review — ${esc(g.rr||"")}</div>`;
  h+=evidenceTrail(g);
  return h;
}
const GRP_LABEL={decision:"decision databases (set C1)",domain:"domain / family signatures",special:"specialised callers"};
function isWin(name,src){ if(!src)return false; const a=name.toUpperCase().replace(/[^A-Z]/g,""),b=src.toUpperCase().replace(/[^A-Z]/g,"");
  return a&&b&&(b.indexOf(a)>=0||a.indexOf(b)>=0); }
function evidenceTrail(g){
  const ev=g.ev||[];
  if(!ev.length) return `<div class="trail"><h4>evidence trail</h4><div class="cnt">no per-tool record for this feature.</div></div>`;
  const infN=ev.filter(r=>r[3]).length;
  let h=`<div class="trail"><h4>evidence trail — every database's call (EC numbers kept)</h4>`
       +`<div class="cnt">${infN} of ${ev.length} databases returned an informative name`
       +`${g.src?` · chosen: <b>${esc(g.src)}</b>`:""}</div>`;
  if(g.c4!=null && g.c4<1)
    h+=`<div class="ecflag">⚑ EC conflict — C4 = ${pct(g.c4)}${g.ecs?" ("+esc(g.ecs)+")":""}. `
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
       +`<span class="ptag" style="background:${OPERON}">${idx.length} genes · operon</span></div>`
       +`<button class="dlbtn" data-dl="${id}" title="Download this operon map as an image">⤓ map</button></div>`;
  h+=`<div class="kv">`;                                        // structural facts only (no recomputed scores)
  h+=`<b>span</b><span>${(span/1000).toFixed(1)} kb</span>`;
  h+=`<b>raised by operon context</b><span>${raised} / ${idx.length} gene${idx.length>1?"s":""}</span>`;
  h+=`<b>flagged for review</b><span>${flagged} / ${idx.length}</span>`;
  h+=`</div>`;
  h+=operonArrowMap(sorted);
  h+=`<div class="cnt">arrows point 5′→3′, coloured by tier · per-gene C1–C4 / final adj·hyb in the downloaded map · click an arrow to open a gene</div>`;
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
      +`<title>${esc((k+1)+". "+(g.nm||"(unnamed)"))} · ${g.st>0?"+":"−"} · final ${pct(g.fin)}${g.rv?" · ⚑":""}</title></polygon>`;
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
    tipShow(`<b>${g.op}</b> · ${operons[g.op].length} genes<br>${esc((g.nm||"").slice(0,60))}`,e.clientX,e.clientY);
  } else {
    const[tn]=tierBadge(g.ti);
    tipShow(`${esc((g.nm||"(unnamed)").slice(0,64))}<br><b>${tn}</b> · final ${pct(g.fin)}${g.rv?" · ⚑":""}`,e.clientX,e.clientY);
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
  const W=1480, MX=32, aw=W-2*MX, FF="Calibri, Arial, Helvetica, sans-serif", INK="#000000";
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
  s+=T(MX,72,`${n} genes  ·  ${((hi-lo)/1000).toFixed(1)} kb region  ·  arrow length ∝ gene length  ·  intergenic distances shown below arrows`,12.5);
  s+=`<line x1="${MX}" y1="${arrY+arrH/2}" x2="${W-MX}" y2="${arrY+arrH/2}" stroke="#cccccc" stroke-width="1.2"/>`;
  const ax=[]; let prevEnd=null;
  sorted.forEach((g,k)=>{
    let x0=sx(g.s), x1=sx(g.e); if(x1-x0<7) x1=x0+7;
    if(prevEnd!=null && x0<prevEnd+2){ const w=x1-x0; x0=prevEnd+2; x1=x0+w; }
    prevEnd=x1; ax.push([x0,x1]);
    const hd=Math.min(12,(x1-x0)*0.5), c=g.ti<0?NONCODE:TIER_COL[g.ti], yt=arrY,yb=arrY+arrH,ym=arrY+arrH/2;
    const pts=g.st>0?`${x0.toFixed(1)},${yt} ${(x1-hd).toFixed(1)},${yt} ${x1.toFixed(1)},${ym} ${(x1-hd).toFixed(1)},${yb} ${x0.toFixed(1)},${yb}`
                    :`${x1.toFixed(1)},${yt} ${(x0+hd).toFixed(1)},${yt} ${x0.toFixed(1)},${ym} ${(x0+hd).toFixed(1)},${yb} ${x1.toFixed(1)},${yb}`;
    s+=`<polygon points="${pts}" fill="${c}" stroke="${g.rv?'#c0143c':'#000000'}" stroke-width="${g.rv?1.8:0.7}"/>`;
    s+=T((x0+x1)/2,arrY-10,String(k+1),11,null,'middle');
  });
  for(let k=0;k<n-1;k++)
    s+=T((ax[k][1]+ax[k+1][0])/2,arrY+arrH+15,igd(sorted[k+1].s-sorted[k].e-1),10,null,'middle');
  // table: # · swatch · gene product · location · C1 · C2 · C3 adj/hyb · C4 · final adj/hyb · tier · review
  const CX={loc:432, c1:710,c2:770,c3:903,c4:965,fin:1090,tier:1180,rev:W-MX};
  s+=`<line x1="${MX}" y1="${headY+7}" x2="${W-MX}" y2="${headY+7}" stroke="#000000" stroke-width="0.8"/>`;
  s+=T(MX,headY,'#',11.5,700);
  s+=T(MX+40,headY,'gene product',11.5,700);
  s+=T(CX.loc,headY,'location (bp)',11.5,700);
  s+=T(CX.c1,headY,'C1',11.5,700,'end');
  s+=T(CX.c2,headY,'C2',11.5,700,'end');
  s+=T(CX.c3,headY,'C3 adj/hyb',11.5,700,'end');
  s+=T(CX.c4,headY,'C4',11.5,700,'end');
  s+=T(CX.fin,headY,'final adj/hyb',11.5,700,'end');
  s+=T(CX.tier,headY,'tier',11.5,700,'end');
  s+=T(CX.rev,headY,'review',11.5,700,'end');
  sorted.forEach((g,k)=>{
    const y=headY+24+rowH*k, c=g.ti<0?NONCODE:TIER_COL[g.ti], tn=g.ti<0?'non-coding':TIER_NAMES[g.ti];
    s+=T(MX,y,String(k+1),11);
    s+=`<rect x="${MX+16}" y="${y-9}" width="11" height="11" fill="${c}" stroke="#000000" stroke-width="0.5"/>`;
    s+=T(MX+40,y,(g.nm||'(unnamed)').slice(0,56),12);
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
  s+=T(lx+4,legY+10,'· red outline = flagged for review',11);
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
    : "blue = operonic · brown = non-operonic · click an operon";
  selKind=null;selVal=null;
  panel.innerHTML=`<div class="empty">${m==="gene"
    ? "Click any gene arc for its full confidence scorecard (C1–C4, preliminary, final, review status)."
    : m==="review"
    ? "Only the "+D.nFlag.toLocaleString()+" genes flagged for review are shown in colour. Click one for its scorecard — if it sits in an operon, use the operon link to jump to that operon."
    : "Click a blue operon to see its member genes and how genome context changed their confidence."}</div>`;
  paint();
}
document.getElementById("mGene").onclick=()=>setMode("gene");
document.getElementById("mOperon").onclick=()=>setMode("operon");
document.getElementById("mReview").onclick=()=>setMode("review");
document.getElementById("showFlags").onchange=paint;

// ---- init ----
document.getElementById("org").textContent=D.short;   // genome identifier / filename, verbatim
document.getElementById("sub").textContent=
  "confidence genome viewer  ·  "+(D.totLen/1e6).toFixed(2)+" Mb  ·  "
  +D.genes.length.toLocaleString()+" genes  ·  "+D.contigs.length+" replicon"+(D.contigs.length>1?"s":"")
  +"  ·  "+D.nOperons.toLocaleString()+" operons  ·  "+D.nFlag.toLocaleString()+" flagged";
setMode("gene");
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
