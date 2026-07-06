#!/usr/bin/env python3
"""score-c1-tool-coverage.py — margie_sb phase11 (scoring), metric C1:
tool coverage.

Reads labeled-genes.tsv and labeled-genes-cluster-agreement.tsv (both
phase10, READ-ONLY) and computes, per gene, how many of 7 INDEPENDENT
evidence sources gave an INFORMATIVE hit -- not just any hit. A bare
"Domain of unknown function" counts as a tool finding nothing useful,
same is_uninformative() gate assign-canonical-label.py already applies
when picking the best hit per tool.

C1 = informative_source_count / 7

7 independent sources, justified empirically on the 1,097-gene
calibration set:

  RAST, KEGG, EGGNOG, COG, PFAM, UNIPROT — each independently curated
  databases with distinct methodologies and zero mutual dependency.

  TIGRFAM_CLUSTER — PGAP, TIGRFAM, and NCBIfam collapsed into one slot.
  NCBIfam absorbed TIGRFAM outright; PGAP's HMM library is built from
  NCBIfam models. Empirically: 0 genes where TIGRFAM hits but NCBIFAM
  misses; only 1 gene where PGAP hits but NCBIFAM misses. All three are
  one independent signal.

Excluded from the denominator:
  HAMAP, PIRSF — empirically 0 unique hits beyond COG+PFAM+EGGNOG+NCBIFAM
  on the calibration set; their models overlap completely with those tools.
  CDD — a meta-database aggregating PFAM/TIGRFAM/COG models; by
  construction not independent.
  GENEPROP — uses TIGRFAM HMMs internally; every GENEPROP hit is already
  a TIGRFAM/NCBIFAM hit. Contributes pathway-level context (relevant to
  C3) but not independent gene-level coverage.
  MEROPS/TCDB/DBCAN — narrow specialist DBs that structurally cannot hit
  most genes; excluded to avoid capping non-specialist genes below 1.0.

Output (labeled-genes-c1-tool-coverage.tsv): identity columns, c1_score,
c1_informative_tool_count, c1_total_tools_considered, c1_informative_tools
(which specific sources counted, for traceability -- not just the bare
number), and c1_formula (the literal arithmetic as text, e.g.
"6/7 = 0.8571", so the score is auditable from this column alone).
"""
import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# Standalone, mutually-independent decision tools -- everything except
# the TIGRFAM/PGAP/NCBIfam cluster, which collapses to one slot below.
STANDALONE_DECISION_TOOLS = ["RAST", "COG", "PFAM", "KEGG", "EGGNOG", "UNIPROT"]
TIGRFAM_CLUSTER_SLOT = "TIGRFAM_CLUSTER"
DECISION_TOOLS = STANDALONE_DECISION_TOOLS + [TIGRFAM_CLUSTER_SLOT]

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]

# ─────────────────────────────────────────────────────────────────────────────
# UNINFORMATIVE HIT CATEGORY  (the single, reproducible exclusion list)
# ─────────────────────────────────────────────────────────────────────────────
# A tool's best hit counts toward C1 ONLY if it names a real function. The rule
# groups below are the *complete, deterministic* category of hit descriptions
# that convey NO function and are therefore EXCLUDED: a tool whose only hit
# matches any group here contributes 0 to the /7, not 1. Every exclusion is
# reproducible from source -- no external list, no ordering dependence.
#
# Kept identical to the gate in labeling/assign-canonical-label.py; if you edit
# one, edit both (assign-canonical-label uses it to pick the best-informative
# hit per tool, this file uses it to count). Groups 4-7 were added after an
# audit found ~205 hits per genome-set leaking in as "informative" -- almost all
# eggNOG "Belongs to the UPF#### family" (UPF = Uncharacterized Protein Family),
# bare DUF tags, and mid-string "of unknown function" phrasings that the
# start-anchored Groups 1-2/8 alone could not catch.
#
# NOTE on the "no function word" guard (Groups 3 & 7): only a molecular/cellular
# ACTIVITY word (…ase, transport, kinase, regulator, hydrolase, …) rescues a hit
# that otherwise carries an unknown-function marker. Pure localization words
# (membrane, secreted, periplasmic) are NOT function words -- so "secreted repeat
# of unknown function" stays excluded while "alpha/beta hydrolase of unknown
# function" is kept. Per project decision, topology-only ("predicted membrane
# protein (DUF2238)") and fold-only ("cupin superfamily") hits carry no
# unknown-function marker at all, so they are NOT in this category and still
# count toward the /7.

# Group 1 — exact null / boilerplate tokens.
_UNINFORMATIVE = frozenset({
    "", "-", ".", "na", "n/a", "none", "null",
    "unknown", "uncharacterized", "uncharacterised", "putative", "predicted",
    "hypothetical protein", "conserved hypothetical protein", "conserved protein",
    "conserved domain protein", "predicted protein", "function unknown",
    "domain of unknown function", "general function prediction only",
    "poorly characterized", "open reading frame",
})
# Group 2 — leading-phrase families (unknown-function nouns, DUF/UPF, ORF, …).
_UNINFORMATIVE_PREFIXES = (
    "domain of unknown function", "protein of unknown function",
    "family of unknown function", "region of unknown function",
    "repeat of unknown function", "module of unknown function",
    "duf", "upf", "uncharacteri", "putative uncharacteri",
    "conserved hypothetical", "hypothetical", "unknown protein",
    "unknown function", "orf", "pfam uncharacteri",
)
_UNKNOWN_FUNC_RE = re.compile(
    r'^\s*(?:(?:bacterial|viral|archaeal|eukaryotic|fungal|plant|marine|'
    r'transmembrane|integral membrane|membrane)\s+)?'
    r'(?:domain|protein|family|repeat|region|module)\s+of\s+unknown\s+function',
    re.IGNORECASE,
)
# Group 8 — DB-id-prefixed hypothetical ("FIG####: … hypothetical").
_DB_ID_HYPOTHETICAL_RE = re.compile(
    r'^(?:fig\d+|tigr\d+)[:\s].*(?:hypothetical|conserved hypothetical)', re.IGNORECASE,
)
# Group 4 — eggNOG "Belongs to the UPF#### family" (Uncharacterized Protein Family).
_UPF_ONLY_RE = re.compile(r'^\s*belongs to the upf\d+', re.IGNORECASE)
# Group 5 — a description that is nothing but a DUF tag ("Pfam:DUF955", "DUF955 family").
_BARE_DUF_RE = re.compile(r'^\s*(?:pfam:)?\(?duf\d+\)?(?:\s+(?:family|domain))?\s*$', re.IGNORECASE)
# Group 6 — "protein containing domains DUF###" (RAST multi-DUF stubs).
_PROTEIN_DOMAINS_DUF_RE = re.compile(r'^\s*protein containing domains?\s+duf', re.IGNORECASE)
_HYPOTHETICAL_ANYWHERE_RE = re.compile(r'\bhypothetical\b', re.IGNORECASE)
_PROTEIN_CONSERVED_IN_BACTERIA_RE = re.compile(r'\bprotein\s+conserved\s+in\s+bacteria\b', re.IGNORECASE)
_INTEGRAL_MEMBRANE_PROTEIN_RE = re.compile(r'^\s*integral\s+membrane\s+protein\s*$', re.IGNORECASE)
# Function-word guard for Groups 3 & 7: presence of a molecular/cellular ACTIVITY
# word means the annotation is informative despite an "unknown/uncharacterized"
# qualifier. Localization-only words are deliberately absent.
_FUNCTION_SIGNAL_RE = re.compile(
    r'(ase\b|transport|permease|pump|export|import|channel|carrier|symport|antiport|'
    r'bind|synth|kinas|reductas|hydrolas|transferas|isomeras|ligas|lyas|mutas|oxidas|'
    r'dehydrogen|regulat|repressor|activator|\bfactor\b|receptor|sensor|subunit|ribosom|'
    r'polymeras|oxidoreduct|enzyme|proteas|peptidas|nucleas|phosphatas|efflux|resistance|'
    r'virulence|toxin|flippase|homeostasis|assembly|motility|adhesin|chaperone|'
    r'helicas|topoisomeras|gyrase|recombinas|integras|transposas|methylas|glycosyl|'
    r'dismutas|catalas|peroxidas|cytochrome|ferredoxin|cytoskelet|cell division|'
    r'degradation|tolerance|translation|utilization)', re.IGNORECASE,
)
_UNCHARACTERIZED_RE = re.compile(r'\buncharacteri[sz]ed\b', re.IGNORECASE)
# Group 7 helper: what may follow "conserved protein" and still be uninformative
# -- a bare locus tag (e.g. "YqhG", "CreA") or a DUF/UPF tag. Anything more (a
# named domain/fold or a functional clause) keeps the hit.
_LOCUS_OR_TAG_RE = re.compile(r'^\(?(?:duf\d+|upf\d+|[a-z]{1,5}\d{0,4}[a-z]?\d{0,4})\)?$', re.IGNORECASE)


def is_uninformative(val: str) -> bool:
    """True iff ``val`` falls in the UNINFORMATIVE HIT CATEGORY above."""
    v = val.strip().lower()
    # Group 1
    if not v or v in _UNINFORMATIVE:
        return True
    # Group 2
    for prefix in _UNINFORMATIVE_PREFIXES:
        if v.startswith(prefix):
            return True
    # Groups 8, 4, 5, 6
    if (_UNKNOWN_FUNC_RE.match(v) or _DB_ID_HYPOTHETICAL_RE.match(v)
            or _UPF_ONLY_RE.match(v) or _BARE_DUF_RE.match(v)
            or _HYPOTHETICAL_ANYWHERE_RE.search(v)
            or _PROTEIN_CONSERVED_IN_BACTERIA_RE.search(v)
            or _INTEGRAL_MEMBRANE_PROTEIN_RE.match(v)
            or _PROTEIN_DOMAINS_DUF_RE.match(v)):
        return True
    # Groups 3 & 7 — an "unknown-function" / "uncharacterized" marker anywhere,
    # only when no molecular-function word rescues it.
    if not _FUNCTION_SIGNAL_RE.search(v):
        if "of unknown function" in v:
            return True
        if _UNCHARACTERIZED_RE.search(v):
            return True
        if v.startswith("conserved protein"):
            rest = v[len("conserved protein"):].strip(" ,.;:-")
            if not rest or _BARE_DUF_RE.match(rest) or _LOCUS_OR_TAG_RE.match(rest):
                return True
    return False


def compute_c1(labeled_row, cluster_row):
    informative_tools = []
    for tool in STANDALONE_DECISION_TOOLS:
        desc = labeled_row.get("RAST_description", "") if tool == "RAST" \
            else labeled_row.get(f"{tool}_best_hit_description", "")
        if desc and not is_uninformative(desc):
            informative_tools.append(tool)
    if cluster_row.get("tigrfam_cluster_status", "no_evidence") != "no_evidence":
        informative_tools.append(TIGRFAM_CLUSTER_SLOT)
    return len(informative_tools) / len(DECISION_TOOLS), informative_tools


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True, help="labeled-genes.tsv")
    parser.add_argument("--cluster-agreement-input", required=True,
                        help="add-cluster-agreement.py's output TSV (labeled-genes-cluster-agreement.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    cluster_path = Path(args.cluster_agreement_input)
    if not labeled_path.is_file():
        print(f"[score-c1-tool-coverage] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)
    if not cluster_path.is_file():
        print(f"[score-c1-tool-coverage] ERROR: input not found: {cluster_path}", file=sys.stderr)
        raise SystemExit(1)

    cluster_by_gene = {}
    with open(cluster_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                cluster_by_gene[fid] = row

    out_columns = _IDENTITY_COLUMNS + [
        "c1_score", "c1_informative_tool_count", "c1_total_tools_considered",
        "c1_informative_tools", "c1_formula",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    score_sum = 0.0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            cluster_row = cluster_by_gene.get(fid, {})
            c1, informative_tools = compute_c1(row, cluster_row)
            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c1_score"] = f"{c1:.4f}"
            out_row["c1_informative_tool_count"] = str(len(informative_tools))
            out_row["c1_total_tools_considered"] = str(len(DECISION_TOOLS))
            out_row["c1_informative_tools"] = ";".join(informative_tools)
            out_row["c1_formula"] = f"{len(informative_tools)}/{len(DECISION_TOOLS)} = {c1:.4f}"
            writer.writerow(out_row)
            score_sum += c1
            n += 1

    print(f"[score-c1-tool-coverage] Wrote {n} genes → {output_path}")
    print(f"    mean C1 = {score_sum / n:.4f}" if n else "    no genes scored")


if __name__ == "__main__":
    main()
