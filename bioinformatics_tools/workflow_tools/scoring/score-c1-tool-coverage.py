#!/usr/bin/env python3
"""score-c1-tool-coverage.py — scoring stage, metric C1: tool coverage.

C1 = informative sources / 7, over RAST, KEGG, EGGNOG, COG, PFAM, UNIPROT and
one TIGRFAM cluster slot (PGAP/TIGRFAM/NCBIfam share models). A hit counts only
if its description is informative. HAMAP, PIRSF, CDD and GENEPROP overlap the
counted tools, and MEROPS/TCDB/DBCAN cover few genes, so none are counted.
Writes the score, count, contributing tools and a "k/7 = x" formula per gene.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

# Independent decision tools; TIGRFAM/PGAP/NCBIfam count as one slot.
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

# ─── Uninformative hit category ──────────────────────────────────────────────
# Descriptions matching these groups name no function and count 0 toward C1.
# Kept identical to labeling/assign-canonical-label.py. Only an activity word
# (…ase, transport, regulator, …) rescues an unknown-function phrase; location
# words do not.

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
# eggNOG "Psort location ..." descriptions give localisation only, not a product.
_PSORT_ONLY_RE = re.compile(r'^\s*psort location\b', re.IGNORECASE)
# Group 5 — a description that is nothing but a DUF tag ("Pfam:DUF955", "DUF955 family").
_BARE_DUF_RE = re.compile(r'^\s*(?:pfam:)?\(?duf\d+\)?(?:\s+(?:family|domain))?\s*$', re.IGNORECASE)
# Group 6 — "protein containing domains DUF###" (RAST multi-DUF stubs).
_PROTEIN_DOMAINS_DUF_RE = re.compile(r'^\s*protein containing domains?\s+duf', re.IGNORECASE)
_HYPOTHETICAL_ANYWHERE_RE = re.compile(r'\bhypothetical\b', re.IGNORECASE)
_PROTEIN_CONSERVED_IN_BACTERIA_RE = re.compile(r'\bprotein\s+conserved\s+in\s+bacteria\b', re.IGNORECASE)
_INTEGRAL_MEMBRANE_PROTEIN_RE = re.compile(r'^\s*integral\s+membrane\s+protein\s*$', re.IGNORECASE)
# Function-word guard for Groups 3 & 7 (activity words only, no location words).
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
# Group 7 helper: a bare locus tag or DUF/UPF tag after "conserved protein".
_LOCUS_OR_TAG_RE = re.compile(r'^\(?(?:duf\d+|upf\d+|[a-z]{1,5}\d{0,4}[a-z]?\d{0,4})\)?$', re.IGNORECASE)


def is_uninformative(val: str) -> bool:
    """Returns True when the description matches the uninformative groups above (regex and word lists)."""
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
            or _UPF_ONLY_RE.match(v) or _PSORT_ONLY_RE.match(v) or _BARE_DUF_RE.match(v)
            or _HYPOTHETICAL_ANYWHERE_RE.search(v)
            or _PROTEIN_CONSERVED_IN_BACTERIA_RE.search(v)
            or _INTEGRAL_MEMBRANE_PROTEIN_RE.match(v)
            or _PROTEIN_DOMAINS_DUF_RE.match(v)):
        return True
    # Groups 3 & 7 — unknown-function markers without a function word.
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
    """Returns (C1 score, contributing tools) from best-hit descriptions and the TIGRFAM cluster status."""
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
    """Joins labeled genes to cluster agreement with csv and writes the C1 table."""
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
