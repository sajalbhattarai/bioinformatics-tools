"""Genome chat -- evidence-grounded Q&A over a finished run.

Separate from phase15 (workflow_tools/llm/score-genes-llm.py), which *produces*
scores as part of a run using the LoRA-fused model. This only *explains*
results that already exist, and it runs the BASE model: the fuse was trained to
emit a score, and the behaviour needed here is the opposite -- say "the evidence
does not answer that" instead of producing a confident number.

Route:  POST /v1/llm/chat  -> {"answer": str, "context_summary": {...}}
        GET  /v1/llm/status -> whether the inference endpoint is up

Grounding: the browser sends only an IDENTIFIER (job + organism + optional
gene/operon). This module rebuilds the context server-side from that job's own
FINAL table, so a tampered page cannot widen its own grounding -- constraining
the context is the entire point of the feature.

The model is reached via chat_server.py running on a GPU node, discovered
through the advert file it publishes. No advert = chat is offline, which is
reported as a clean 503 rather than an error.
"""
from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import urllib.error
import urllib.request
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel

from bioinformatics_tools.api.auth import get_current_user
from bioinformatics_tools.api.routers.ssh import _build_connection, _resolve_job_work_dir
from bioinformatics_tools.utilities import ssh_sftp

LOGGER = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/llm", tags=["llm"])

ADVERT_PATH = Path(os.getenv(
    "MARGIE_CHAT_ADVERTISE",
    os.path.expanduser("~/.local/share/bsp/chat-endpoint.json"),
))

FINAL_TSV = "FINAL_ANNOTATION_WITH_CONFIDENCE.tsv"
TIERS = ["highest", "high", "medium", "fair", "low"]

# Carried over from score-genes-llm.py's SYSTEM_PROMPT where it still applies.
# Those scoping rules encode real failure modes this pipeline already hit --
# neighbour evidence being attributed to the candidate gene, and EC granularity
# being read as conflict -- so they are repeated verbatim in intent rather than
# reinvented. What is new is the grounding contract and the citation duty.
SYSTEM_PROMPT = """\
You are explaining the results of a MARGIE(SB) genome annotation to the person \
who ran it. You are given that pipeline's own records. Explain what they mean. \
Do not re-score anything.

GROUNDING CONTRACT -- this overrides everything else:
  * Every claim must be traceable to a specific field or row in the records below.
  * You have no outside knowledge. You do not know what any gene "usually" does, \
what any organism is "known for", or what any database "normally" reports. If it \
is not in the records, you do not know it.
  * If the records do not answer the question, say what is missing and stop. Do \
not substitute plausible biology.
  * Absent is not zero, not average, and not "typical". Never infer a missing value.
  * Cite the source of each claim inline -- the column name or the tool name.

SCOPING RULES (these are real error modes, not hypotheticals):
  * The gene asked about is the CANDIDATE. Evidence listed for other genes or \
operon members belongs to THEM. Never attribute a neighbour's EC number, hit or \
description to the candidate.
  * EC specificity differences are NOT conflicts. 2.7.1.2 vs 2.7.1.- is the same \
enzyme at different resolution. Only different reaction classes (2.7.1.2 vs \
1.1.1.27) are genuine conflicts.
  * C1 (tool coverage) and C2 (operon reliability) are arithmetic. Report them as \
given; do not re-derive them.
  * "Different function" is not "contradiction". Uninformative or hypothetical \
hits carry no signal either way -- say so rather than reading them as negative.

VOCABULARY:
  * Tiers, most to least confident: highest, high, medium, fair, low. \
NOT_APPLICABLE_NON_CODING marks non-coding features.
  * Two tier columns exist -- CONFIDENCE_TIER and CONFIDENCE_TIER_hybrid. State \
which one you are quoting.
  * Confidence values are decimals 0-1, never percentages.
  * NEEDS_REVIEW means the pipeline wants human attention. It does NOT mean the \
annotation is wrong.

WHAT YOU ARE LOOKING AT:
  You are given (a) the COLUMN SCHEMA of this run's FINAL_ANNOTATION_WITH_
CONFIDENCE.tsv -- what each field means -- (b) genome-level aggregates, and \
(c) the FULL RECORDS of any genes relevant to the question, when the question \
names one. Study the records before answering. Aggregates describe the WHOLE \
GENOME and say nothing about any individual gene: never read a genome total as \
a statement about one gene. If a gene's record is not present, you do not know \
that gene's tier, score or review status -- say so and ask for the gene name.

ANSWER FORMAT -- every answer uses these headings, and omits any that do not apply:

  **Answer** — the direct answer, first, in one or two sentences.

  **From the evidence** — each supporting fact with the field or tool it came \
from, e.g. "CONFIDENCE_TIER = high", "NEEDS_REVIEW? = yes". Only things \
literally present in the records.

  **Interpretation** — anything you inferred, reasoned about, or drew on \
general biological knowledge to say. State plainly that it is interpretation, \
not a pipeline output. If you used no interpretation at all, omit this heading \
entirely rather than writing "none".

  **Not in the evidence** — what the records cannot answer, named specifically.

This separation is the point: the reader must be able to tell, at a glance, \
which sentences are the pipeline's findings and which are your reading of them. \
Never blend the two.

LENGTH: as long as the question genuinely needs, up to the token budget. A \
one-line question gets a few lines. Do not pad, do not restate the question, \
do not list facts that were not asked for."""


class ChatRequest(BaseModel):
    job_id: str
    organism: str                 # organism folder name inside the run
    question: str
    gene_id: str | None = None    # optional: narrows context to one gene
    max_tokens: int = 5000


def _endpoint() -> dict:
    """Current inference endpoint, or 503 if chat is not running."""
    if not ADVERT_PATH.is_file():
        raise HTTPException(
            status_code=503,
            detail="Genome chat is offline — no inference server is running. "
                   "Start it with: sbatch chat-server.sbatch",
        )
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        return {"url": f"http://{adv['host']}:{adv['port']}", "model": adv.get("model", "")}
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"Chat endpoint unreadable: {exc}")


def _read_final(job_id: str, organism: str, current_user: dict) -> list[dict]:
    """Read this organism's FINAL table off the user's cluster."""
    conn = _build_connection(current_user)
    work_dir = _resolve_job_work_dir(job_id, current_user, conn)
    # Post-reorganize the FINAL table sits at the organism top level; mid-run it
    # is still under scoring/. Try both rather than assuming a layout.
    for rel in (f"{organism}/{FINAL_TSV}", f"{organism}/scoring/{FINAL_TSV}"):
        try:
            raw = b"".join(ssh_sftp.stream_remote_file(f"{work_dir}/{rel}", connection=conn))
            return list(csv.DictReader(io.StringIO(raw.decode("utf-8", "replace")),
                                       delimiter="\t"))
        except FileNotFoundError:
            continue
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Could not read {rel}: {exc}")
    raise HTTPException(status_code=404,
                        detail=f"No {FINAL_TSV} for organism '{organism}' in this job.")


def _col(row: dict, name: str) -> str:
    """FINAL tables carry prefixed headers ('Column-AS: CONFIDENCE_TIER')."""
    if name in row:
        return (row[name] or "").strip()
    for k, v in row.items():
        if k and k.split(":")[-1].strip() == name:
            return (v or "").strip()
    return ""


# Fields worth naming explicitly. The model cannot reason about a table whose
# columns it has never been told the meaning of -- without this it treats every
# header as opaque and falls back on genome totals.
_SCHEMA_HELP = {
    "gene_id": "gene identifier, <accession>_<start><strand><len>",
    "organism_name": "genome this gene belongs to",
    "canonical_label": "the annotation the pipeline settled on",
    "RAST_start": "start coordinate on its replicon",
    "RAST_end": "end coordinate",
    "RAST_strand": "+ or -",
    "CONFIDENCE_TIER": "highest | high | medium | fair | low | NOT_APPLICABLE_NON_CODING",
    "CONFIDENCE_TIER_hybrid": "same scale, operon-context-adjusted variant",
    "NEEDS_REVIEW?": "yes/no — pipeline wants a human to look",
    "IS_IN_OPERON?": "yes/no",
    "UniOP_OPERON_id": "operon identifier, or a NOT_* sentinel",
    "C1_score_tool_coverage": "fraction of tools returning an informative hit",
    "C2_score_operon_probability": "geometric mean of operon-member probabilities",
    "C3_score_operon_context": "operon/neighbourhood coherence",
    "C4_score_ec_agreement": "agreement among EC numbers",
    "confidence_score": "final combined confidence, 0-1",
}


def _schema_block(rows: list[dict]) -> str:
    """Tell the model what this run's FINAL table actually contains."""
    if not rows:
        return ""
    names = []
    for k in rows[0].keys():
        if not k:
            continue
        bare = k.split(":")[-1].strip()
        names.append(f"  {bare}" + (f" — {_SCHEMA_HELP[bare]}" if bare in _SCHEMA_HELP else ""))
    return ("FINAL_ANNOTATION_WITH_CONFIDENCE.tsv — one row per gene, "
            f"{len(rows)} rows, {len(names)} columns:\n" + "\n".join(names))


_STOPWORDS = {
    "what", "which", "does", "have", "this", "that", "gene", "genes", "genome",
    "about", "there", "any", "can", "you", "tell", "give", "more", "into",
    "details", "detail", "look", "the", "and", "for", "with", "from", "its",
    "it", "is", "are", "has", "how", "many", "much", "please", "called",
}


def _search_genes(rows: list[dict], question: str, limit: int = 6) -> list[dict]:
    """Rows whose label/id matches meaningful words in the question.

    Without this the model is handed genome totals and nothing else, and when
    asked about a named gene it invents a per-gene answer from an aggregate --
    the exact failure this is here to remove.
    """
    words = {w for w in re.findall(r"[A-Za-z0-9_.\-]{3,}", question.lower())
             if w not in _STOPWORDS}
    if not words:
        return []
    scored = []
    for r in rows:
        hay = " ".join(filter(None, (
            _col(r, "gene_id"), _col(r, "canonical_label"),
            _col(r, "UniOP_OPERON_id"),
        ))).lower()
        if not hay:
            continue
        hits = sum(1 for w in words if w in hay)
        if hits:
            scored.append((hits, r))
    scored.sort(key=lambda x: -x[0])
    return [r for _, r in scored[:limit]]


def _render_gene(r: dict) -> str:
    out = []
    for k, v in r.items():
        if not k:
            continue
        v = (v or "").strip()
        if not v:
            continue
        bare = k.split(":")[-1].strip()
        if bare == "RAST_na_sequence":
            v = f"<{len(v)} nt, omitted>"
        out.append(f"    {bare}: {v}")
    return "\n".join(out)


def _genome_context(rows: list[dict], organism: str) -> tuple[str, dict]:
    """Whole-genome scope: aggregates only. Per-gene records cannot fit."""
    tally = {t: 0 for t in TIERS}
    noncoding = flagged = operonic = 0
    for r in rows:
        t = _col(r, "CONFIDENCE_TIER")
        if t in tally:
            tally[t] += 1
        elif t.startswith("NOT_"):
            noncoding += 1
        if _col(r, "NEEDS_REVIEW?").lower() == "yes":
            flagged += 1
        if _col(r, "IS_IN_OPERON?").lower() == "yes":
            operonic += 1
    lines = [f"GENOME: {organism}", f"total genes: {len(rows)}",
             "confidence tier counts (CONFIDENCE_TIER):"]
    lines += [f"  {t}: {tally[t]}" for t in TIERS]
    lines += [f"  NOT_APPLICABLE_NON_CODING: {noncoding}",
              f"genes flagged NEEDS_REVIEW: {flagged}",
              f"genes in an operon (IS_IN_OPERON?): {operonic}"]
    summary = {"scope": "genome", "genes": len(rows), "flagged": flagged,
               "tiers": tally, "non_coding": noncoding}
    return "\n".join(lines), summary


def _gene_context(rows: list[dict], gene_id: str, organism: str) -> tuple[str, dict]:
    """Single-gene scope: that gene's full FINAL record, verbatim."""
    match = next((r for r in rows if _col(r, "gene_id") == gene_id), None)
    if match is None:
        raise HTTPException(status_code=404, detail=f"gene_id '{gene_id}' not in this organism.")
    lines = [f"GENOME: {organism}", f"GENE: {gene_id}", "",
             "This gene's complete record from FINAL_ANNOTATION_WITH_CONFIDENCE.tsv.",
             "Field names are given exactly as the pipeline wrote them:", ""]
    for k, v in match.items():
        if k is None:
            continue
        v = (v or "").strip()
        if not v:
            continue
        # The nucleotide sequence is long and never the subject of a question.
        if k.split(":")[-1].strip() == "RAST_na_sequence":
            v = f"<{len(v)} nt, omitted>"
        lines.append(f"  {k.split(':')[-1].strip()}: {v}")
    return "\n".join(lines), {"scope": "gene", "gene_id": gene_id,
                              "fields": sum(1 for v in match.values() if v)}


def _resolve_model_path(current_user: dict) -> str:
    """<db_root>/llm/base, where db_root is the database folder the user set in
    the GUI. Never hardcoded: the whole point of the Profile setting is that the
    database lives wherever the user put it.

    The BASE model deliberately, not fused-model — see this module's docstring.
    """
    conn = _build_connection(current_user)
    home = current_user["home_dir"]
    cfg = {}
    try:
        cfg = ssh_sftp.read_remote_yaml(
            f"{home}/.config/bioinformatics-tools/config.yaml", connection=conn) or {}
    except Exception as exc:
        LOGGER.warning("Could not read config for db_root: %s", exc)

    # Same precedence the workflow uses: an explicit db.llm wins, else
    # <workflow>.db_root/llm. db_path() itself is Snakemake-side, so the two
    # keys are read directly here rather than importing the workflow config.
    explicit = (cfg.get("db") or {}).get("llm")
    root = None
    if explicit:
        root = str(explicit)
    else:
        for wf in ("margie_sb", "margie"):
            r = (cfg.get(wf) or {}).get("db_root")
            if r:
                root = f"{str(r).rstrip('/')}/llm"
                break
    if not root:
        root = (current_user.get("db_root") or "").rstrip("/")
        root = f"{root}/llm" if root else ""
    if not root:
        raise HTTPException(
            status_code=400,
            detail="No database folder configured. Set it in Profile → database "
                   "path (the folder containing llm/base) before starting chat.",
        )
    return f"{root.rstrip('/')}/base"


@router.post("/start")
def start_chat(current_user: dict = Depends(get_current_user)):
    """Submit the GPU job that hosts the chat model — the 'Start chat' button.

    Chat needs a resident model, and a SLURM allocation is the only way to hold
    one on this cluster. Returns immediately with the job id; the model takes a
    minute or two to load, during which /status reports offline.
    """
    if ADVERT_PATH.is_file():
        st = chat_status(current_user)
        if st.get("online"):
            return {"started": False, "already_running": True, **st}

    model = _resolve_model_path(current_user)
    conn = _build_connection(current_user)
    ssh = conn.connect()
    try:
        _, out, err = ssh.exec_command(
            f"test -d {model!r} && echo OK || echo MISSING")
        if out.read().decode().strip() != "OK":
            raise HTTPException(
                status_code=400,
                detail=f"Base model not found at {model}. Check the database "
                       "folder set in Profile.",
            )
        # Backend location is wherever this package lives on the cluster.
        backend = str(Path(__file__).resolve().parents[3])
        sbatch = f"{backend}/bioinformatics_tools/workflow_tools/llm/chat-server.sbatch"
        cmd = (f"MARGIE_LLM_MODEL={model!r} MARGIE_BACKEND={backend!r} "
               f"sbatch --parsable {sbatch!r}")
        _, out, err = ssh.exec_command(cmd)
        job = out.read().decode().strip()
        stderr = err.read().decode().strip()
        if not job.isdigit():
            raise HTTPException(status_code=500,
                                detail=f"sbatch failed: {stderr or job or 'no job id'}")
    finally:
        ssh.close()

    LOGGER.info("chat server submitted as job %s (model %s)", job, model)
    return {"started": True, "job_id": job, "model": model,
            "detail": "Chat server queued. It reports online once the model has loaded "
                      "(usually 1–3 minutes, longer if the GPU queue is busy)."}


@router.post("/stop")
def stop_chat(current_user: dict = Depends(get_current_user)):
    """End interactive mode — called when the map page is closed.

    Asks the server to exit (which frees the GPU and removes its advert). This
    is the fast path only: it never arrives if the browser crashed, the laptop
    slept, or the network dropped, so chat_server also exits on its own idle
    timeout. Never rely on this alone to release a GPU.
    """
    if not ADVERT_PATH.is_file():
        return {"stopped": False, "detail": "chat was not running"}
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        req = urllib.request.Request(
            f"http://{adv['host']}:{adv['port']}/shutdown", data=b"{}",
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=10):
            pass
        LOGGER.info("chat server on %s asked to stop", adv.get("host"))
        return {"stopped": True}
    except Exception as exc:
        # Unreachable means it is already gone; clear the stale advert so the
        # UI stops pointing at a dead endpoint.
        try:
            ADVERT_PATH.unlink()
        except Exception:
            pass
        return {"stopped": True, "detail": f"endpoint already gone ({exc}); advert cleared"}


@router.get("/status")
def chat_status(current_user: dict = Depends(get_current_user)):
    """Whether the chat backend is reachable — lets the UI show an honest state."""
    if not ADVERT_PATH.is_file():
        return {"online": False, "detail": "no inference server running"}
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        url = f"http://{adv['host']}:{adv['port']}/health"
        with urllib.request.urlopen(url, timeout=5) as r:
            return {"online": True, "model": json.loads(r.read()).get("model", ""),
                    "host": adv["host"]}
    except Exception as exc:
        # Advert present but unreachable => the job died without cleaning up.
        return {"online": False, "detail": f"endpoint not responding: {exc}"}


@router.post("/chat")
def chat(body: ChatRequest, current_user: dict = Depends(get_current_user)):
    question = (body.question or "").strip()
    if not question:
        raise HTTPException(status_code=400, detail="question is required")

    ep = _endpoint()
    rows = _read_final(body.job_id, body.organism, current_user)

    if body.gene_id:
        context, summary = _gene_context(rows, body.gene_id, body.organism)
        matched = []
    else:
        context, summary = _genome_context(rows, body.organism)
        # Pull in the full record of any gene the question names. Aggregates
        # alone caused the model to answer per-gene questions from genome
        # totals ("the gene is not flagged" from a count of 1212 flagged).
        matched = _search_genes(rows, question)
        if matched:
            context += ("\n\nGENES MATCHING THIS QUESTION — full records, "
                        "one block each:\n")
            for r in matched:
                context += (f"\n  [{_col(r, 'gene_id') or '?'}] "
                            f"{_col(r, 'canonical_label')}\n{_render_gene(r)}\n")
            summary["matched_genes"] = [_col(r, "gene_id") for r in matched]
        else:
            context += ("\n\nNo gene in this genome matched the wording of the "
                        "question, so no per-gene record is included. Only the "
                        "genome-level aggregates above are available.")

    context = _schema_block(rows) + "\n\n" + context

    prompt = (f"{context}\n\n"
              f"----\nQUESTION: {question}\n\n"
              "Study the records above, then answer using the required headings. "
              "Cite the field or tool behind every fact, and keep anything you "
              "inferred under Interpretation so it cannot be mistaken for a "
              "pipeline output.")

    try:
        req = urllib.request.Request(
            f"{ep['url']}/chat",
            data=json.dumps({"system": SYSTEM_PROMPT, "prompt": prompt,
                             "max_tokens": body.max_tokens}).encode(),
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=180) as r:
            answer = json.loads(r.read()).get("text", "")
    except urllib.error.URLError as exc:
        raise HTTPException(status_code=503,
                            detail=f"Chat backend unreachable: {exc}. It may have hit walltime.")
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Chat failed: {exc}")

    return {"answer": answer, "context_summary": summary, "model": ep["model"]}
