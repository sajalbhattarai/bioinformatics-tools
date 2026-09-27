"""
Counts a genome's report numbers from its FINAL confidence table: genes,
operons, support, confidence tiers and the Map's ring glyph.
Mirrors margie-fe's lib/workspace/final-summary.ts (FinalTally) so both sides
give the same numbers; counting here avoids sending 10-30 MB tables to the browser.
"""
from __future__ import annotations

import re
from typing import Iterable

GLYPH_TIERS = ['highest', 'high', 'medium', 'fair', 'low', 'none']
GLYPH_BINS = 240

_PLAIN = re.compile(r'^Column-[A-Z]+:\s*')
_GENE_CONTIG = re.compile(r'^(.*)_\d+[+-]\d+$')
_FEATURE_CONTIG = re.compile(r'^(.*)_\d+$')
_YES = re.compile(r'^yes$', re.I)
_NO_HITS = re.compile(r'^no db hits$', re.I)


def plain_name(c: str) -> str:
    """Strips the "Column-A: " prefix from a FINAL table column name."""
    return _PLAIN.sub('', c)


def _number(text: str) -> float | None:
    """Parses a cell like JavaScript's Number(): '' is 0, unparsable text is None."""
    text = text.strip()
    if not text:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return None


def tally(lines: Iterable[str]) -> dict:
    """Counts genes, operons, support and tiers from FINAL table lines and builds the ring glyph."""
    columns: list[str] = []
    idx: dict[str, int] = {}
    genes = in_operons = with_support = 0
    envelope = ''
    operons: set[str] = set()
    tiers: dict[str, int] = {}
    placed: list[tuple[str, float, float, str, int]] = []

    for raw in lines:
        line = raw.rstrip('\r\n')
        if not columns:
            columns = line.split('\t')
            idx = {plain_name(c): i for i, c in enumerate(columns)}
            continue
        if not line:
            continue
        cells = line.split('\t')

        def get(name: str) -> str:
            i = idx.get(name)
            return cells[i].strip() if i is not None and i < len(cells) else ''

        genes += 1
        if _YES.match(get('IS_IN_OPERON?')):
            in_operons += 1
        op = get('UniOP_OPERON_id')
        if op.startswith('operon'):
            operons.add(op)
        desc = get('best_consensus_product_descriptor')
        if desc and not _NO_HITS.match(desc):
            with_support += 1
        tier = get('CONFIDENCE_TIER')
        if tier:
            tiers[tier] = tiers.get(tier, 0) + 1
        envelope = envelope or get('ENVELOPE') or get('envelope')

        start = _number(get('RAST_start') or get('gene_start'))
        end = _number(get('RAST_end') or get('gene_end'))
        if start is not None and end is not None and start > 0 and end > 0:
            m = _GENE_CONTIG.match(get('gene_id'))
            if m:
                contig = m.group(1)
            else:
                m = _FEATURE_CONTIG.match(get('RAST_feature_id'))
                contig = m.group(1) if m else ''
            t = (get('CONFIDENCE_TIER_hybrid') or tier).lower()
            ti = GLYPH_TIERS.index(t) if t in GLYPH_TIERS else -1
            placed.append((contig, min(start, end), max(start, end), get('RAST_strand'), 5 if ti < 0 else ti))

    return {
        'columns': [plain_name(c) for c in columns],
        'genes': genes,
        'operons': len(operons),
        'inOperons': in_operons,
        'withSupport': with_support,
        'envelope': envelope,
        'tiers': tiers,
        'glyph': _glyph(placed),
    }


def _glyph(genes: list[tuple[str, float, float, str, int]]) -> dict | None:
    """Bins a genome's genes around one ring, contigs end to end in the order they first appear."""
    if not genes:
        return None
    lengths: dict[str, float] = {}
    for contig, _, end, _, _ in genes:
        lengths[contig] = max(lengths.get(contig, 0), end)
    offset: dict[str, float] = {}
    total = 0.0
    for contig, length in lengths.items():
        offset[contig] = total
        total += length
    counts = {side: [[0] * len(GLYPH_TIERS) for _ in range(GLYPH_BINS)] for side in ('plus', 'minus')}
    for contig, start, end, strand, tier in genes:
        at = offset.get(contig, 0) + (start + end) / 2
        b = min(GLYPH_BINS - 1, int((at / total) * GLYPH_BINS)) if total else 0
        counts['minus' if strand == '-' else 'plus'][b][tier] += 1

    def dominant(c: list[int]) -> int:
        return -1 if not any(c) else c.index(max(c))

    return {
        'length': int(total) if float(total).is_integer() else total,
        'contigs': len(lengths),
        'plus': [dominant(c) for c in counts['plus']],
        'minus': [dominant(c) for c in counts['minus']],
    }
