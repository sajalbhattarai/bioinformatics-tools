#!/usr/bin/env python3
"""make-final-excel.py — convert the FINAL annotation TSV to a coloured Excel workbook.

Each data row is filled edge to edge with its CONFIDENCE_TIER_HYBRID tier colour
(highest=blue, high=green, medium=yellow, fair=orange, low=red; non-coding=grey),
so every row reads as a single tier band. Rows flagged NEEDS_REVIEW? = yes also
get a box border around the whole row. A second "Legend" sheet explains both.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

csv.field_size_limit(10_000_000)

# ── FINAL-file coloring — MATCHES the operon-diagram figures
#    (reportfig_lib.CONF_TIER_COLOR). Every row is tinted edge to edge with its
#    CONFIDENCE_TIER_HYBRID tier colour; review rows also get a box border.
#    Kept in sync with api/routers/ssh.py. ─────────────────────────────────────
_TIER_BRIGHT = {
    "highest": "1F77FF",   # blue
    "high":    "00B84D",   # green
    "medium":  "FFCC00",   # yellow
    "fair":    "FF8C00",   # orange
    "low":     "EE2233",   # red
}
_ROW_NONCODING_BG = "F2F2F2"   # non-coding / unscored rows (no confidence tier)
_ROW_NONCODING_FG = "8A8A8A"
_ROW_TINT = 0.72               # how light the tier colour is across each row
_REVIEW_SIDE = Side(style="medium", color="000000")  # box border on review rows


def _tint_hex(h: str, toward_white: float = 0.86) -> str:
    h = h.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r = round(r + (255 - r) * toward_white)
    g = round(g + (255 - g) * toward_white)
    b = round(b + (255 - b) * toward_white)
    return f"{r:02X}{g:02X}{b:02X}"

# ── header ─────────────────────────────────────────────────────────────────
_HDR_BG = "203864"
_HDR_FG = "FFFFFF"

# ── column widths (characters) ──────────────────────────────────────────────
_COL_WIDTHS: dict[str, int] = {
    "organism_name":   35,
    "feature_id":      30,
    "gene_id":         28,
    "na_seq":           8,
    "aa_seq":           8,
    "na_length":       10,
    "aa_length":       10,
    "gene_start":      12,
    "gene_end":        12,
    "RAST_feature_type": 16,
    "RAST_strand":     10,
    "operon_id":       16,
    "operon_member_count": 14,
    "operon_gene_position_in_operon": 20,
    "gene_fingerprint_with_hash": 40,
    "best_consensus_product_descriptor": 36,
    "product_descriptor_source":    20,
    "product_descriptor_source_id": 24,
    "product_descriptor_hierarchy": 22,
    "product_descriptor_audit_trail": 36,
    "product_descriptor_confirmatory_summary_specialized_tools": 38,
    "product_descriptor_hierarchy_tier_name": 32,
    "c1_score_database_coverage": 14,
    "c1_score_reasoning": 40,
    "c2_score_operon_probability": 14,
    "c2_score_reasoning": 40,
    "c3_score_operon_context": 14,
    "c3_reasoning": 40,
    "c4_score_EC_agreement": 14,
    "c4_reasoning": 40,
    "c4_ec_agreement_status": 20,
    "preliminary_confidence_c1_c4": 16,
    "final_confidence_operon_context": 18,
    "does_context_improve_confidence?": 20,
    "confidence_tier": 14,
    "needs_review": 12,
    "needs_review_reason": 44,
    "gene_id_copy": 28,
    "gene_start_copy": 12,
    "gene_end_copy": 12,
    "best_consensus_product_descriptor_copy": 36,
}
_DEFAULT_WIDTH = 16


_FILL_CACHE: dict[str, PatternFill] = {}
_FONT_CACHE: dict[tuple[str, bool], Font] = {}


def _fill(hex_color: str) -> PatternFill:
    # reuse one PatternFill object per colour -- openpyxl then dedups by identity,
    # which turns ~200k per-cell style assignments from minutes into seconds.
    f = _FILL_CACHE.get(hex_color)
    if f is None:
        f = PatternFill("solid", fgColor=hex_color)
        _FILL_CACHE[hex_color] = f
    return f


def _font(hex_color: str, bold: bool = False) -> Font:
    key = (hex_color, bold)
    f = _FONT_CACHE.get(key)
    if f is None:
        f = Font(color=hex_color, bold=bold)
        _FONT_CACHE[key] = f
    return f


def _norm_col(name: str) -> str:
    # FINAL_ANNOTATION_WITH_CONFIDENCE supports prefixed headers like:
    #   "[AN]-NEEDS_REVIEW?"  (legacy)
    #   "Column-AN: NEEDS_REVIEW?"  (current)
    return re.sub(r"^(?:\[[A-Z]+\]-|Column-[A-Z]+:\s*)", "", str(name or "").strip(), flags=re.IGNORECASE).strip().lower()


def _row_get(row: dict, *candidate_names: str) -> str:
    normalized_targets = {_norm_col(n) for n in candidate_names}
    for k, v in row.items():
        if _norm_col(k) in normalized_targets:
            return v
    return ""


def _row_tint(row: dict) -> tuple[str, str]:
    """(bg, fg) whole-row colour = the row's CONFIDENCE_TIER_HYBRID tier colour,
    tinted and applied across the ENTIRE row so each row reads as one tier band.
    Rows with no scored tier -- empty or NOT_APPLICABLE_NON_CODING (rna /
    prophage) -- get neutral grey. This is the only colouring; no per-cell
    accents."""
    tier = str(_row_get(row, "confidence_tier_hybrid", "CONFIDENCE_TIER_hybrid")).strip().lower()
    if tier not in _TIER_BRIGHT:
        return _ROW_NONCODING_BG, _ROW_NONCODING_FG
    return _tint_hex(_TIER_BRIGHT[tier], _ROW_TINT), "000000"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True,
                        help="FINAL-scored-labeled-genes-annotated.tsv")
    parser.add_argument("--output", required=True,
                        help="Output .xlsx path")
    args = parser.parse_args()

    in_path  = Path(args.input)
    out_path = Path(args.output)

    if not in_path.is_file():
        print(f"[make-final-excel] ERROR: input not found: {in_path}", file=sys.stderr)
        raise SystemExit(1)

    wb = Workbook()
    ws = wb.active
    ws.title = "Annotation Results"

    with open(in_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        headers = list(reader.fieldnames or [])

        # ── header row ──────────────────────────────────────────────────────
        hdr_fill = _fill(_HDR_BG)
        hdr_font = _font(_HDR_FG, bold=True)
        for ci, h in enumerate(headers, 1):
            cell = ws.cell(row=1, column=ci, value=h)
            cell.fill = hdr_fill
            cell.font = hdr_font
            cell.alignment = Alignment(horizontal="center", vertical="center",
                                       wrap_text=True)

        # ── data rows ───────────────────────────────────────────────────────
        # ONE shared alignment object for every data cell (creating a fresh
        # Alignment per cell is what made a 4.6k-row sheet take minutes).
        data_align = Alignment(wrap_text=False, vertical="top")
        ncol = len(headers)
        n = 0
        review_n = 0
        band_counts: dict[str, int] = {}
        for ri, row in enumerate(reader, 2):
            row_bg, row_fg = _row_tint(row)
            row_fill = _fill(row_bg)
            row_font = _font(row_fg)
            review = str(_row_get(row, "needs_review?", "NEEDS_REVIEW?",
                                  "needs_review")).strip().lower() == "yes"

            for ci, h in enumerate(headers, 1):
                cell = ws.cell(row=ri, column=ci, value=row.get(h, ""))
                cell.alignment = data_align
                cell.fill = row_fill
                cell.font = row_font
                if review:                          # box border around the whole review row
                    cell.border = Border(
                        top=_REVIEW_SIDE, bottom=_REVIEW_SIDE,
                        left=_REVIEW_SIDE if ci == 1 else None,
                        right=_REVIEW_SIDE if ci == ncol else None)

            # tally by the row's CONFIDENCE_TIER_HYBRID band (the whole-row colour)
            band = str(_row_get(row, "confidence_tier_hybrid",
                                "CONFIDENCE_TIER_hybrid")).strip().lower()
            if band not in _TIER_BRIGHT:
                band = "non-coding"
            band_counts[band] = band_counts.get(band, 0) + 1
            review_n += 1 if review else 0
            n += 1

    # ── freeze panes and row height ─────────────────────────────────────────
    ws.freeze_panes = "A2"
    ws.row_dimensions[1].height = 40

    # ── column widths ───────────────────────────────────────────────────────
    for ci, h in enumerate(headers, 1):
        ws.column_dimensions[get_column_letter(ci)].width = _COL_WIDTHS.get(h, _DEFAULT_WIDTH)

    # ── Legend sheet: what each row colour and the box border mean ───────────
    lg = wb.create_sheet("Legend")
    lg.column_dimensions["A"].width = 12
    lg.column_dimensions["B"].width = 14
    lg.column_dimensions["C"].width = 66

    def _legend_row(r: int, swatch_bg: str, label: str, meaning: str, boxed: bool = False) -> None:
        sw = lg.cell(r, 1)
        sw.fill = _fill(swatch_bg)
        if boxed:
            sw.border = Border(top=_REVIEW_SIDE, bottom=_REVIEW_SIDE,
                               left=_REVIEW_SIDE, right=_REVIEW_SIDE)
        lg.cell(r, 2, label).font = _font("000000", bold=True)
        lg.cell(r, 3, meaning)

    t = lg.cell(1, 1, "How to read this workbook")
    t.font = _font(_HDR_FG, bold=True)
    for c in (1, 2, 3):
        lg.cell(1, c).fill = _fill(_HDR_BG)

    lg.cell(3, 1, "ROW COLOUR = confidence tier (CONFIDENCE_TIER_HYBRID)").font = _font("000000", bold=True)
    tiers = [
        ("highest", "highest confidence  (final > 0.90)"),
        ("high",    "high confidence  (0.70 – 0.90)"),
        ("medium",  "medium confidence  (0.50 – 0.70)"),
        ("fair",    "fair confidence  (0.30 – 0.50)"),
        ("low",     "low confidence  (final <= 0.30)"),
    ]
    r = 4
    for tier, meaning in tiers:
        _legend_row(r, _tint_hex(_TIER_BRIGHT[tier], _ROW_TINT), tier, meaning)
        r += 1
    _legend_row(r, _ROW_NONCODING_BG, "non-coding", "no confidence score (rna / prophage / unscored)")

    r += 2
    lg.cell(r, 1, "ROW BORDER = review status").font = _font("000000", bold=True)
    _legend_row(r + 1, "FFFFFF", "boxed row", "flagged for manual review (NEEDS_REVIEW? = yes)", boxed=True)
    _legend_row(r + 2, "FFFFFF", "no box", "does not need review (NEEDS_REVIEW? = no)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out_path))
    print(f"[make-final-excel] {n} genes → {out_path}")
    print(f"  {len(headers)} columns; whole-row colour = CONFIDENCE_TIER_HYBRID band")
    for band in ("highest", "high", "medium", "fair", "low", "non-coding"):
        if band_counts.get(band):
            print(f"    {band:11s}: {band_counts[band]}")
    print(f"  {review_n} rows boxed (need review); + Legend sheet")


if __name__ == "__main__":
    main()
