"""Append envelope classification context to a phase8 tool's results.tsv.

Envelope's decision (envelope_type, inference_basis, evidence_json) is
genome-level, not per-protein, so the same three values get repeated on
every row -- giving anyone reading deepsig/psortb/signalp4's own output
table direct visibility into *why* the -k/-organism gram-class flag was
set the way it was (e.g. a real detected call vs. the conservative
tie-breaker default used for wall-less organisms like Mycoplasma), without
cross-referencing envelope's separate output file.

Usage:
    python enrich_with_envelope.py --input <tool>_results.tsv \
        --envelope-summary envelope_summary.tsv --output <enriched>.tsv
"""
import argparse
import csv


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--envelope-summary", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    with open(args.envelope_summary, newline="", encoding="utf-8") as f:
        envelope_row = next(csv.DictReader(f, delimiter="\t"))

    extra = {
        "ENVELOPE_envelope_type": envelope_row.get("envelope_type", ""),
        "ENVELOPE_inference_basis": envelope_row.get("inference_basis", ""),
        "ENVELOPE_evidence_json": envelope_row.get("evidence_json", ""),
    }

    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames) + list(extra.keys())
        rows = list(reader)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row.update(extra)
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
