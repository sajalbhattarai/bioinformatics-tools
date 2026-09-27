"""Appends the genome-level envelope classification to a tool's results.tsv.

The three ENVELOPE_ values repeat on every row, so deepsig/psortb/signalp4
tables show which gram-class flag was used and why.

Usage:
    python enrich_with_envelope.py --input <tool>_results.tsv \
        --envelope-summary envelope_summary.tsv --output <enriched>.tsv
"""
import argparse
import csv


def main() -> int:
    """Reads the envelope summary row and writes the tool table with ENVELOPE_ columns added, using csv."""
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
