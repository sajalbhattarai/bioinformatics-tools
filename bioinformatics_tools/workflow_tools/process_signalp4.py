#!/usr/bin/env python3
"""Parse SignalP 4.1's `-f short` stdout into the project's standard
per-protein results.tsv shape (margie_sb phase8, envelope-dependent).

signalp4 has no build-here container/entrypoint (HPC envmodule wrapping
a pre-built image we don't own) -- this script is the missing
"processing" step that deepsig/psortb get for free from their own
containers' /usr/local/bin/run. Invoked directly from margie_sb.smk's
run_signalp4 rule; writes into signalp4's own processed/ dir first, then
run_signalp4 cp's the result into the final Snakemake output, same as
every other tool. Envelope enrichment (ENVELOPE_*) is added afterwards
by enrich_with_envelope.py in load_signalp4_to_db, same as deepsig/psortb.

Usage:
    python process_signalp4.py --input signalp4_out.txt \
        --output signalp4_results.tsv --organism-name <genome> \
        --gram-class <gram+/gram-> --command-used "<cmd>" \
        --database-used "<db note>" --input-path <faa> --output-path <raw dir>
"""
import argparse
import csv

# SignalP 4.1's "-f short" data-row columns, in order, after the name
# column (which is the FASTA header up to the first whitespace, same
# truncation RASTtk's fig|...peg.N ids already have -- no split needed).
_DATA_COLUMNS = (
    "signalp4_cmax", "signalp4_cmax_pos", "signalp4_ymax", "signalp4_ymax_pos",
    "signalp4_smax", "signalp4_smax_pos", "signalp4_smean", "signalp4_d",
    "signalp4_prediction", "signalp4_d_cutoff", "signalp4_networks_used",
)


def _parse_predictions(path: str) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 1 + len(_DATA_COLUMNS):
                continue
            row = dict(zip(_DATA_COLUMNS, fields[1:1 + len(_DATA_COLUMNS)]))
            row["feature_id"] = fields[0]
            row["signalp4_has_signal_peptide"] = 1 if row["signalp4_prediction"] == "Y" else 0
            rows.append(row)
    return rows


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--organism-name", required=True)
    p.add_argument("--gram-class", required=True)
    p.add_argument("--command-used", required=True)
    p.add_argument("--database-used", required=True)
    p.add_argument("--input-path", required=True)
    p.add_argument("--output-path", required=True)
    args = p.parse_args()

    rows = _parse_predictions(args.input)

    fieldnames = [
        "organism_name", "feature_id", *_DATA_COLUMNS,
        "signalp4_has_signal_peptide", "signalp4_gram_class",
        "signalp4_command_used", "signalp4_database_used",
        "input_path", "output_path",
    ]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row["organism_name"] = args.organism_name
            row["signalp4_gram_class"] = args.gram_class
            row["signalp4_command_used"] = args.command_used
            row["signalp4_database_used"] = args.database_used
            row["input_path"] = args.input_path
            row["output_path"] = args.output_path
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
