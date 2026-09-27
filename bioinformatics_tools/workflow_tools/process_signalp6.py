#!/usr/bin/env python3
"""Parses SignalP 6.0 `--format none` prediction_results.txt into the standard results.tsv.

Called by margie_sb.smk's run_signalp6 rule, since the SignalP 6.0 module has no
processing entrypoint of its own.

Usage:
    python process_signalp6.py --input prediction_results.txt \
        --output signalp6_results.tsv --organism-name <genome> \
        --tool-used "SignalP 6.0" --command-used "<cmd>" \
        --database-used "<db note>" --input-path <faa> --output-path <raw dir>
"""
import argparse
import csv

# Tab-separated data columns after ID: Prediction, OTHER, SP, LIPO, TAT, TATLIPO, PILIN, CS Position.
_DATA_COLUMNS = (
    "signalp6_prediction", "signalp6_prob_other", "signalp6_prob_sp_sec_spi",
    "signalp6_prob_lipo_sec_spii", "signalp6_prob_tat_spi",
    "signalp6_prob_tatlipo_sec_spii", "signalp6_prob_pilin_sec_spiii",
    "signalp6_cs_position",
)


def _parse_predictions(path: str) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            fields += [""] * (1 + len(_DATA_COLUMNS) - len(fields))
            # feature_id is the first word of the RASTtk FAA header.
            feature_id = fields[0].split()[0] if fields[0].strip() else ""
            row = dict(zip(_DATA_COLUMNS, (v.strip() for v in fields[1:])))
            row["feature_id"] = feature_id
            row["signalp6_has_signal_peptide"] = 0 if row["signalp6_prediction"] == "OTHER" else 1
            rows.append(row)
    return rows


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--organism-name", required=True)
    p.add_argument("--tool-used", required=True)
    p.add_argument("--command-used", required=True)
    p.add_argument("--database-used", required=True)
    p.add_argument("--input-path", required=True)
    p.add_argument("--output-path", required=True)
    args = p.parse_args()

    rows = _parse_predictions(args.input)

    fieldnames = [
        "organism_name", "feature_id", *_DATA_COLUMNS,
        "signalp6_has_signal_peptide", "signalp6_tool_used",
        "signalp6_command_used", "signalp6_database_used",
        "input_path", "output_path",
    ]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row["organism_name"] = args.organism_name
            row["signalp6_tool_used"] = args.tool_used
            row["signalp6_command_used"] = args.command_used
            row["signalp6_database_used"] = args.database_used
            row["input_path"] = args.input_path
            row["output_path"] = args.output_path
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
