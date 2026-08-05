#!/usr/bin/env python3
"""Parse SignalP 6.0's --format none prediction_results.txt into the
project's standard per-protein results.tsv shape (margie_sb phase6).

signalp6 has no build-here container/entrypoint (HPC envmodule wrapping
a pre-built image we don't own) -- this script is the missing
"processing" step that phobius/tmbed get for free from their own
containers' /usr/local/bin/run. Invoked directly from margie_sb.smk's
run_signalp6 rule; writes into signalp6's own processed/ dir first, then
run_signalp6 cp's the result into the final Snakemake output, same as
every other tool.

Usage:
    python process_signalp6.py --input prediction_results.txt \
        --output signalp6_results.tsv --organism-name <genome> \
        --tool-used "SignalP 6.0" --command-used "<cmd>" \
        --database-used "<db note>" --input-path <faa> --output-path <raw dir>
"""
import argparse
import csv

# --format none's two header lines start with "#"; data rows are tab-
# separated: ID  Prediction  OTHER  SP(Sec/SPI)  LIPO(Sec/SPII)
# TAT(Tat/SPI)  TATLIPO(Sec/SPII)  PILIN(Sec/SPIII)  CS Position
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
            # RASTtk FAA headers are "fig|...peg.N description [...]" --
            # feature_id is just the first whitespace-delimited token.
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
