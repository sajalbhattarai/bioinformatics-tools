"""
Build a GFF3 file from the mechanical FINAL-scored-labeled-genes-annotated.tsv.

Each CDS feature carries these attributes (GFF3 percent-encoded):
    ID               – feature_id  (fig|..peg.N)
    gene_id          – gene_id     (e.g. NC_000908.2_686+1143)
    product_descriptor – winning functional descriptor
    fingerprint      – fingerprint_hash (16-char hex)
    confidence_score – final blended confidence score
    flagging         – confidence_flag; flag_reason (if any)

seqname (column 1) is looked up from rast.gff via the ID= attribute,
because rast.gff is the authoritative source for per-feature contig IDs.

Usage:
    python make-gff.py \\
        --final   <FINAL-scored-labeled-genes-annotated.tsv> \\
        --rast-gff <rasttk/rast.gff> \\
        --output  <annotation.gff3>
"""
import argparse
import csv
import sys
from urllib.parse import quote

csv.field_size_limit(10_000_000)

_GFF3_SAFE = " ,:/.-_@|#"  # chars kept as-is inside GFF3 attribute values


def _encode(value: str) -> str:
    """Percent-encode a GFF3 attribute value, preserving common safe chars."""
    return quote(value, safe=_GFF3_SAFE)


def _build_seqname_lookup(rast_gff_path):
    """Return {feature_id: seqname} from rast.gff CDS lines."""
    lookup = {}
    with open(rast_gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqname, _, feat_type = cols[0], cols[1], cols[2]
            if feat_type != "CDS":
                continue
            attrs = cols[8]
            for part in attrs.split(";"):
                if part.startswith("ID="):
                    fid = part[3:].strip()
                    lookup[fid] = seqname
                    break
    return lookup


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--final",    required=True, help="FINAL-scored-labeled-genes-annotated.tsv")
    ap.add_argument("--rast-gff", required=True, help="rasttk/rast.gff")
    ap.add_argument("--output",   required=True, help="output .gff3 path")
    args = ap.parse_args()

    seqname_map = _build_seqname_lookup(args.rast_gff)

    with open(args.final, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    missing_seqnames = 0
    out_lines = ["##gff-version 3"]

    for row in rows:
        fid            = row.get("feature_id", "").strip()
        gene_id        = row.get("gene_id", "").strip()
        gene_start_raw = row.get("gene_start", "").strip()
        gene_end_raw   = row.get("gene_end", "").strip()
        strand         = row.get("RAST_strand", ".").strip() or "."
        descriptor     = row.get("best_consensus_product_descriptor", "").strip()
        fp_hash        = row.get("fingerprint_hash", "").strip()
        conf_score     = row.get("confidence_score", ".").strip() or "."
        conf_flag      = row.get("confidence_flag", "").strip()
        flag_reason    = row.get("flag_reason", "").strip()

        seqname = seqname_map.get(fid)
        if seqname is None:
            missing_seqnames += 1
            seqname = "."

        try:
            start = int(gene_start_raw)
            end   = int(gene_end_raw)
        except ValueError:
            print(f"[make-gff] skipping {fid}: invalid start/end ({gene_start_raw}/{gene_end_raw})",
                  file=sys.stderr)
            continue

        flagging = conf_flag
        if flag_reason:
            flagging = f"{conf_flag}; {flag_reason}"

        attrs = (
            f"ID={_encode(fid)}"
            f";gene_id={_encode(gene_id)}"
            f";product_descriptor={_encode(descriptor)}"
            f";fingerprint={_encode(fp_hash)}"
            f";confidence_score={_encode(conf_score)}"
            f";flagging={_encode(flagging)}"
        )

        line = "\t".join([
            seqname,
            "margie",
            "CDS",
            str(start),
            str(end),
            ".",
            strand,
            ".",
            attrs,
        ])
        out_lines.append(line)

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")

    if missing_seqnames:
        print(f"[make-gff] WARNING: {missing_seqnames} feature(s) had no seqname in rast.gff "
              f"(written as '.')", file=sys.stderr)
    print(f"[make-gff] wrote {len(out_lines) - 1} CDS features → {args.output}")


if __name__ == "__main__":
    main()
