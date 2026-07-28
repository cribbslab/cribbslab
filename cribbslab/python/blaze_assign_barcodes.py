#!/usr/bin/env python3
"""
Assign reads to BLAZE cell barcodes (BLAZE step 3) for downstream BAM tagging.

Reads putative_bc.csv and whitelist.csv produced by BLAZE, applies the same
read-to-whitelist matching as blaze.read_assignment (subsequence edit distance,
max ED 2 by default), and writes a TSV of corrected assignments:

    read_id<TAB>CB<TAB>UB

Only successfully assigned reads are written. Unassigned and ambiguous reads
are omitted (matching BLAZE demultiplex behaviour).
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from blaze.config import DEFAULT_ASSIGNMENT_ED
from blaze.read_assignment import _match_bc_row


def load_whitelist(path: str) -> set[str]:
    barcodes: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # BLAZE whitelist.csv and Cell Ranger exports may include GEM suffix.
            token = line.split(",")[0].strip()
            barcodes.append(token.split("-")[0])
    return set(barcodes)


def assign_reads(
    putative_bc_csv: str,
    whitelist_csv: str,
    output_tsv: str,
    max_ed: int = DEFAULT_ASSIGNMENT_ED,
    min_q: int = 0,
) -> tuple[int, int]:
    df = pd.read_csv(putative_bc_csv)
    df = df.fillna("")
    whitelist = load_whitelist(whitelist_csv)

    assigned = 0
    with open(output_tsv, "w", encoding="utf-8") as out:
        for row in df.itertuples(index=False):
            bc, umi, _strand = _match_bc_row(row, whitelist, max_ed, min_q)
            if bc:
                out.write(f"{row.read_id}\t{bc}\t{umi}\n")
                assigned += 1

    return assigned, len(df)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assign reads to BLAZE whitelist barcodes for BAM CB/UB tagging."
    )
    parser.add_argument(
        "--putative-bc",
        required=True,
        help="BLAZE putative_bc.csv from step 1",
    )
    parser.add_argument(
        "--whitelist",
        required=True,
        help="BLAZE whitelist.csv from step 2",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV: read_id, corrected CB, corrected UMI",
    )
    parser.add_argument(
        "--max-ed",
        type=int,
        default=DEFAULT_ASSIGNMENT_ED,
        help=f"Maximum subsequence edit distance (default: {DEFAULT_ASSIGNMENT_ED})",
    )
    parser.add_argument(
        "--min-q",
        type=int,
        default=0,
        help="Minimum putative barcode minQ for assignment (default: 0, match BLAZE demux)",
    )
    args = parser.parse_args(argv)

    assigned, total = assign_reads(
        args.putative_bc,
        args.whitelist,
        args.output,
        max_ed=args.max_ed,
        min_q=args.min_q,
    )
    print(
        f"BLAZE read assignment: {assigned}/{total} reads assigned to whitelist barcodes",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
