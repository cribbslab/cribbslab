#!/usr/bin/env python3
"""
ctat_fusion_format.py -- Join ctat-LR-fusion output to cell barcodes from tagged BAM.

Mirrors wf-single-cell ``format_ctat_output`` (workflow-glue) using the tagged
BAM CB/UB tags instead of read_summary_tags.tsv.

Usage
-----
    python ctat_fusion_format.py \\
        --predictions fusions/sample/ctat-LR-fusion.fusion_predictions.tsv \\
        --bam tagged/sample.tagged.bam \\
        --sample sample \\
        --out-read fusions/sample.ctat-LR-fusion.fusion_predictions_per-read.tsv \\
        --out-fusion fusions/sample.ctat-LR-fusion.fusion_predictions_per-fusion.tsv
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import pandas as pd
import pysam


def _read_fusion_read_map(predictions_path: str, ctat_dir: str) -> dict[str, set[str]]:
    """
    Map fusion name -> supporting read names.

    CTAT may list reads in fusion_predictions.tsv (``#FusionName`` + read list
    columns) or in companion ``*fusion_reads*`` / ``*supporting_reads*`` files.
    """
    fusion_reads: dict[str, set[str]] = defaultdict(set)

    # Companion read-level files in the ctat output directory
    for pattern in (
        "*fusion_reads*.tsv",
        "*supporting_reads*.tsv",
        "*fusion_reads*.txt",
    ):
        for path in glob.glob(os.path.join(ctat_dir, pattern)):
            try:
                df = pd.read_csv(path, sep="\t", comment="#")
            except Exception:
                continue
            read_col = None
            fusion_col = None
            for c in df.columns:
                cl = str(c).lower()
                if read_col is None and ("read" in cl or cl == "readname"):
                    read_col = c
                if fusion_col is None and "fusion" in cl:
                    fusion_col = c
            if read_col is None or fusion_col is None:
                continue
            for _, row in df.iterrows():
                fusion_reads[str(row[fusion_col])].add(str(row[read_col]))

    if not os.path.exists(predictions_path):
        return fusion_reads

    pred = pd.read_csv(predictions_path, sep="\t", comment="#")
    if pred.empty:
        return fusion_reads

    name_col = None
    for c in pred.columns:
        if str(c).lower() in ("fusionname", "#fusionname", "fusion_name"):
            name_col = c
            break
    if name_col is None:
        name_col = pred.columns[0]

    read_cols = [
        c
        for c in pred.columns
        if re.search(r"read|junction|split", str(c), re.I)
    ]

    for _, row in pred.iterrows():
        fusion = str(row[name_col])
        for c in read_cols:
            val = row[c]
            if pd.isna(val) or str(val).strip() in ("", ".", "NA"):
                continue
            for token in re.split(r"[,;\s]+", str(val)):
                token = token.strip()
                if token and token not in (".", "NA"):
                    fusion_reads[fusion].add(token)

    return fusion_reads


def _bam_read_tags(bam_path: str) -> dict[str, tuple[str, str]]:
    """Return read_name -> (CB, UB) from tagged BAM."""
    tags: dict[str, tuple[str, str]] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            cb = read.get_tag("CB") if read.has_tag("CB") else ""
            ub = read.get_tag("UB") if read.has_tag("UB") else ""
            tags[read.query_name] = (str(cb), str(ub))
    return tags


def format_ctat_output(
    predictions_path: str,
    bam_path: str,
    sample: str,
    out_read: str,
    out_fusion: str,
) -> None:
    """Write per-read and per-fusion summary tables."""
    ctat_dir = os.path.dirname(predictions_path) or "."
    fusion_reads = _read_fusion_read_map(predictions_path, ctat_dir)
    read_tags = _bam_read_tags(bam_path)

    pred = pd.DataFrame()
    if os.path.exists(predictions_path):
        pred = pd.read_csv(predictions_path, sep="\t", comment="#")

    name_col = None
    left_col = right_col = None
    if not pred.empty:
        for c in pred.columns:
            cl = str(c).lower()
            if name_col is None and "fusion" in cl and "name" in cl:
                name_col = c
            if left_col is None and "left" in cl and "gene" in cl:
                left_col = c
            if right_col is None and "right" in cl and "gene" in cl:
                right_col = c
        if name_col is None:
            name_col = pred.columns[0]

    meta = {}
    if not pred.empty and name_col is not None:
        for _, row in pred.iterrows():
            fusion = str(row[name_col])
            meta[fusion] = {
                "left_gene": str(row[left_col]) if left_col else "",
                "right_gene": str(row[right_col]) if right_col else "",
            }

    per_read_rows = []
    per_fusion_rows = []

    all_fusions = set(fusion_reads.keys()) | set(meta.keys())
    for fusion in sorted(all_fusions):
        reads = fusion_reads.get(fusion, set())
        m = meta.get(fusion, {"left_gene": "", "right_gene": ""})
        cells: set[str] = set()
        umis: set[str] = set()
        n_reads = 0
        for rid in reads:
            cb, ub = read_tags.get(rid, ("", ""))
            if not cb:
                continue
            n_reads += 1
            cells.add(cb)
            if ub:
                umis.add(f"{cb}:{ub}")
            per_read_rows.append(
                {
                    "sample": sample,
                    "fusion_name": fusion,
                    "read_id": rid,
                    "barcode": cb,
                    "umi": ub,
                    "left_gene": m["left_gene"],
                    "right_gene": m["right_gene"],
                }
            )
        per_fusion_rows.append(
            {
                "sample": sample,
                "fusion_name": fusion,
                "left_gene": m["left_gene"],
                "right_gene": m["right_gene"],
                "n_reads": n_reads,
                "n_cells": len(cells),
                "n_umis": len(umis),
            }
        )

    os.makedirs(os.path.dirname(out_read) or ".", exist_ok=True)
    pd.DataFrame(per_read_rows).to_csv(out_read, sep="\t", index=False)
    pd.DataFrame(per_fusion_rows).to_csv(out_fusion, sep="\t", index=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-read", required=True)
    parser.add_argument("--out-fusion", required=True)
    args = parser.parse_args(argv)

    format_ctat_output(
        args.predictions,
        args.bam,
        args.sample,
        args.out_read,
        args.out_fusion,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
