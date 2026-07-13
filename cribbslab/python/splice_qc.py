#!/usr/bin/env python3
"""
splice_qc.py -- Per-barcode splicing proportion QC with MultiQC output.

Concatenates per-sample barcode_qc.tsv tables produced by
isoquant_matrices.py, flags outliers on unspliced_fraction using the median
absolute deviation (MAD) method, cross-tabulates the flag against total_umis
and n_genes bins, and writes a MultiQC custom-content file.

Usage
-----
    python splice_qc.py --qc-tables s1.barcode_qc.tsv s2.barcode_qc.tsv
                        --mad-multiplier 3.0
                        --outfile qc_splice/splice_proportion_mqc.tsv
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# MAD outlier detection
# ---------------------------------------------------------------------------

def mad_outlier_flag(series, multiplier):
    """
    Flag values where |x - median| > multiplier * MAD.

    When MAD is zero (all values identical) no outliers are flagged.

    Returns a boolean Series, True = outlier.
    """
    median = series.median()
    mad = (series - median).abs().median()
    if mad == 0:
        return pd.Series(False, index=series.index)
    return (series - median).abs() > multiplier * mad


# ---------------------------------------------------------------------------
# MultiQC custom-content helpers
# ---------------------------------------------------------------------------

_MQC_HEADER = """\
# id: splice_proportion_qc
# section_name: 'Splicing Proportion QC'
# description: >
#   Per-barcode unspliced fraction flagged by MAD outlier detection.
#   Derived from IsoQuant single-cell molecule assignments.
#   Outlier column: 1 = flagged, 0 = pass.
# format: 'tsv'
# plot_type: 'table'
# pconfig:
#     id: 'splice_proportion_table'
#     title: 'Splice Proportion QC'
"""


def write_multiqc_tsv(df, outfile):
    """
    Write a MultiQC custom-content TSV.

    The file header uses MultiQC's comment-block format so the table is
    picked up automatically when multiqc scans the working directory.
    """
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    with open(outfile, "w") as fh:
        fh.write(_MQC_HEADER)
        df.to_csv(fh, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Cross-tabulation helpers
# ---------------------------------------------------------------------------

def _umi_bin(series, bins=(0, 100, 500, 2000, np.inf)):
    labels = [
        "<{:g}".format(bins[i + 1]) if np.isfinite(bins[i + 1])
        else ">={:g}".format(bins[i])
        for i in range(len(bins) - 1)
    ]
    return pd.cut(series, bins=list(bins), labels=labels, right=False)


def _gene_bin(series, bins=(0, 50, 200, 1000, np.inf)):
    labels = [
        "<{:g}".format(bins[i + 1]) if np.isfinite(bins[i + 1])
        else ">={:g}".format(bins[i])
        for i in range(len(bins) - 1)
    ]
    return pd.cut(series, bins=list(bins), labels=labels, right=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Per-barcode splicing-proportion QC with MultiQC output."
    )
    parser.add_argument(
        "--qc-tables", nargs="+", required=True,
        help="Per-sample barcode_qc.tsv files from isoquant_matrices.py.",
    )
    parser.add_argument(
        "--mad-multiplier", type=float, default=3.0,
        help="MAD multiplier for unspliced_fraction outlier detection "
             "(default 3.0).",
    )
    parser.add_argument(
        "--outfile", required=True,
        help="Output MultiQC custom-content TSV (filename must end in _mqc.tsv).",
    )
    args = parser.parse_args(argv)

    # ---- Load and concatenate per-sample QC tables --------------------------
    dfs = []
    for path in args.qc_tables:
        if not os.path.exists(path):
            print(
                "WARNING: QC table not found, skipping: {}".format(path),
                file=sys.stderr,
            )
            continue
        df = pd.read_csv(path, sep="\t")
        # Infer sample name from filename: {sample}.barcode_qc.tsv
        sample_name = os.path.basename(path).replace(".barcode_qc.tsv", "")
        df.insert(0, "sample", sample_name)
        dfs.append(df)

    if not dfs:
        print("ERROR: No QC tables loaded.", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(dfs, ignore_index=True)
    print(
        "Loaded {} barcodes from {} sample(s).".format(
            len(combined), len(dfs)),
        file=sys.stderr,
    )

    # ---- Flag outliers on unspliced_fraction --------------------------------
    flagged_per_sample = []
    for sample, grp in combined.groupby("sample"):
        flags = mad_outlier_flag(grp["unspliced_fraction"], args.mad_multiplier)
        flagged_per_sample.append(flags)
    combined["outlier"] = pd.concat(flagged_per_sample).astype(int)

    n_flagged = combined["outlier"].sum()
    n_total = len(combined)
    print(
        "Flagged {} / {} barcodes as unspliced_fraction outliers "
        "(MAD multiplier={}).".format(n_flagged, n_total, args.mad_multiplier),
        file=sys.stderr,
    )

    # ---- Cross-tabulation: outlier vs UMI bins and gene bins ----------------
    combined["umi_bin"] = _umi_bin(combined["total_umis"])
    combined["gene_bin"] = _gene_bin(combined["n_genes"])

    xtab_umi = (
        combined.groupby(["sample", "umi_bin"])["outlier"]
        .agg(n_barcodes="count", n_outliers="sum")
        .reset_index()
    )
    xtab_umi["outlier_pct"] = (
        100.0 * xtab_umi["n_outliers"] / xtab_umi["n_barcodes"]
    ).round(2)

    xtab_gene = (
        combined.groupby(["sample", "gene_bin"])["outlier"]
        .agg(n_barcodes="count", n_outliers="sum")
        .reset_index()
    )
    xtab_gene["outlier_pct"] = (
        100.0 * xtab_gene["n_outliers"] / xtab_gene["n_barcodes"]
    ).round(2)

    # Print cross-tabs to stderr for logging
    print("\nOutlier rate by UMI bin:\n{}\n".format(
        xtab_umi.to_string(index=False)), file=sys.stderr)
    print("Outlier rate by gene bin:\n{}\n".format(
        xtab_gene.to_string(index=False)), file=sys.stderr)

    # ---- Build MultiQC output table -----------------------------------------
    # Summarise at sample level for the MultiQC table; per-barcode detail
    # would be too large for the report.
    summary = (
        combined.groupby("sample")
        .agg(
            n_barcodes=("barcode", "count"),
            median_total_umis=("total_umis", "median"),
            median_unspliced_fraction=("unspliced_fraction", "median"),
            n_outliers=("outlier", "sum"),
        )
        .reset_index()
    )
    summary["outlier_pct"] = (
        100.0 * summary["n_outliers"] / summary["n_barcodes"]
    ).round(2)
    summary["median_unspliced_fraction"] = summary[
        "median_unspliced_fraction"
    ].round(4)

    write_multiqc_tsv(summary, args.outfile)
    print("Written MultiQC table: {}".format(args.outfile), file=sys.stderr)

    # Write full per-barcode table alongside for downstream use
    detail_path = args.outfile.replace("_mqc.tsv", "_detail.tsv")
    out_cols = [
        "sample", "barcode", "total_umis", "spliced_umis", "unspliced_umis",
        "ambiguous_umis", "unspliced_fraction", "n_genes", "outlier",
    ]
    available = [c for c in out_cols if c in combined.columns]
    combined[available].to_csv(detail_path, sep="\t", index=False)
    print("Written per-barcode detail: {}".format(detail_path), file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
