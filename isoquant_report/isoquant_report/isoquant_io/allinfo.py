"""Parse IsoQuant allinfo (post-dedup molecule) files."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

ALLINFO_COLUMNS = [
    "read_id",
    "gene_id",
    "cell_type",
    "barcode",
    "umi",
    "introns",
    "TSS",
    "polyA",
    "exons",
    "read_type",
    "intron_count",
    "transcript_id",
    "transcript_type",
]


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _scan_tsv_header(fh) -> tuple[list[str], int]:
    """Skip comment lines; return column names and lines skipped."""
    for n, raw_line in enumerate(fh, start=1):
        parts = raw_line.lstrip("#").strip().split("\t")
        if len(parts) > 1:
            if parts[0].lower() in ("read_id", "barcode") or "read_id" in [
                p.lower() for p in parts
            ]:
                return [p.lstrip("#").strip() for p in parts], n
            # Data row without header
            if n == 1:
                n_cols = len(parts)
                cols = (
                    list(ALLINFO_COLUMNS)
                    if n_cols == len(ALLINFO_COLUMNS)
                    else list(ALLINFO_COLUMNS[:5])
                    + [f"col{i}" for i in range(5, n_cols)]
                )
                return cols, 0
    raise ValueError("No tab-separated content in allinfo file")


def load_allinfo_barcodes(path: str) -> pd.DataFrame:
    """
    Load read_id, barcode, umi from an IsoQuant UMI_filtered allinfo file.

    The file often has no header row.
    """
    with _open(path) as fh:
        col_names, skip = _scan_tsv_header(fh)

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=col_names,
        skiprows=skip if skip > 0 else None,
        dtype=str,
        low_memory=False,
        comment="#",
    )
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    rename = {}
    for src, dst in [
        ("read_id", "read_id"),
        ("barcode", "barcode"),
        ("umi", "umi"),
        ("UMI", "umi"),
    ]:
        if src in df.columns:
            rename[src] = dst
    df = df.rename(columns=rename)

    for col in ["read_id", "barcode", "umi"]:
        if col not in df.columns:
            df[col] = np.nan
        else:
            df[col] = df[col].astype(str)

    return df[["read_id", "barcode", "umi"]]
