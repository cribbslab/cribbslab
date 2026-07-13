"""Parse IsoQuant read assignment tables."""

from __future__ import annotations

import gzip
import re
from typing import Optional

import numpy as np
import pandas as pd

from isoquant_report.isoquant_io.allinfo import load_allinfo_barcodes
from isoquant_report.isoquant_io.discover import IsoQuantPaths

SQANTI_MAP = {
    "full_splice_match": "FSM",
    "incomplete_splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNIC",
    "genic": "genic",
    "antisense": "antisense",
    "intergenic": "intergenic",
}


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _scan_header(path: str) -> tuple[list[str], int]:
    with _open(path) as fh:
        for n, raw_line in enumerate(fh, start=1):
            parts = raw_line.lstrip("#").strip().split("\t")
            if len(parts) > 1:
                return [p.lstrip("#").strip() for p in parts], n
    raise ValueError(f"No header in {path}")


def _extract_classification(row: pd.Series) -> Optional[str]:
    """Parse SQANTI-style category from additional_info or classification."""
    for col in ("classification", "additional_info", "additional"):
        if col not in row.index or pd.isna(row[col]):
            continue
        val = str(row[col]).lower()
        for key, label in SQANTI_MAP.items():
            if key in val:
                return label
    return None


def load_assignments(paths: IsoQuantPaths) -> pd.DataFrame:
    """
    Load per-read assignments with barcode when available.

    Joins read_assignments with allinfo for barcode/umi when needed.
    """
    read_path = paths.read_info or paths.read_assignments
    if not read_path:
        raise FileNotFoundError("No read assignment file found")

    col_names, skip = _scan_header(read_path)
    df = pd.read_csv(
        read_path,
        sep="\t",
        header=None,
        names=col_names,
        skiprows=skip,
        dtype=str,
        low_memory=False,
    )
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    # Normalise column names
    rename = {}
    if "isoform_assignment_type" in df.columns:
        rename["isoform_assignment_type"] = "assignment_type"
    if "isoform_assignment_events" in df.columns:
        rename["isoform_assignment_events"] = "assignment_events"
    df = df.rename(columns=rename)

    if "assignment_events" not in df.columns:
        df["assignment_events"] = ""

    # Barcode join via allinfo
    if "barcode" not in df.columns and paths.allinfo:
        allinfo = load_allinfo_barcodes(paths.allinfo)
        if "read_id" in df.columns and not allinfo.empty:
            df["read_id"] = df["read_id"].astype(str)
            allinfo["read_id"] = allinfo["read_id"].astype(str)
            df = df.merge(
                allinfo[["read_id", "barcode", "umi"]],
                on="read_id",
                how="left",
            )
    elif "groups" in df.columns and "barcode" not in df.columns:
        # groups column often holds barcode in SC mode
        df["barcode"] = df["groups"].str.split(",").str[0]

    df["structural_category"] = df.apply(_extract_classification, axis=1)

    # Novel vs known from isoform_id pattern or read_type
    def _novel_flag(row: pd.Series) -> str:
        if "read_type" in row.index and str(row.get("read_type", "")) == "novel":
            return "novel"
        iso = str(row.get("isoform_id", ""))
        if iso and (iso.startswith("novel") or "novel" in iso.lower()):
            return "novel"
        if row.get("structural_category") in ("NIC", "NNIC"):
            return "novel"
        return "known"

    df["novel_known"] = df.apply(_novel_flag, axis=1)
    return df
