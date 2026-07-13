"""Load derived AnnData matrices and barcode QC tables."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd

from isoquant_report.config import ReportConfig

BACKED_THRESHOLD_BYTES = 500 * 1024 * 1024  # 500 MB


def load_h5ad(path: str, backed: bool = True) -> Optional[ad.AnnData]:
    """Load AnnData, using backed mode for large files."""
    if not path or not os.path.exists(path):
        return None
    size = os.path.getsize(path)
    mode = "r" if backed and size > BACKED_THRESHOLD_BYTES else None
    return ad.read_h5ad(path, backed=mode)


def load_barcode_qc(path: str) -> Optional[pd.DataFrame]:
    """Load per-barcode QC table from pipeline splice_matrices output."""
    if not path or not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def load_matrices(config: ReportConfig) -> dict:
    """
    Load all derived matrix artefacts referenced in config.

    Returns dict with keys: gene_adata, transcript_adata, barcode_qc (optional).
    """
    return {
        "gene_adata": load_h5ad(config.gene_h5ad),
        "transcript_adata": load_h5ad(config.transcript_h5ad),
        "barcode_qc": load_barcode_qc(config.barcode_qc),
    }
