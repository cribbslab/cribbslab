"""Isoform-level statistics and gene/isoform comparison metrics."""

from __future__ import annotations

from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def isoforms_per_gene(
    assignments: pd.DataFrame,
    gtf_models: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Count distinct isoforms per gene from assignments or GTF."""
    if gtf_models is not None and not gtf_models.empty:
        return (
            gtf_models.groupby("gene_id")["transcript_id"]
            .nunique()
            .reset_index()
            .rename(columns={"transcript_id": "n_isoforms"})
        )
    if "gene_id" not in assignments.columns or "isoform_id" not in assignments.columns:
        return pd.DataFrame()
    df = assignments.dropna(subset=["gene_id", "isoform_id"])
    return (
        df.groupby("gene_id")["isoform_id"]
        .nunique()
        .reset_index()
        .rename(columns={"isoform_id": "n_isoforms"})
    )


def gene_isoform_detection(
    gene_adata: ad.AnnData,
    tx_adata: ad.AnnData,
) -> pd.DataFrame:
    """Per-cell genes detected vs isoforms detected."""
    gene_det = (gene_adata.X > 0).sum(axis=1)
    if hasattr(gene_det, "A1"):
        gene_det = gene_det.A1
    else:
        gene_det = np.asarray(gene_det).flatten()
    tx_det = (tx_adata.X > 0).sum(axis=1)
    if hasattr(tx_det, "A1"):
        tx_det = tx_det.A1
    else:
        tx_det = np.asarray(tx_det).flatten()
    return pd.DataFrame(
        {
            "barcode": list(gene_adata.obs_names),
            "genes_detected": gene_det,
            "isoforms_detected": tx_det,
        }
    )


def tpm_agreement(
    native_totals: np.ndarray,
    derived_totals: np.ndarray,
    barcodes: list[str],
) -> pd.DataFrame:
    """Scatter data for native vs derived count agreement."""
    return pd.DataFrame(
        {
            "barcode": barcodes,
            "native_total": native_totals,
            "derived_total": derived_totals,
        }
    )
