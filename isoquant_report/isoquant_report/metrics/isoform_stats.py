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


def _n_features_detected(adata: ad.AnnData) -> np.ndarray:
    """Per-cell count of features with >0 counts."""
    det = (adata.X > 0).sum(axis=1)
    if hasattr(det, "A1"):
        return det.A1
    return np.asarray(det).flatten()


def gene_isoform_detection(
    gene_adata: ad.AnnData,
    tx_adata: ad.AnnData,
) -> pd.DataFrame:
    """
    Per-cell genes detected vs isoforms detected.

    Aligns on the intersection of barcodes present in both matrices.
    Gene and transcript h5ad files may differ in cell sets when spliced
    transcript counts omit cells with no spliced molecules.
    """
    common = sorted(set(gene_adata.obs_names) & set(tx_adata.obs_names))
    if not common:
        return pd.DataFrame(
            columns=["barcode", "genes_detected", "isoforms_detected"]
        )

    gene_sub = gene_adata[common, :]
    tx_sub = tx_adata[common, :]
    return pd.DataFrame(
        {
            "barcode": common,
            "genes_detected": _n_features_detected(gene_sub),
            "isoforms_detected": _n_features_detected(tx_sub),
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
