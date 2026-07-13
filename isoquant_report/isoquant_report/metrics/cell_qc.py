"""Single-cell count matrix QC metrics."""

from __future__ import annotations

from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def compute_qc_metrics(
    adata: ad.AnnData,
    mito_prefix: str = "MT-",
    ribo_prefixes: Optional[list[str]] = None,
) -> ad.AnnData:
    """Add standard scanpy QC metrics to adata.obs."""
    ribo_prefixes = ribo_prefixes or ["RPS", "RPL"]
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    ribo_mask = False
    for p in ribo_prefixes:
        ribo_mask = ribo_mask | adata.var_names.str.startswith(p)
    adata.var["ribo"] = ribo_mask
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def apply_qc_filters(
    adata: ad.AnnData,
    min_counts: int = 500,
    min_genes: int = 200,
    max_mito: float = 0.2,
    max_unspliced_fraction: float = 1.0,
    barcode_qc: Optional[pd.DataFrame] = None,
) -> np.ndarray:
    """
    Return boolean mask of cells passing QC thresholds.

    unspliced_fraction filter uses barcode_qc if provided.
    """
    mask = np.ones(adata.n_obs, dtype=bool)
    if "total_counts" in adata.obs.columns:
        mask &= adata.obs["total_counts"].values >= min_counts
    if "n_genes_by_counts" in adata.obs.columns:
        mask &= adata.obs["n_genes_by_counts"].values >= min_genes
    if "pct_counts_mt" in adata.obs.columns:
        mask &= adata.obs["pct_counts_mt"].values <= max_mito * 100

    if barcode_qc is not None and "unspliced_fraction" in barcode_qc.columns:
        bc_map = barcode_qc.set_index("barcode")["unspliced_fraction"]
        uf = adata.obs_names.to_series().map(bc_map).fillna(0).values
        mask &= uf <= max_unspliced_fraction

    return mask


def knee_plot_data(adata: ad.AnnData) -> pd.DataFrame:
    """Barcode rank vs total counts for knee plot."""
    if "total_counts" not in adata.obs.columns:
        return pd.DataFrame()
    counts = adata.obs["total_counts"].sort_values(ascending=False).values
    return pd.DataFrame({"rank": np.arange(1, len(counts) + 1), "counts": counts})
