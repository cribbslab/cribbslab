"""Reconcile IsoQuant native grouped counts against derived AnnData matrices."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from isoquant_report.isoquant_io.counts import load_grouped_matrix
from isoquant_report.isoquant_io.discover import IsoQuantPaths


@dataclass
class ReconcileResult:
    """Summary of agreement between native and derived matrices."""

    n_cells_native: int = 0
    n_cells_derived: int = 0
    n_cells_overlap: int = 0
    cell_jaccard: float = 0.0
    n_features_native: int = 0
    n_features_derived: int = 0
    n_features_overlap: int = 0
    feature_jaccard: float = 0.0
    cells_only_native: list[str] = None
    cells_only_derived: list[str] = None
    features_only_native: list[str] = None
    features_only_derived: list[str] = None
    total_count_correlation: float = float("nan")
    mean_abs_rel_diff: float = float("nan")

    def __post_init__(self):
        self.cells_only_native = self.cells_only_native or []
        self.cells_only_derived = self.cells_only_derived or []
        self.features_only_native = self.features_only_native or []
        self.features_only_derived = self.features_only_derived or []


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def reconcile_gene_matrix(
    paths: IsoQuantPaths,
    gene_adata: Optional[ad.AnnData],
    layer: str = "total",
) -> ReconcileResult:
    """
    Compare IsoQuant native gene grouped counts to derived AnnData.

    Uses the specified layer if present, else X.
    """
    result = ReconcileResult()
    if gene_adata is None:
        return result

    try:
        native_mat, native_features, native_barcodes = load_grouped_matrix(
            paths, "gene"
        )
    except FileNotFoundError:
        return result

    derived_barcodes = list(gene_adata.obs_names)
    derived_features = list(gene_adata.var_names)

    set_native_bc = set(native_barcodes)
    set_derived_bc = set(derived_barcodes)
    set_native_feat = set(native_features)
    set_derived_feat = set(derived_features)

    result.n_cells_native = len(set_native_bc)
    result.n_cells_derived = len(set_derived_bc)
    result.n_cells_overlap = len(set_native_bc & set_derived_bc)
    result.cell_jaccard = _jaccard(set_native_bc, set_derived_bc)
    result.cells_only_native = sorted(set_native_bc - set_derived_bc)[:50]
    result.cells_only_derived = sorted(set_derived_bc - set_native_bc)[:50]

    result.n_features_native = len(set_native_feat)
    result.n_features_derived = len(set_derived_feat)
    result.n_features_overlap = len(set_native_feat & set_derived_feat)
    result.feature_jaccard = _jaccard(set_native_feat, set_derived_feat)
    result.features_only_native = sorted(set_native_feat - set_derived_feat)[:50]
    result.features_only_derived = sorted(set_derived_feat - set_native_feat)[:50]

    overlap_bc = sorted(set_native_bc & set_derived_bc)
    if not overlap_bc:
        return result

    bc_native_idx = {b: i for i, b in enumerate(native_barcodes)}
    bc_derived_idx = {b: i for i, b in enumerate(derived_barcodes)}

    if layer in gene_adata.layers:
        derived_mat = gene_adata.layers[layer].T.tocsr()
    else:
        derived_mat = gene_adata.X.T.tocsr()

    native_totals = []
    derived_totals = []
    for bc in overlap_bc:
        ni = bc_native_idx[bc]
        di = bc_derived_idx[bc]
        native_totals.append(float(native_mat[:, ni].sum()))
        derived_totals.append(float(derived_mat[:, di].sum()))

    native_totals = np.array(native_totals)
    derived_totals = np.array(derived_totals)
    if len(native_totals) > 1 and native_totals.std() > 0:
        result.total_count_correlation = float(
            np.corrcoef(native_totals, derived_totals)[0, 1]
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(native_totals - derived_totals) / np.maximum(native_totals, 1)
        result.mean_abs_rel_diff = float(np.nanmean(rel))

    return result
