"""Load IsoQuant grouped count matrices."""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from isoquant_report.isoquant_io.discover import IsoQuantPaths


def _load_mtx(mtx_path: str) -> sp.csr_matrix:
    return mmread(mtx_path).tocsr()


def load_grouped_matrix(
    paths: IsoQuantPaths,
    feature: str = "gene",
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """
    Load a barcode-grouped count matrix from discovered IsoQuant paths.

    Returns (features x barcodes CSR matrix, feature ids, barcode ids).
    """
    if feature == "gene":
        mtx_path = paths.gene_grouped_mtx
        feat_path = paths.gene_grouped_features
        bc_path = paths.gene_grouped_barcodes
        linear_path = paths.gene_grouped_linear
    else:
        mtx_path = paths.transcript_grouped_mtx
        feat_path = paths.transcript_grouped_features
        bc_path = paths.transcript_grouped_barcodes
        linear_path = paths.transcript_grouped_linear

    if mtx_path and feat_path and bc_path:
        mat = _load_mtx(mtx_path)
        features = pd.read_csv(feat_path, sep="\t", header=None)[0].tolist()
        barcodes = pd.read_csv(bc_path, sep="\t", header=None)[0].tolist()
        return mat.tocsr(), features, barcodes

    if linear_path:
        df = pd.read_csv(
            linear_path,
            sep="\t",
            header=None,
            names=["feature_id", "group_id", "count"],
            comment="#",
        )
        features = sorted(df["feature_id"].unique())
        barcodes = sorted(df["group_id"].unique())
        feat_idx = {f: i for i, f in enumerate(features)}
        bc_idx = {b: i for i, b in enumerate(barcodes)}
        rows = df["feature_id"].map(feat_idx).values
        cols = df["group_id"].map(bc_idx).values
        data = df["count"].values.astype(np.float32)
        mat = sp.csr_matrix(
            (data, (rows, cols)),
            shape=(len(features), len(barcodes)),
        )
        return mat, features, barcodes

    raise FileNotFoundError(
        f"No grouped {feature} count matrix found in {paths.isoquant_dir}"
    )


def load_grouped_tpm(
    paths: IsoQuantPaths,
    feature: str = "gene",
) -> Optional[tuple[sp.csr_matrix, list[str], list[str]]]:
    """Load grouped TPM matrix if present."""
    mtx_path = (
        paths.gene_grouped_tpm_mtx
        if feature == "gene"
        else paths.transcript_grouped_tpm_mtx
    )
    if not mtx_path:
        return None
    feat_path = mtx_path.replace(".matrix.mtx", ".features.tsv")
    bc_path = mtx_path.replace(".matrix.mtx", ".barcodes.tsv")
    if not (feat_path and bc_path):
        return None
    mat = _load_mtx(mtx_path)
    features = pd.read_csv(feat_path, sep="\t", header=None)[0].tolist()
    barcodes = pd.read_csv(bc_path, sep="\t", header=None)[0].tolist()
    return mat.tocsr(), features, barcodes
