"""Saturation curve metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd


def isoform_saturation_curve(
    assignments: pd.DataFrame,
    n_steps: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Downsample reads and count distinct isoforms detected.

    Uses isoform_id column from assignments.
    """
    if "isoform_id" not in assignments.columns:
        return pd.DataFrame()
    df = assignments.dropna(subset=["isoform_id"])
    n = len(df)
    if n == 0:
        return pd.DataFrame()

    rng = np.random.default_rng(seed)
    indices = rng.permutation(n)
    fracs = np.linspace(0.05, 1.0, n_steps)
    rows = []
    for frac in fracs:
        k = max(1, int(n * frac))
        sub = df.iloc[indices[:k]]
        rows.append(
            {
                "fraction_reads": frac,
                "n_reads": k,
                "n_isoforms": sub["isoform_id"].nunique(),
            }
        )
    return pd.DataFrame(rows)


def umi_saturation_curve(allinfo: pd.DataFrame, n_steps: int = 20) -> pd.DataFrame:
    """Downsample UMIs and count distinct genes detected."""
    if allinfo.empty or "umi" not in allinfo.columns:
        return pd.DataFrame()
    df = allinfo.dropna(subset=["umi", "gene_id"]) if "gene_id" in allinfo.columns else allinfo
    n = len(df)
    if n == 0:
        return pd.DataFrame()
    fracs = np.linspace(0.05, 1.0, n_steps)
    rows = []
    for frac in fracs:
        k = max(1, int(n * frac))
        sub = df.iloc[:k]
        rows.append(
            {
                "fraction_umis": frac,
                "n_umis": k,
                "n_genes": sub["gene_id"].nunique() if "gene_id" in sub.columns else 0,
            }
        )
    return pd.DataFrame(rows)
