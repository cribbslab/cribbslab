"""Differential isoform usage between two cell groups."""

from __future__ import annotations

from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats


def differential_isoform_usage(
    tx_adata: ad.AnnData,
    group_a: list[str],
    group_b: list[str],
    tx_to_gene: Optional[pd.DataFrame] = None,
    min_cells: int = 3,
) -> pd.DataFrame:
    """
    Compare isoform proportions between two barcode groups.

    Uses chi-squared test per gene with >=2 isoforms. Returns volcano-style table.
    """
    if tx_adata is None:
        return pd.DataFrame()

    cells_a = [c for c in group_a if c in tx_adata.obs_names]
    cells_b = [c for c in group_b if c in tx_adata.obs_names]
    if len(cells_a) < min_cells or len(cells_b) < min_cells:
        return pd.DataFrame()

    idx_a = [tx_adata.obs_names.get_loc(c) for c in cells_a]
    idx_b = [tx_adata.obs_names.get_loc(c) for c in cells_b]

    X = tx_adata.X
    if hasattr(X, "toarray"):
        # Avoid densifying full matrix; use sparse slices
        counts_a = X[idx_a, :].sum(axis=0)
        counts_b = X[idx_b, :].sum(axis=0)
        if hasattr(counts_a, "A1"):
            counts_a = counts_a.A1
            counts_b = counts_b.A1
        else:
            counts_a = np.asarray(counts_a).flatten()
            counts_b = np.asarray(counts_b).flatten()
    else:
        counts_a = X[idx_a, :].sum(axis=0)
        counts_b = X[idx_b, :].sum(axis=0)

    if tx_to_gene is None:
        tx_to_gene = pd.DataFrame(
            {"transcript_id": tx_adata.var_names, "gene_id": tx_adata.var_names}
        )
    else:
        tx_to_gene = tx_to_gene.rename(
            columns={"isoform_id": "transcript_id"}
            if "isoform_id" in tx_to_gene.columns
            else {}
        )

    gene_map = tx_to_gene.set_index("transcript_id")["gene_id"].to_dict()
    rows = []
    genes = set(gene_map.get(t, t) for t in tx_adata.var_names)

    for gene in genes:
        iso_indices = [
            i
            for i, tx in enumerate(tx_adata.var_names)
            if gene_map.get(tx, tx) == gene
        ]
        if len(iso_indices) < 2:
            continue
        ca = counts_a[iso_indices]
        cb = counts_b[iso_indices]
        if ca.sum() < 1 and cb.sum() < 1:
            continue
        table = np.vstack([ca, cb])
        try:
            chi2, p, _, _ = stats.chi2_contingency(table + 1e-9)
        except ValueError:
            continue
        frac_a = ca / (ca.sum() + 1e-9)
        frac_b = cb / (cb.sum() + 1e-9)
        effect = float(np.max(np.abs(frac_a - frac_b)))
        rows.append(
            {
                "gene_id": gene,
                "n_isoforms": len(iso_indices),
                "chi2": chi2,
                "pvalue": p,
                "effect_size": effect,
            }
        )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["neg_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    # Benjamini-Hochberg FDR
    p = df["pvalue"].values
    n = len(p)
    order = np.argsort(p)
    fdr = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = p[order[i]] * n / rank
        prev = min(prev, val)
        fdr[order[i]] = prev
    df["fdr"] = fdr
    return df.sort_values("pvalue")
