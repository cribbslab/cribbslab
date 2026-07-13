"""Compute QC metrics from assignments and AnnData."""

from __future__ import annotations

import pandas as pd


def assignment_type_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    """Count reads per assignment_type."""
    if "assignment_type" not in assignments.columns:
        return pd.DataFrame()
    vc = assignments["assignment_type"].value_counts().reset_index()
    vc.columns = ["assignment_type", "count"]
    return vc


def assignment_type_per_cell(assignments: pd.DataFrame) -> pd.DataFrame:
    """Assignment type counts per barcode."""
    if "barcode" not in assignments.columns:
        return pd.DataFrame()
    df = assignments.dropna(subset=["barcode"])
    return (
        df.groupby(["barcode", "assignment_type"])
        .size()
        .reset_index(name="count")
    )


def novel_known_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    """Novel vs known transcript model proportions."""
    if "novel_known" not in assignments.columns:
        return pd.DataFrame()
    vc = assignments["novel_known"].value_counts(normalize=True).reset_index()
    vc.columns = ["class", "proportion"]
    return vc


def structural_category_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    """SQANTI-style category proportions."""
    if "structural_category" not in assignments.columns:
        return pd.DataFrame()
    df = assignments.dropna(subset=["structural_category"])
    if df.empty:
        return pd.DataFrame()
    vc = df["structural_category"].value_counts(normalize=True).reset_index()
    vc.columns = ["category", "proportion"]
    return vc


def per_cell_assignment_rate(assignments: pd.DataFrame) -> pd.DataFrame:
    """Per-cell assigned read counts and assignment rate."""
    if "barcode" not in assignments.columns:
        return pd.DataFrame()
    df = assignments.dropna(subset=["barcode"])
    total = df.groupby("barcode").size().reset_index(name="assigned_reads")
    if "assignment_type" in df.columns:
        assigned = df[
            ~df["assignment_type"].isin(["intergenic", "noninformative", "ambiguous"])
        ]
        ok = assigned.groupby("barcode").size().reset_index(name="informative_reads")
        total = total.merge(ok, on="barcode", how="left").fillna(0)
        total["assignment_rate"] = (
            total["informative_reads"] / total["assigned_reads"]
        ).fillna(0)
    return total
