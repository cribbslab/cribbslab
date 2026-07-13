"""Plotly figure builders for IsoQuant QC report."""

from __future__ import annotations

from typing import Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def bar_counts(df: pd.DataFrame, x: str, y: str, title: str) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.bar(df, x=x, y=y, title=title)
    fig.update_layout(template="plotly_white")
    return fig


def pie_proportions(df: pd.DataFrame, names: str, values: str, title: str) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.pie(df, names=names, values=values, title=title)
    fig.update_layout(template="plotly_white")
    return fig


def stacked_bar_per_cell(
    df: pd.DataFrame, title: str, max_cells: int = 50
) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    top = df.groupby("barcode")["count"].sum().nlargest(max_cells).index
    sub = df[df["barcode"].isin(top)]
    fig = px.bar(
        sub,
        x="barcode",
        y="count",
        color="assignment_type",
        title=title,
        barmode="stack",
    )
    fig.update_layout(template="plotly_white", xaxis_tickangle=-45)
    return fig


def scatter_xy(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    color: Optional[str] = None,
) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.scatter(df, x=x, y=y, color=color, title=title, opacity=0.6)
    fig.update_layout(template="plotly_white")
    return fig


def violin_metric(df: pd.DataFrame, y: str, title: str) -> go.Figure:
    if df.empty or y not in df.columns:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.violin(df, y=y, box=True, points="outliers", title=title)
    fig.update_layout(template="plotly_white")
    return fig


def knee_plot(df: pd.DataFrame, title: str = "Barcode rank plot") -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.line(df, x="rank", y="counts", title=title, log_x=True, log_y=True)
    fig.update_layout(template="plotly_white")
    return fig


def saturation_line(df: pd.DataFrame, x: str, y: str, title: str) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.line(df, x=x, y=y, title=title, markers=True)
    fig.update_layout(template="plotly_white")
    return fig


def volcano(df: pd.DataFrame, title: str = "Differential isoform usage") -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.scatter(
        df,
        x="effect_size",
        y="neg_log10_p",
        hover_data=["gene_id", "fdr"],
        title=title,
        labels={"effect_size": "Max isoform fraction delta", "neg_log10_p": "-log10(p)"},
    )
    fig.update_layout(template="plotly_white")
    return fig


def histogram(df: pd.DataFrame, x: str, title: str) -> go.Figure:
    if df.empty:
        return go.Figure().update_layout(title=title + " (no data)")
    fig = px.histogram(df, x=x, title=title)
    fig.update_layout(template="plotly_white")
    return fig
