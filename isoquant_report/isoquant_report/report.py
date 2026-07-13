"""Core report data loading and figure assembly."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd

from isoquant_report.config import ReportConfig
from isoquant_report.isoquant_io.allinfo import load_allinfo_barcodes
from isoquant_report.isoquant_io.assignments import load_assignments
from isoquant_report.isoquant_io.counts import load_grouped_matrix, load_grouped_tpm
from isoquant_report.isoquant_io.discover import IsoQuantPaths, discover_isoquant
from isoquant_report.isoquant_io.gtf import parse_transcript_models
from isoquant_report.matrices import load_matrices
from isoquant_report.metrics.assignment_qc import (
    assignment_type_per_cell,
    assignment_type_summary,
    novel_known_summary,
    per_cell_assignment_rate,
    structural_category_summary,
)
from isoquant_report.metrics.cell_qc import (
    apply_qc_filters,
    compute_qc_metrics,
    knee_plot_data,
)
from isoquant_report.metrics.differential_isoform import differential_isoform_usage
from isoquant_report.metrics.isoform_stats import (
    gene_isoform_detection,
    isoforms_per_gene,
)
from isoquant_report.metrics.saturation import isoform_saturation_curve, umi_saturation_curve
from isoquant_report.plots import (
    bar_counts,
    histogram,
    knee_plot,
    pie_proportions,
    saturation_line,
    scatter_xy,
    stacked_bar_per_cell,
    violin_metric,
    volcano,
)
from isoquant_report.reconcile import ReconcileResult, reconcile_gene_matrix
import plotly.graph_objects as go


@dataclass
class ReportData:
    """All loaded artefacts and precomputed tables for one sample."""

    config: ReportConfig
    paths: IsoQuantPaths
    assignments: Optional[pd.DataFrame] = None
    gtf_models: Optional[pd.DataFrame] = None
    allinfo: Optional[pd.DataFrame] = None
    gene_adata: Optional[ad.AnnData] = None
    transcript_adata: Optional[ad.AnnData] = None
    barcode_qc: Optional[pd.DataFrame] = None
    reconcile: Optional[ReconcileResult] = None
    gene_qc_adata: Optional[ad.AnnData] = None
    errors: list[str] = field(default_factory=list)


def load_report_data(config: ReportConfig) -> ReportData:
    """Load all inputs referenced in config, tolerating missing files."""
    data = ReportData(config=config, paths=discover_isoquant(
        config.isoquant_dir, config.sample
    ))

    try:
        data.assignments = load_assignments(data.paths)
    except Exception as exc:
        data.errors.append(f"assignments: {exc}")

    if data.paths.allinfo:
        try:
            data.allinfo = load_allinfo_barcodes(data.paths.allinfo)
        except Exception as exc:
            data.errors.append(f"allinfo: {exc}")

    gtf_path = data.paths.transcript_models_gtf or data.paths.extended_annotation_gtf
    if gtf_path:
        try:
            data.gtf_models = parse_transcript_models(gtf_path)
        except Exception as exc:
            data.errors.append(f"gtf: {exc}")

    mats = load_matrices(config)
    data.gene_adata = mats["gene_adata"]
    data.transcript_adata = mats["transcript_adata"]
    data.barcode_qc = mats["barcode_qc"]

    if data.gene_adata is not None:
        try:
            data.gene_qc_adata = compute_qc_metrics(
                data.gene_adata,
                mito_prefix=config.mito_prefix,
                ribo_prefixes=config.ribo_prefixes,
            )
            data.reconcile = reconcile_gene_matrix(
                data.paths, data.gene_adata, layer="total"
            )
        except Exception as exc:
            data.errors.append(f"reconcile: {exc}")

    return data


def build_all_figures(data: ReportData) -> dict[str, go.Figure]:
    """Build Plotly figures for static HTML export."""
    figs: dict[str, go.Figure] = {}
    cfg = data.config

    if data.assignments is not None:
        at = assignment_type_summary(data.assignments)
        figs["assignment_type"] = bar_counts(
            at.rename(columns={"assignment_type": "assignment_type"}),
            "assignment_type",
            "count",
            "Assignment type distribution",
        )
        sc = structural_category_summary(data.assignments)
        if not sc.empty:
            figs["structural_category"] = pie_proportions(
                sc, "category", "proportion", "Structural categories"
            )
        nk = novel_known_summary(data.assignments)
        if not nk.empty:
            figs["novel_known"] = pie_proportions(
                nk, "class", "proportion", "Novel vs known"
            )
        pc = assignment_type_per_cell(data.assignments)
        figs["assignment_per_cell"] = stacked_bar_per_cell(
            pc, "Assignment type per cell (top 50)"
        )
        rate = per_cell_assignment_rate(data.assignments)
        figs["assignment_rate"] = scatter_xy(
            rate, "assigned_reads", "assignment_rate", "Per-cell assignment rate"
        )
        sat = isoform_saturation_curve(data.assignments)
        figs["isoform_saturation"] = saturation_line(
            sat, "fraction_reads", "n_isoforms", "Isoform detection saturation"
        )

    if data.gene_qc_adata is not None:
        obs = data.gene_qc_adata.obs.copy()
        mask = apply_qc_filters(
            data.gene_qc_adata,
            cfg.min_counts,
            cfg.min_genes,
            cfg.max_mito,
            cfg.max_unspliced_fraction,
            data.barcode_qc,
        )
        obs["pass_qc"] = mask
        figs["total_counts_violin"] = violin_metric(
            obs.reset_index(), "total_counts", "Total counts per cell"
        )
        figs["genes_detected_violin"] = violin_metric(
            obs.reset_index(), "n_genes_by_counts", "Genes detected per cell"
        )
        figs["mito_violin"] = violin_metric(
            obs.reset_index(), "pct_counts_mt", "Mitochondrial fraction"
        )
        knee = knee_plot_data(data.gene_qc_adata)
        figs["knee_plot"] = knee_plot(knee)
        figs["depth_vs_genes"] = scatter_xy(
            obs.reset_index(), "total_counts", "n_genes_by_counts", "Depth vs genes"
        )

    if data.gene_adata is not None and data.transcript_adata is not None:
        det = gene_isoform_detection(data.gene_adata, data.transcript_adata)
        figs["gene_vs_isoform_detection"] = scatter_xy(
            det, "genes_detected", "isoforms_detected", "Genes vs isoforms detected"
        )

    if data.assignments is not None:
        ipg = isoforms_per_gene(data.assignments, data.gtf_models)
        if not ipg.empty:
            figs["isoforms_per_gene"] = histogram(
                ipg, "n_isoforms", "Isoforms per gene"
            )

    if data.gtf_models is not None and not data.gtf_models.empty:
        figs["tx_length"] = histogram(
            data.gtf_models, "length", "Transcript model length"
        )
        figs["exon_count"] = histogram(
            data.gtf_models, "exon_count", "Exon count per transcript"
        )

    if data.allinfo is not None:
        sat_u = umi_saturation_curve(data.allinfo)
        figs["umi_saturation"] = saturation_line(
            sat_u, "fraction_umis", "n_genes", "UMI saturation"
        )

    return figs


def headline_metrics(data: ReportData) -> dict[str, Any]:
    """Summary metrics for overview cards."""
    m: dict[str, Any] = {"sample": data.config.sample}
    if data.gene_adata is not None:
        m["n_cells"] = data.gene_adata.n_obs
        m["n_genes"] = data.gene_adata.n_vars
    if data.transcript_adata is not None:
        m["n_isoforms"] = data.transcript_adata.n_vars
    if data.barcode_qc is not None and not data.barcode_qc.empty:
        m["median_umis"] = float(data.barcode_qc["total_umis"].median())
        m["median_unspliced_fraction"] = float(
            data.barcode_qc["unspliced_fraction"].median()
        )
    if data.assignments is not None and "assignment_type" in data.assignments.columns:
        unique_frac = (
            data.assignments["assignment_type"] == "unique"
        ).mean()
        m["pct_unique_assignments"] = round(100 * unique_frac, 1)
    if data.reconcile is not None:
        m["cell_jaccard"] = round(data.reconcile.cell_jaccard, 3)
        m["count_correlation"] = round(data.reconcile.total_count_correlation, 3)
    return m
