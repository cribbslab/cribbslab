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
        nk_gtf = _novelty_from_gtf(data.gtf_models)
        if not nk_gtf.empty and nk_gtf["count"].sum() > 0:
            figs["model_novelty"] = pie_proportions(
                nk_gtf, "class", "count", "Transcript models: novel vs known"
            )

    if data.allinfo is not None:
        sat_u = umi_saturation_curve(data.allinfo)
        figs["umi_saturation"] = saturation_line(
            sat_u, "fraction_umis", "n_genes", "UMI saturation"
        )

    # ---- Splice / unspliced QC ------------------------------------------
    if data.barcode_qc is not None and not data.barcode_qc.empty:
        bq = data.barcode_qc
        if "unspliced_fraction" in bq.columns:
            figs["unspliced_fraction_hist"] = histogram(
                bq, "unspliced_fraction", "Unspliced fraction per cell"
            )
        if {"spliced_umis", "unspliced_umis"} <= set(bq.columns):
            figs["spliced_vs_unspliced"] = scatter_xy(
                bq,
                "spliced_umis",
                "unspliced_umis",
                "Spliced vs unspliced UMIs per cell",
            )
        comp_cols = [
            c
            for c in ("spliced_umis", "unspliced_umis", "ambiguous_umis")
            if c in bq.columns
        ]
        if comp_cols:
            comp = pd.DataFrame(
                {
                    "class": [c.replace("_umis", "") for c in comp_cols],
                    "umis": [float(bq[c].sum()) for c in comp_cols],
                }
            )
            figs["splice_composition"] = pie_proportions(
                comp, "class", "umis", "Global splice composition (UMIs)"
            )

    # ---- Top expressed features -----------------------------------------
    if data.gene_adata is not None:
        tg = _top_features(data.gene_adata, n=20)
        if not tg.empty:
            figs["top_genes"] = bar_counts(
                tg, "feature", "total_counts", "Top 20 genes by total UMIs"
            )
    if data.transcript_adata is not None:
        ti = _top_features(data.transcript_adata, n=20)
        if not ti.empty:
            figs["top_isoforms"] = bar_counts(
                ti, "feature", "total_counts", "Top 20 isoforms by total UMIs"
            )

    return figs


def _top_features(adata: ad.AnnData, n: int = 20) -> pd.DataFrame:
    """Return the top-n features by summed counts across cells."""
    try:
        X = adata.X
        totals = np.asarray(X.sum(axis=0)).flatten()
    except Exception:
        return pd.DataFrame(columns=["feature", "total_counts"])
    if totals.size == 0:
        return pd.DataFrame(columns=["feature", "total_counts"])
    order = np.argsort(totals)[::-1][:n]
    return pd.DataFrame(
        {
            "feature": np.asarray(adata.var_names)[order],
            "total_counts": totals[order],
        }
    )


def _novelty_from_gtf(gtf_models: pd.DataFrame) -> pd.DataFrame:
    """Classify transcript models as known (Ensembl) vs novel by ID prefix."""
    if gtf_models is None or gtf_models.empty:
        return pd.DataFrame(columns=["class", "count"])
    tid = gtf_models["transcript_id"].astype(str)
    is_known = tid.str.upper().str.startswith("ENS")
    n_known = int(is_known.sum())
    n_novel = int((~is_known).sum())
    return pd.DataFrame(
        {"class": ["known", "novel"], "count": [n_known, n_novel]}
    )


def _fmt(value: Any) -> Any:
    """Format numbers for display tables."""
    if isinstance(value, float):
        if np.isnan(value):
            return "n/a"
        if abs(value) >= 1000:
            return f"{value:,.0f}"
        return f"{value:,.3g}"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    return value


def build_summary_tables(data: ReportData) -> dict[str, pd.DataFrame]:
    """Build human-readable summary tables to accompany the figures."""
    tables: dict[str, pd.DataFrame] = {}
    cfg = data.config

    # --- Assignment types ------------------------------------------------
    if data.assignments is not None:
        try:
            at = assignment_type_summary(data.assignments)
            if not at.empty:
                at = at.copy()
                at["percent"] = (100 * at["count"] / at["count"].sum()).round(2)
                tables["Assignment types"] = at
        except Exception as exc:
            data.errors.append(f"table assignment_types: {exc}")
        try:
            sc = structural_category_summary(data.assignments)
            if not sc.empty:
                sc = sc.copy()
                sc["percent"] = (100 * sc["proportion"]).round(2)
                tables["Structural categories"] = sc[["category", "percent"]]
        except Exception:
            pass

    # --- Cell QC summary -------------------------------------------------
    if data.gene_qc_adata is not None:
        try:
            obs = data.gene_qc_adata.obs
            mask = apply_qc_filters(
                data.gene_qc_adata,
                cfg.min_counts,
                cfg.min_genes,
                cfg.max_mito,
                cfg.max_unspliced_fraction,
                data.barcode_qc,
            )
            n_total = int(len(mask))
            n_pass = int(mask.sum())
            rows = [
                ("Cells (called barcodes)", _fmt(n_total)),
                (
                    "Cells passing QC filters",
                    f"{n_pass:,} ({100 * n_pass / max(n_total, 1):.1f}%)",
                ),
                ("Cells failing QC filters", _fmt(n_total - n_pass)),
            ]
            if "total_counts" in obs:
                rows.append(("Median UMIs / cell", _fmt(float(obs["total_counts"].median()))))
            if "n_genes_by_counts" in obs:
                rows.append(("Median genes / cell", _fmt(float(obs["n_genes_by_counts"].median()))))
            if "pct_counts_mt" in obs:
                rows.append(("Median % mito", _fmt(float(obs["pct_counts_mt"].median()))))
            if "pct_counts_ribo" in obs:
                rows.append(("Median % ribo", _fmt(float(obs["pct_counts_ribo"].median()))))
            rows += [
                ("Threshold: min UMIs", _fmt(cfg.min_counts)),
                ("Threshold: min genes", _fmt(cfg.min_genes)),
                ("Threshold: max mito fraction", _fmt(cfg.max_mito)),
                ("Threshold: max unspliced fraction", _fmt(cfg.max_unspliced_fraction)),
            ]
            tables["Cell QC summary"] = pd.DataFrame(rows, columns=["metric", "value"])
        except Exception as exc:
            data.errors.append(f"table cell_qc: {exc}")

    # --- Splice summary --------------------------------------------------
    if data.barcode_qc is not None and not data.barcode_qc.empty:
        try:
            bq = data.barcode_qc
            rows = []
            for col, label in [
                ("spliced_umis", "Total spliced UMIs"),
                ("unspliced_umis", "Total unspliced UMIs"),
                ("ambiguous_umis", "Total ambiguous UMIs"),
                ("total_umis", "Total UMIs"),
            ]:
                if col in bq.columns:
                    rows.append((label, _fmt(float(bq[col].sum()))))
            if "unspliced_fraction" in bq.columns:
                rows.append(
                    ("Median unspliced fraction", _fmt(float(bq["unspliced_fraction"].median())))
                )
                rows.append(
                    ("Mean unspliced fraction", _fmt(float(bq["unspliced_fraction"].mean())))
                )
            if rows:
                tables["Splice summary"] = pd.DataFrame(rows, columns=["metric", "value"])
        except Exception as exc:
            data.errors.append(f"table splice: {exc}")

    # --- Isoform structure ----------------------------------------------
    if data.gtf_models is not None and not data.gtf_models.empty:
        try:
            gm = data.gtf_models
            nk = _novelty_from_gtf(gm)
            if not nk.empty and nk["count"].sum() > 0:
                nk = nk.copy()
                nk["percent"] = (100 * nk["count"] / nk["count"].sum()).round(2)
                tables["Transcript models (novel vs known)"] = nk
            rows = [
                ("Transcript models", _fmt(int(gm["transcript_id"].nunique()))),
                ("Genes with models", _fmt(int(gm["gene_id"].nunique()))),
                ("Median transcript length (bp)", _fmt(float(gm["length"].median()))),
                ("Median exons / transcript", _fmt(float(gm["exon_count"].median()))),
                ("Mono-exonic transcripts", _fmt(int((gm["exon_count"] == 1).sum()))),
            ]
            tables["Isoform structure"] = pd.DataFrame(rows, columns=["metric", "value"])
        except Exception as exc:
            data.errors.append(f"table isoform_structure: {exc}")

    # --- Reconciliation --------------------------------------------------
    if data.reconcile is not None:
        try:
            r = data.reconcile
            rows = [
                ("Cells (native / derived)", f"{r.n_cells_native:,} / {r.n_cells_derived:,}"),
                ("Cell overlap", _fmt(r.n_cells_overlap)),
                ("Cell Jaccard", _fmt(round(r.cell_jaccard, 3))),
                ("Features (native / derived)", f"{r.n_features_native:,} / {r.n_features_derived:,}"),
                ("Feature Jaccard", _fmt(round(r.feature_jaccard, 3))),
                ("Total-count correlation", _fmt(round(r.total_count_correlation, 3))),
                ("Mean abs. relative diff", _fmt(round(r.mean_abs_rel_diff, 3))),
            ]
            tables["Matrix reconciliation"] = pd.DataFrame(rows, columns=["metric", "value"])
        except Exception as exc:
            data.errors.append(f"table reconcile: {exc}")

    return tables


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
