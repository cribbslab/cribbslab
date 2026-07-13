"""Streamlit session state helpers."""

from __future__ import annotations

from typing import Optional

import streamlit as st

from isoquant_report.config import ReportConfig
from isoquant_report.report import ReportData, load_report_data


def init_state(config: ReportConfig) -> None:
    """Initialise session state keys if absent."""
    defaults = {
        "config": config,
        "min_counts": config.min_counts,
        "min_genes": config.min_genes,
        "max_mito": config.max_mito,
        "max_unspliced_fraction": config.max_unspliced_fraction,
        "cluster_col": None,
        "selected_cluster": None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


@st.cache_data(show_spinner="Loading IsoQuant data...")
def cached_report_data(config_dict: dict) -> ReportData:
    """Cache report data load by config dict."""
    config = ReportConfig.from_dict(config_dict)
    return load_report_data(config)


def get_report_data() -> ReportData:
    """Return cached report data from session config."""
    cfg: ReportConfig = st.session_state["config"]
    return cached_report_data(cfg.to_dict())


def sidebar_filters() -> ReportConfig:
    """Render sidebar path overrides and QC sliders; return updated config."""
    st.sidebar.header("Configuration")
    cfg: ReportConfig = st.session_state.get("config", ReportConfig())

    cfg.isoquant_dir = st.sidebar.text_input(
        "IsoQuant directory", value=cfg.isoquant_dir
    )
    cfg.sample = st.sidebar.text_input("Sample", value=cfg.sample)
    cfg.gene_h5ad = st.sidebar.text_input("Gene h5ad", value=cfg.gene_h5ad)
    cfg.transcript_h5ad = st.sidebar.text_input(
        "Transcript h5ad", value=cfg.transcript_h5ad
    )
    cfg.barcode_qc = st.sidebar.text_input(
        "Barcode QC TSV (optional)", value=cfg.barcode_qc
    )

    st.sidebar.header("QC thresholds")
    st.session_state["min_counts"] = st.sidebar.slider(
        "Min counts", 0, 10000, st.session_state.get("min_counts", 500)
    )
    st.session_state["min_genes"] = st.sidebar.slider(
        "Min genes", 0, 5000, st.session_state.get("min_genes", 200)
    )
    st.session_state["max_mito"] = st.sidebar.slider(
        "Max mito fraction", 0.0, 1.0, st.session_state.get("max_mito", 0.2)
    )
    st.session_state["max_unspliced_fraction"] = st.sidebar.slider(
        "Max unspliced fraction",
        0.0,
        1.0,
        st.session_state.get("max_unspliced_fraction", 0.95),
    )

    cfg.min_counts = st.session_state["min_counts"]
    cfg.min_genes = st.session_state["min_genes"]
    cfg.max_mito = st.session_state["max_mito"]
    cfg.max_unspliced_fraction = st.session_state["max_unspliced_fraction"]
    st.session_state["config"] = cfg
    return cfg
