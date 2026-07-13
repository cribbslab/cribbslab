"""
IsoQuant single-cell long-read QC report.

Run from the isoquant_report directory:

    streamlit run app.py

Or with a pipeline-written config:

    streamlit run app.py -- --config isoquant_report/SAMPLE/config.yaml
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

from isoquant_report.config import ReportConfig
from isoquant_report.state import init_state, sidebar_filters

st.set_page_config(
    page_title="IsoQuant QC Report",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Optional --config from CLI
config_path = None
if "--config" in sys.argv:
    idx = sys.argv.index("--config")
    if idx + 1 < len(sys.argv):
        config_path = sys.argv[idx + 1]

if config_path and Path(config_path).exists():
    config = ReportConfig.from_yaml(config_path)
else:
    config = ReportConfig()

init_state(config)
sidebar_filters()

st.title("IsoQuant single-cell long-read QC")
st.markdown(
    "Use the sidebar to set paths and QC thresholds. "
    "Navigate pages via the left menu."
)
st.info(
    "Pages: Overview, Read assignment QC, Counts matrix QC, "
    "Gene vs isoform, Splicing structural QC, Cluster exploration."
)
