"""Counts matrix QC."""

import streamlit as st

from isoquant_report.metrics.cell_qc import apply_qc_filters, knee_plot_data
from isoquant_report.plots import knee_plot, scatter_xy, violin_metric
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()
cfg = st.session_state["config"]

st.header("Counts matrix QC")
if data.gene_qc_adata is None:
    st.warning("Gene AnnData not loaded.")
    st.stop()

mask = apply_qc_filters(
    data.gene_qc_adata,
    cfg.min_counts,
    cfg.min_genes,
    cfg.max_mito,
    cfg.max_unspliced_fraction,
    data.barcode_qc,
)
n_pass = int(mask.sum())
n_fail = int((~mask).sum())
st.metric("Cells passing QC", f"{n_pass} / {n_pass + n_fail}")

obs = data.gene_qc_adata.obs.copy()
obs["pass_qc"] = mask

c1, c2 = st.columns(2)
with c1:
    st.plotly_chart(
        violin_metric(obs.reset_index(), "total_counts", "Total counts"),
        use_container_width=True,
    )
with c2:
    st.plotly_chart(
        violin_metric(obs.reset_index(), "n_genes_by_counts", "Genes detected"),
        use_container_width=True,
    )

st.plotly_chart(
    violin_metric(obs.reset_index(), "pct_counts_mt", "Mitochondrial %"),
    use_container_width=True,
)
st.plotly_chart(knee_plot(knee_plot_data(data.gene_qc_adata)), use_container_width=True)
st.plotly_chart(
    scatter_xy(
        obs.reset_index(),
        "total_counts",
        "n_genes_by_counts",
        "Depth vs genes",
        color="pass_qc",
    ),
    use_container_width=True,
)

if data.transcript_adata is not None:
    tx_obs = data.transcript_adata.obs.copy()
    st.plotly_chart(
        violin_metric(
            tx_obs.reset_index().assign(
                isoforms_detected=(data.transcript_adata.X > 0).sum(axis=1).A1
                if hasattr(data.transcript_adata.X, "A1")
                else (data.transcript_adata.X > 0).sum(axis=1)
            ),
            "isoforms_detected",
            "Isoforms detected",
        ),
        use_container_width=True,
    )
