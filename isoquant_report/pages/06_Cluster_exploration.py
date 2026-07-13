"""Per-cell / per-cluster exploration and differential isoform usage."""

import streamlit as st

from isoquant_report.metrics.differential_isoform import differential_isoform_usage
from isoquant_report.plots import volcano
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()

st.header("Cluster exploration")

gene_adata = data.gene_adata
tx_adata = data.transcript_adata
if gene_adata is None:
    st.warning("Gene AnnData required.")
    st.stop()

# Cluster column selector
obs_cols = [c for c in gene_adata.obs.columns if gene_adata.obs[c].dtype.name in ("category", "object", "int64", "str")]
cluster_col = st.selectbox(
    "Cluster / cell-type column (optional)",
    options=["— none —"] + obs_cols,
)

if cluster_col != "— none —":
    groups = sorted(gene_adata.obs[cluster_col].astype(str).unique())
    selected = st.selectbox("Filter to group", options=["All"] + groups)
    if selected != "All":
        cells = gene_adata.obs_names[gene_adata.obs[cluster_col].astype(str) == selected]
        st.write(f"Showing {len(cells)} cells in **{selected}**")

st.subheader("Differential isoform usage")
if tx_adata is None:
    st.warning("Transcript AnnData required for DIU.")
    st.stop()

if cluster_col == "— none —":
    st.info(
        "No cluster labels in AnnData. Add a column to obs (e.g. leiden) "
        "or use QC pass/fail groups."
    )
    if data.barcode_qc is not None and "outlier" in data.barcode_qc.columns:
        group_a = data.barcode_qc[data.barcode_qc["outlier"] == 0]["barcode"].tolist()
        group_b = data.barcode_qc[data.barcode_qc["outlier"] == 1]["barcode"].tolist()
        st.caption("Using non-outlier vs outlier cells from splice QC.")
    else:
        st.stop()
else:
    g1 = st.selectbox("Group A", options=groups)
    g2 = st.selectbox("Group B", options=[g for g in groups if g != g1])
    group_a = gene_adata.obs_names[gene_adata.obs[cluster_col].astype(str) == g1].tolist()
    group_b = gene_adata.obs_names[gene_adata.obs[cluster_col].astype(str) == g2].tolist()

tx_to_gene = data.gtf_models[["transcript_id", "gene_id"]] if data.gtf_models is not None and not data.gtf_models.empty else None
diu = differential_isoform_usage(tx_adata, group_a, group_b, tx_to_gene)

if diu.empty:
    st.warning("No differential isoform results (need genes with >=2 isoforms).")
else:
    st.plotly_chart(volcano(diu), use_container_width=True)
    st.dataframe(diu.head(200), use_container_width=True)
    st.download_button(
        "Download DIU results",
        diu.to_csv(index=False),
        file_name="differential_isoform_usage.csv",
    )
