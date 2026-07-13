"""Gene vs isoform comparison."""

import numpy as np
import streamlit as st

from isoquant_report.isoquant_io.counts import load_grouped_matrix
from isoquant_report.metrics.isoform_stats import gene_isoform_detection, isoforms_per_gene
from isoquant_report.plots import histogram, scatter_xy
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()

st.header("Gene vs isoform comparison")
if data.gene_adata is None or data.transcript_adata is None:
    st.warning("Both gene and transcript h5ad files are required.")
    st.stop()

det = gene_isoform_detection(data.gene_adata, data.transcript_adata)
st.plotly_chart(
    scatter_xy(det, "genes_detected", "isoforms_detected", "Detection per cell"),
    use_container_width=True,
)

ipg = isoforms_per_gene(data.assignments, data.gtf_models) if data.assignments is not None else isoforms_per_gene(None, data.gtf_models)
if not ipg.empty:
    st.plotly_chart(
        histogram(ipg, "n_isoforms", "Isoforms per gene"),
        use_container_width=True,
    )
    st.dataframe(
        ipg.sort_values("n_isoforms", ascending=False).head(100),
        use_container_width=True,
    )
    st.download_button(
        "Download isoform diversity table",
        ipg.to_csv(index=False),
        file_name="isoforms_per_gene.csv",
    )

# TPM / count agreement
try:
    native_mat, native_features, native_barcodes = load_grouped_matrix(
        data.paths, "gene"
    )
    overlap = sorted(set(native_barcodes) & set(data.gene_adata.obs_names))
    if overlap:
        bc_n = {b: i for i, b in enumerate(native_barcodes)}
        bc_d = {b: i for i, b in enumerate(data.gene_adata.obs_names)}
        layer = (
            data.gene_adata.layers["total"]
            if "total" in data.gene_adata.layers
            else data.gene_adata.X
        )
        derived_mat = layer.T if hasattr(layer, "T") else layer
        native_totals = [float(native_mat[:, bc_n[b]].sum()) for b in overlap]
        derived_totals = [float(derived_mat[:, bc_d[b]].sum()) for b in overlap]
        agree = {
            "barcode": overlap,
            "native_total": native_totals,
            "derived_total": derived_totals,
        }
        import pandas as pd
        adf = pd.DataFrame(agree)
        st.plotly_chart(
            scatter_xy(adf, "native_total", "derived_total", "Count agreement"),
            use_container_width=True,
        )
        if len(native_totals) > 1:
            corr = np.corrcoef(native_totals, derived_totals)[0, 1]
            st.metric("Pearson correlation", f"{corr:.3f}")
except Exception as exc:
    st.info(f"Native grouped matrix agreement unavailable: {exc}")
