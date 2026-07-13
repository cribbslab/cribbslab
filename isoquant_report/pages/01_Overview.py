"""Overview / summary page."""

import streamlit as st

from isoquant_report.export_html import export_html
from isoquant_report.report import headline_metrics
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()
metrics = headline_metrics(data)

st.header("Overview")
cols = st.columns(4)
labels = [
    ("Cells", metrics.get("n_cells", "—")),
    ("Genes", metrics.get("n_genes", "—")),
    ("Isoforms", metrics.get("n_isoforms", "—")),
    ("Median UMIs", metrics.get("median_umis", "—")),
]
for col, (label, val) in zip(cols, labels):
    col.metric(label, val)

st.subheader("Data completeness")
comp = data.paths.completeness()
st.table({k: "yes" if v else "no" for k, v in comp.items()})

if data.reconcile:
    st.subheader("Reconciliation (native vs derived)")
    r = data.reconcile
    st.write(
        f"Cell Jaccard: **{r.cell_jaccard:.3f}** | "
        f"Feature Jaccard: **{r.feature_jaccard:.3f}** | "
        f"Count correlation: **{r.total_count_correlation:.3f}**"
    )
    if r.cells_only_native:
        st.caption(f"Cells only in native (first 10): {r.cells_only_native[:10]}")
    if r.cells_only_derived:
        st.caption(f"Cells only in derived (first 10): {r.cells_only_derived[:10]}")

if data.errors:
    st.warning("Load warnings: " + "; ".join(data.errors))

if st.button("Export static HTML"):
    out = f"report_export_{data.config.sample}.html"
    export_html(data, out)
    st.success(f"Wrote {out}")
