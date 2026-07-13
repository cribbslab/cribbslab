"""IsoQuant read-assignment QC."""

import streamlit as st

from isoquant_report.metrics.assignment_qc import (
    assignment_type_per_cell,
    assignment_type_summary,
    novel_known_summary,
    per_cell_assignment_rate,
    structural_category_summary,
)
from isoquant_report.plots import bar_counts, pie_proportions, scatter_xy, stacked_bar_per_cell
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()

st.header("Read assignment QC")
if data.assignments is None:
    st.warning("Read assignments not loaded.")
    st.stop()

at = assignment_type_summary(data.assignments)
st.plotly_chart(
    bar_counts(at, "assignment_type", "count", "Assignment types"),
    use_container_width=True,
)

sc = structural_category_summary(data.assignments)
if sc.empty:
    st.info("Structural categories unavailable in this IsoQuant version.")
else:
    st.plotly_chart(
        pie_proportions(sc, "category", "proportion", "Structural categories"),
        use_container_width=True,
    )

nk = novel_known_summary(data.assignments)
if not nk.empty:
    st.plotly_chart(
        pie_proportions(nk, "class", "proportion", "Novel vs known"),
        use_container_width=True,
    )

pc = assignment_type_per_cell(data.assignments)
st.plotly_chart(
    stacked_bar_per_cell(pc, "Assignment type per cell"),
    use_container_width=True,
)

rate = per_cell_assignment_rate(data.assignments)
st.plotly_chart(
    scatter_xy(rate, "assigned_reads", "assignment_rate", "Assignment rate"),
    use_container_width=True,
)

st.download_button(
    "Download assignments (CSV)",
    data.assignments.to_csv(index=False),
    file_name="assignments.csv",
)
