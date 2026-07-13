"""Splicing / isoform structural QC."""

import streamlit as st

from isoquant_report.metrics.saturation import isoform_saturation_curve, umi_saturation_curve
from isoquant_report.plots import histogram, saturation_line
from isoquant_report.state import get_report_data, sidebar_filters

sidebar_filters()
data = get_report_data()

st.header("Splicing / structural QC")

if data.gtf_models is not None and not data.gtf_models.empty:
    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(
            histogram(data.gtf_models, "length", "Transcript length"),
            use_container_width=True,
        )
    with c2:
        st.plotly_chart(
            histogram(data.gtf_models, "exon_count", "Exon count"),
            use_container_width=True,
        )
    mono = data.gtf_models[data.gtf_models["exon_count"] == 1]
    st.metric("Mono-exonic transcripts", len(mono))

if data.assignments is not None and "assignment_events" in data.assignments.columns:
    ir = data.assignments["assignment_events"].str.contains("intron_retention", na=False)
    st.metric("Reads with intron retention", int(ir.sum()))

if data.assignments is not None:
    sat = isoform_saturation_curve(data.assignments)
    st.plotly_chart(
        saturation_line(sat, "fraction_reads", "n_isoforms", "Isoform saturation"),
        use_container_width=True,
    )

if data.allinfo is not None:
    sat_u = umi_saturation_curve(data.allinfo)
    st.plotly_chart(
        saturation_line(sat_u, "fraction_umis", "n_genes", "UMI saturation"),
        use_container_width=True,
    )

if data.assignments is not None and "novel_known" in data.assignments.columns:
    nk = data.assignments["novel_known"].value_counts(normalize=True)
    st.bar_chart(nk)
