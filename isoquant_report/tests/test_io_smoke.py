"""Smoke tests for IsoQuant report package."""

from pathlib import Path

import pytest
import pandas as pd

from isoquant_report.config import ReportConfig
from isoquant_report.export_html import export_html
from isoquant_report.isoquant_io.discover import discover_isoquant
from isoquant_report.report import build_all_figures, headline_metrics, load_report_data

from tests.make_fixtures import SAMPLE, _write_fixtures


@pytest.fixture(scope="module")
def fixture_root(tmp_path_factory):
    root = _write_fixtures()
    return root


def test_discover(fixture_root):
    iq_dir = fixture_root / "isoquant" / SAMPLE / SAMPLE
    paths = discover_isoquant(str(iq_dir), SAMPLE)
    comp = paths.completeness()
    assert comp["gene_grouped_counts"]
    assert comp["read_assignments"]


def test_load_and_figures(fixture_root, tmp_path):
    iq_dir = fixture_root / "isoquant" / SAMPLE / SAMPLE
    splice = fixture_root / "splice_matrices"
    config = ReportConfig(
        sample=SAMPLE,
        isoquant_dir=str(iq_dir),
        gene_h5ad=str(splice / f"{SAMPLE}.gene.h5ad"),
        transcript_h5ad=str(splice / f"{SAMPLE}.transcript.spliced.h5ad"),
        barcode_qc=str(splice / f"{SAMPLE}.barcode_qc.tsv"),
    )
    data = load_report_data(config)
    assert data.gene_adata is not None
    assert data.assignments is not None
    metrics = headline_metrics(data)
    assert metrics["n_cells"] == 3
    figs = build_all_figures(data)
    assert len(figs) > 0
    out = tmp_path / "report.html"
    export_html(data, out)
    assert out.exists()
    assert out.stat().st_size > 100


def test_gene_isoform_detection_mismatched_cells():
    """Spliced transcript matrix may have fewer cells than gene matrix."""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    from isoquant_report.metrics.isoform_stats import gene_isoform_detection

    gene = ad.AnnData(
        X=sp.csr_matrix([[1, 0], [0, 1], [1, 1]]),
        obs=pd.DataFrame(index=["A", "B", "C"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    tx = ad.AnnData(
        X=sp.csr_matrix([[1, 0], [0, 1]]),
        obs=pd.DataFrame(index=["A", "B"]),
        var=pd.DataFrame(index=["t1", "t2"]),
    )
    det = gene_isoform_detection(gene, tx)
    assert len(det) == 2
    assert set(det["barcode"]) == {"A", "B"}
