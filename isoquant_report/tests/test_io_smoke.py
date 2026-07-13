"""Smoke tests for IsoQuant report package."""

from pathlib import Path

import pytest

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
