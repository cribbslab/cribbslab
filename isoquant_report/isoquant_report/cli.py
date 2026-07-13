"""Headless CLI for pipeline static HTML export."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from isoquant_report.config import ReportConfig
from isoquant_report.export_html import export_html
from isoquant_report.report import load_report_data


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate static IsoQuant QC HTML report."
    )
    parser.add_argument("--isoquant-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gene-h5ad", default="")
    parser.add_argument("--transcript-h5ad", default="")
    parser.add_argument("--barcode-qc", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--write-config", action="store_true")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    config = ReportConfig(
        sample=args.sample,
        isoquant_dir=args.isoquant_dir,
        gene_h5ad=args.gene_h5ad,
        transcript_h5ad=args.transcript_h5ad,
        barcode_qc=args.barcode_qc,
    )

    data = load_report_data(config)
    report_path = outdir / "report.html"
    export_html(data, report_path)

    if args.write_config:
        config.write_yaml(outdir / "config.yaml")

    print(f"Wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
