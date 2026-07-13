# IsoQuant QC Report

Interactive Streamlit QC and exploration report for single-cell long-read
RNA-seq data processed with IsoQuant, plus a headless static HTML export
invoked from `pipeline_sclong`.

## Install

```bash
cd isoquant_report
pip install -e ".[dev]"
```

Or use the `sclong` conda environment (includes Streamlit, plotly, scanpy).

## Interactive app

```bash
streamlit run app.py
```

With a pipeline-written config:

```bash
streamlit run app.py -- --config isoquant_report/SAMPLE/config.yaml
```

## Headless HTML (pipeline / CLI)

```bash
python -m isoquant_report \
  --isoquant-dir isoquant/SAMPLE/SAMPLE \
  --sample SAMPLE \
  --gene-h5ad splice_matrices/SAMPLE.gene.h5ad \
  --transcript-h5ad splice_matrices/SAMPLE.transcript.spliced.h5ad \
  --barcode-qc splice_matrices/SAMPLE.barcode_qc.tsv \
  --outdir isoquant_report/SAMPLE \
  --write-config
```

## Expected inputs

### IsoQuant native (directory `isoquant/{sample}/{sample}/`)

- `{sample}.gene_grouped_barcode_counts.{matrix.mtx,linear.tsv,...}`
- `{sample}.transcript_grouped_barcode_counts.*`
- `{sample}.read_assignments.tsv.gz` or `read_info.tsv.gz`
- `{sample}.UMI_filtered.ED*.allinfo.gz` (no header; 13 columns)
- `{sample}.transcript_models.gtf`

### Derived matrices (from `build_splice_matrices`)

- `splice_matrices/{sample}.gene.h5ad` (layers: total, spliced, unspliced)
- `splice_matrices/{sample}.transcript.spliced.h5ad`
- `splice_matrices/{sample}.barcode_qc.tsv` (optional)

## Tests

```bash
pytest tests/
```
