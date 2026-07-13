"""Generate minimal synthetic IsoQuant fixtures for smoke tests."""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

SAMPLE = "test_sample"
FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "synthetic_isoquant"


def _write_fixtures() -> Path:
    iq = FIXTURE_ROOT / "isoquant" / SAMPLE / SAMPLE
    splice = FIXTURE_ROOT / "splice_matrices"
    iq.mkdir(parents=True, exist_ok=True)
    splice.mkdir(parents=True, exist_ok=True)

    barcodes = ["AAACCC", "TTTGGG", "GGGAAA"]
    genes = ["ENSG00001", "ENSG00002", "ENSG00003"]
    txs = ["ENST00001", "ENST00002", "ENST00003"]

    # Grouped gene counts linear
    rows = []
    for g in genes:
        for b in barcodes:
            rows.append([g, b, np.random.randint(1, 20)])
    pd.DataFrame(rows, columns=["feature_id", "group_id", "count"]).to_csv(
        iq / f"{SAMPLE}.gene_grouped_barcode_counts.linear.tsv",
        sep="\t",
        header=False,
        index=False,
    )

    # Transcript grouped linear
    for t in txs:
        for b in barcodes:
            rows.append([t, b, np.random.randint(1, 10)])
    pd.DataFrame(rows[-9:], columns=["feature_id", "group_id", "count"]).to_csv(
        iq / f"{SAMPLE}.transcript_grouped_barcode_counts.linear.tsv",
        sep="\t",
        header=False,
        index=False,
    )

    # read_assignments
    assign_rows = [
        ["read1", "chr1", "+", "ENST00001", "ENSG00001", "unique", "", "1", "", "AAACCC"],
        ["read2", "chr1", "+", "ENST00002", "ENSG00002", "inconsistent", "intron_retention", "2", "", "TTTGGG"],
        ["read3", "chr1", "-", ".", ".", "intergenic", "", "", "", "GGGAAA"],
    ]
    with gzip.open(iq / f"{SAMPLE}.read_assignments.tsv.gz", "wt") as fh:
        fh.write("#read_id\tchr\tstrand\tisoform_id\tgene_id\tassignment_type\tassignment_events\texons\tadditional_info\tgroups\n")
        for r in assign_rows:
            fh.write("\t".join(r) + "\n")

    # allinfo (no header)
    with gzip.open(iq / f"{SAMPLE}.UMI_filtered.ED3.allinfo.gz", "wt") as fh:
        fh.write("read1\tENSG00001\tNone\tAAACCC\tAGACGATGTAAA\t.\tNoTSS\tNoPolyA\t.\tknown\t1\tENST00001\tknown\n")
        fh.write("read2\tENSG00002\tNone\tTTTGGG\tCGACGATGTAAA\t.\tNoTSS\tNoPolyA\t.\tnovel\t2\tENST00002\tunknown_type\n")

    # minimal GTF
    (iq / f"{SAMPLE}.transcript_models.gtf").write_text(
        'chr1\tIsoQuant\texon\t100\t200\t.\t+\t.\tgene_id "ENSG00001"; transcript_id "ENST00001";\n'
        'chr1\tIsoQuant\texon\t300\t400\t.\t+\t.\tgene_id "ENSG00002"; transcript_id "ENST00002";\n'
    )

    # h5ad
    Xg = sp.csr_matrix(np.random.poisson(3, (len(barcodes), len(genes))).astype(float))
    gene = ad.AnnData(X=Xg, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
    gene.layers["total"] = Xg.copy()
    gene.layers["spliced"] = Xg.copy()
    gene.layers["unspliced"] = sp.csr_matrix(Xg.shape)
    gene.write_h5ad(splice / f"{SAMPLE}.gene.h5ad")

    Xt = sp.csr_matrix(np.random.poisson(2, (len(barcodes), len(txs))).astype(float))
    tx = ad.AnnData(X=Xt, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=txs))
    tx.write_h5ad(splice / f"{SAMPLE}.transcript.spliced.h5ad")

    pd.DataFrame(
        {
            "barcode": barcodes,
            "total_umis": [10, 8, 5],
            "spliced_umis": [7, 4, 2],
            "unspliced_umis": [2, 3, 2],
            "ambiguous_umis": [1, 1, 1],
            "unspliced_fraction": [0.2, 0.375, 0.4],
            "n_genes": [2, 2, 1],
        }
    ).to_csv(splice / f"{SAMPLE}.barcode_qc.tsv", sep="\t", index=False)

    return FIXTURE_ROOT


if __name__ == "__main__":
    print(_write_fixtures())
