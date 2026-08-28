#!/usr/bin/env python3
"""
longshot_sc.py -- Per-cell longshot SNV workflow (wf-single-cell port).

Subcommands mirror epi2me-labs/wf-single-cell ``subworkflows/snv.nf``:
  call1            Round-1 per-cell dedup, SplitNCigarReads, longshot
  bulk             Bulk longshot on merged exon-split BAMs (per contig)
  merge_candidates Merge bulk + round-1 VCFs into candidate sites
  genotype2        Round-2 per-cell genotyping at candidate sites
  merge_matrix     Merge cell VCFs + build sparse genotype MEX matrix
  clip_depth       Clip read depth at high-coverage regions (genotype2 helper)

Usage
-----
    python longshot_sc.py call1 --sample-dir snv/sample --ref snv/ref.fa ...
"""

from __future__ import annotations

import argparse
import glob
import gzip
import os
import random
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _run(cmd: list[str], cwd: str | None = None) -> None:
    """Run a shell command, raising on failure."""
    print("+", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, cwd=cwd)


def _glob_files(pattern: str) -> list[str]:
    return sorted(glob.glob(pattern))


def _barcode_from_bam(path: str) -> str:
    return Path(path).stem


# ---------------------------------------------------------------------------
# clip_depth (wf workflow-glue clip_depth)
# ---------------------------------------------------------------------------

def clip_depth(bed_path: str, bam_in: str, bam_out: str, target_depth: int = 200) -> None:
    """
    Reduce per-read depth in high-coverage windows by randomly subsampling reads.

    Approximates wf-single-cell ``workflow-glue clip_depth`` when that tool is
    unavailable: for each (read_id, window) in the BED, keep the read with
    probability target_depth / observed_depth (capped at 1).
    """
    import pysam

    windows: dict[str, list[tuple[str, int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            windows.setdefault(chrom, []).append((chrom, start, end))

    keep_reads: set[str] = set()
    with pysam.AlignmentFile(bam_in, "rb") as bam:
        for chrom, regs in windows.items():
            if chrom not in bam.references:
                continue
            for chrom, start, end in regs:
                reads_in_window = []
                for read in bam.fetch(chrom, start, end):
                    if read.is_unmapped:
                        continue
                    reads_in_window.append(read.query_name)
                if not reads_in_window:
                    continue
                depth = len(reads_in_window)
                if depth <= target_depth:
                    keep_reads.update(reads_in_window)
                else:
                    random.shuffle(reads_in_window)
                    keep_reads.update(reads_in_window[:target_depth])

    with pysam.AlignmentFile(bam_in, "rb") as inp, pysam.AlignmentFile(
        bam_out, "wb", template=inp
    ) as out:
        for read in inp.fetch(until_eof=True):
            if read.query_name in keep_reads or read.is_unmapped:
                out.write(read)


# ---------------------------------------------------------------------------
# call1
# ---------------------------------------------------------------------------

def _call1_one(
    cell_bam: str,
    ref_fa: str,
    out_dir: str,
    min_alt_count: int,
    min_cov: int,
) -> tuple[str, str, str]:
    """Process one cell BAM; return (barcode, vcf_gz, exon_split_bam)."""
    bc = _barcode_from_bam(cell_bam)
    work = os.path.join(out_dir, bc)
    os.makedirs(work, exist_ok=True)
    tmp = os.path.join(work, "tmp_dir")
    os.makedirs(tmp, exist_ok=True)

    dedup = os.path.join(work, "dedup.bam")
    exon_tmp = os.path.join(work, "exon_split_tmp.bam")
    exon_split = os.path.join(out_dir, "exon_split_bams", f"{bc}.bam")
    os.makedirs(os.path.dirname(exon_split), exist_ok=True)
    vcf = os.path.join(out_dir, "round1_vcfs", f"{bc}.vcf")
    os.makedirs(os.path.dirname(vcf), exist_ok=True)

    _run(
        [
            "umi_tools", "dedup",
            "--umi-tag", "UB",
            "--temp-dir", tmp,
            "--cell-tag", "CB",
            "--extract-umi-method=tag",
            "--method", "unique",
            "--per-cell",
            "-I", cell_bam,
            "-S", dedup,
        ],
        cwd=work,
    )
    _run(
        [
            "gatk", "SplitNCigarReads",
            "-R", ref_fa,
            "--input", dedup,
            "--output", exon_tmp,
            "--tmp-dir", tmp,
            "--java-options", "-Xmx8g",
            "--max-reads-in-memory", "100000",
        ],
        cwd=work,
    )
    with open(exon_split, "wb") as outfh:
        subprocess.run(
            ["samtools", "reheader", "-c", "grep -v ^@PG", exon_tmp],
            stdout=outfh,
            check=True,
        )
    _run(["samtools", "index", exon_split])
    if os.path.exists(exon_tmp):
        os.remove(exon_tmp)

    _run(
        [
            "longshot",
            "--bam", exon_split,
            "--ref", ref_fa,
            "--min_alt_count", str(min_alt_count),
            "--min_cov", str(min_cov),
            "--sample_id", bc,
            "--out", vcf,
        ],
        cwd=work,
    )
    vcf_gz = vcf + ".gz"
    _run(["bgzip", "-f", vcf])
    _run(["tabix", vcf_gz])
    return bc, vcf_gz, exon_split


def run_call1(
    sample_dir: str,
    ref_fa: str,
    threads: int,
    min_alt_count: int = 1,
    min_cov: int = 1,
) -> None:
    per_cell = os.path.join(sample_dir, "per_cell_bams")
    out_dir = os.path.join(sample_dir, "call1")
    os.makedirs(out_dir, exist_ok=True)
    cell_bams = _glob_files(os.path.join(per_cell, "*.bam"))
    if not cell_bams:
        raise SystemExit(f"No cell BAMs in {per_cell}")

    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(_call1_one, bam, ref_fa, out_dir, min_alt_count, min_cov)
            for bam in cell_bams
        ]
        for fut in as_completed(futures):
            fut.result()

    sentinel = os.path.join(sample_dir, "call1.sentinel")
    Path(sentinel).write_text("ok\n")


# ---------------------------------------------------------------------------
# bulk
# ---------------------------------------------------------------------------

def run_bulk(
    sample_dir: str,
    ref_fa: str,
    min_alt_count: int = 2,
    min_cov: int = 2,
) -> None:
    call1_dir = os.path.join(sample_dir, "call1", "exon_split_bams")
    exon_bams = _glob_files(os.path.join(call1_dir, "*.bam"))
    if not exon_bams:
        raise SystemExit(f"No exon-split BAMs in {call1_dir}")

    merged = os.path.join(sample_dir, "bulk", "merged.bam")
    os.makedirs(os.path.dirname(merged), exist_ok=True)
    _run(
        [
            "samtools", "merge", "--no-PG",
            "-o", merged,
            *exon_bams,
        ]
    )
    _run(["samtools", "index", merged])

    contigs = []
    proc = subprocess.run(
        ["samtools", "idxstats", merged],
        capture_output=True,
        text=True,
        check=True,
    )
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and int(parts[2]) > 0:
            contigs.append(parts[0])

    bulk_dir = os.path.join(sample_dir, "bulk", "chunks")
    os.makedirs(bulk_dir, exist_ok=True)
    vcfs = []
    for i, chrom in enumerate(contigs):
        vcf = os.path.join(bulk_dir, f"snv_bulk_{i}.vcf")
        _run(
            [
                "longshot",
                "--bam", merged,
                "--ref", ref_fa,
                "--min_alt_count", str(min_alt_count),
                "--min_cov", str(min_cov),
                "--out", vcf,
                "--region", chrom,
            ]
        )
        vcf_gz = vcf + ".gz"
        _run(["bgzip", "-f", vcf])
        _run(["tabix", vcf_gz])
        vcfs.append(vcf_gz)

    list_path = os.path.join(sample_dir, "bulk", "concat_list.txt")
    with open(list_path, "w") as fh:
        for v in vcfs:
            fh.write(v + "\n")
    out_vcf = os.path.join(sample_dir, "bulk_merged.vcf.gz")
    _run(
        [
            "bcftools", "concat",
            "--file-list", list_path,
            "--output-type", "z",
            "-o", out_vcf,
        ]
    )
    _run(["tabix", out_vcf])


# ---------------------------------------------------------------------------
# merge_candidates
# ---------------------------------------------------------------------------

def run_merge_candidates(sample_dir: str, merge_threads: int) -> None:
    bulk_vcf = os.path.join(sample_dir, "bulk_merged.vcf.gz")
    round1 = _glob_files(os.path.join(sample_dir, "call1", "round1_vcfs", "*.vcf.gz"))
    if not os.path.exists(bulk_vcf):
        raise SystemExit(f"Missing bulk VCF: {bulk_vcf}")

    merge_list = os.path.join(sample_dir, "merge_list.txt")
    with open(merge_list, "w") as fh:
        fh.write(bulk_vcf + "\n")
        for v in round1:
            _run(["tabix", v])
            fh.write(v + "\n")

    candidates = os.path.join(sample_dir, "candidates.vcf.gz")
    norm_cmd = (
        f"bcftools merge --file-list {merge_list} "
        f"| bcftools norm --atomize --multiallelics - "
        f"| bcftools +fill-tags - -- -t all "
        f"| bcftools view -s SAMPLE --output-type z -o {candidates}"
    )
    subprocess.run(norm_cmd, shell=True, check=True)
    _run(["tabix", candidates])


# ---------------------------------------------------------------------------
# genotype2
# ---------------------------------------------------------------------------

def _genotype2_one(
    cell_bam: str,
    candidates_vcf: str,
    ref_fa: str,
    ref_fai: str,
    out_dir: str,
    depth_target: int,
    min_alt_count: int,
    min_cov: int,
) -> str:
    bc = _barcode_from_bam(cell_bam)
    work = os.path.join(out_dir, bc)
    os.makedirs(work, exist_ok=True)
    cell_candidates = os.path.join(work, "per_cell_candidates.vcf.gz")
    raw_vcf = os.path.join(work, "per_cell_candidates.vcf")

    with open(raw_vcf, "wb") as outfh:
        subprocess.run(
            [
                "bedtools", "intersect",
                "-header", "-u",
                "-a", candidates_vcf,
                "-b", cell_bam,
            ],
            stdout=outfh,
            check=True,
        )
    _run(["bgzip", "-f", raw_vcf])
    _run(["tabix", cell_candidates])

    n_candidates = int(
        subprocess.run(
            ["bash", "-c", f"zcat {cell_candidates} | grep -vc '^#' || true"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        or "0"
    )

    longshot_input = cell_bam
    if n_candidates > 0:
        chr_sizes = os.path.join(work, "chr.sizes")
        with open(ref_fai) as inf, open(chr_sizes, "w") as outf:
            for line in inf:
                parts = line.split("\t")
                outf.write(f"{parts[0]}\t{parts[1]}\n")

        variant_windows = os.path.join(work, "variant_windows.bed")
        subprocess.run(
            [
                "bash", "-c",
                (
                    f"zgrep -v '^#' {cell_candidates} "
                    r"| awk 'BEGIN {OFS=\"\t\"} {print $1, $2, $2}' | uniq "
                    f"| bedtools slop -i - -g {chr_sizes} -b 100 "
                    f"| bedtools merge > {variant_windows}"
                ),
            ],
            check=True,
        )

        if os.path.getsize(variant_windows) > 0:
            md_prefix = os.path.join(work, "mosdepth", "md")
            os.makedirs(os.path.dirname(md_prefix), exist_ok=True)
            _run(
                [
                    "mosdepth",
                    "--no-per-base",
                    "--by", variant_windows,
                    "--threads", "4",
                    "--fast-mode",
                    md_prefix,
                    cell_bam,
                ]
            )
            regions = md_prefix + ".regions.bed.gz"
            high_cov = os.path.join(work, "high_cov.bed")
            with open(high_cov, "w") as outfh:
                subprocess.run(
                    [
                        "csvtk", "filter",
                        "--tabs", "--no-header-row",
                        "--filter", "4>220",
                        regions,
                    ],
                    stdout=outfh,
                    check=True,
                )
            if os.path.getsize(high_cov) > 0:
                clipped = os.path.join(work, "depth_clipped.bam")
                clip_depth(high_cov, cell_bam, clipped, depth_target)
                _run(["samtools", "index", clipped])
                longshot_input = clipped

    vcf_out = os.path.join(work, f"{bc}.vcf")
    debug_dir = os.path.join(work, "variant_debug_dir")
    _run(
        [
            "longshot",
            "--potential_variants", cell_candidates,
            "--bam", longshot_input,
            "--ref", ref_fa,
            "--min_alt_count", str(min_alt_count),
            "--min_cov", str(min_cov),
            "--sample_id", bc,
            "--out", vcf_out,
            "--force_overwrite",
            "--variant_debug_dir", debug_dir,
        ]
    )
    final_vcf = os.path.join(debug_dir, "4.0.final_genotypes.vcf")
    if not os.path.exists(final_vcf):
        alt = os.path.join(debug_dir, "1.0.potential_SNVs.vcf")
        if os.path.exists(alt):
            shutil.copy(alt, final_vcf)
        else:
            Path(final_vcf).write_text("##fileformat=VCFv4.2\n")

    out_vcf_gz = os.path.join(out_dir, "round2_vcfs", f"{bc}.vcf.gz")
    os.makedirs(os.path.dirname(out_vcf_gz), exist_ok=True)
    _run(["bgzip", "-f", final_vcf])
    shutil.move(final_vcf + ".gz", out_vcf_gz)
    _run(["tabix", out_vcf_gz])
    return out_vcf_gz


def run_genotype2(
    sample_dir: str,
    ref_fa: str,
    ref_fai: str,
    threads: int,
    depth_target: int = 200,
    min_alt_count: int = 1,
    min_cov: int = 1,
) -> None:
    candidates = os.path.join(sample_dir, "candidates.vcf.gz")
    per_cell = os.path.join(sample_dir, "per_cell_bams")
    out_dir = os.path.join(sample_dir, "genotype2")
    os.makedirs(out_dir, exist_ok=True)
    cell_bams = _glob_files(os.path.join(per_cell, "*.bam"))
    if not os.path.exists(candidates):
        raise SystemExit(f"Missing candidates VCF: {candidates}")

    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(
                _genotype2_one,
                bam,
                candidates,
                ref_fa,
                ref_fai,
                out_dir,
                depth_target,
                min_alt_count,
                min_cov,
            )
            for bam in cell_bams
        ]
        for fut in as_completed(futures):
            fut.result()

    Path(os.path.join(sample_dir, "genotype2.sentinel")).write_text("ok\n")


# ---------------------------------------------------------------------------
# merge_matrix / variant_mex
# ---------------------------------------------------------------------------

def _gt_to_code(gt: str) -> int | None:
    """Encode genotype as 0=hom ref, 1=het, 2=hom alt; None if missing."""
    if not gt or gt in (".", "./.", ".|."):
        return None
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return None
    alt_count = sum(1 for a in alleles if a != "0")
    if alt_count == 0:
        return 0
    if alt_count == len(alleles):
        return 2
    return 1


def variant_mex(vcf_gz: str, out_dir: str, report_variants: str | None = None) -> None:
    """Build sparse MEX genotype matrix from merged per-cell VCF."""
    import gzip as gzmod

    os.makedirs(out_dir, exist_ok=True)
    feature_filter: set[str] | None = None
    if report_variants and os.path.exists(report_variants):
        feature_filter = set()
        proc = subprocess.run(
            ["bcftools", "query", "--format", "%CHROM\\_%POS", report_variants],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in proc.stdout.splitlines()[:50]:
            feature_filter.add(line.strip())

    # Parse VCF: one row per (variant, sample) genotype
    fmt = "[%CHROM\\_%POS\\t%SAMPLE\\t%GT\\n]"
    proc = subprocess.run(
        ["bcftools", "query", "-f", fmt, vcf_gz],
        capture_output=True,
        text=True,
        check=True,
    )

    barcodes: list[str] = []
    barcode_idx: dict[str, int] = {}
    feature_idx: dict[str, int] = {}
    entries: list[tuple[int, int, int]] = []

    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        feat = parts[0]
        bc = parts[1]
        gt = parts[2] if len(parts) > 2 else "."
        if feature_filter is not None and feat not in feature_filter:
            continue
        code = _gt_to_code(gt)
        if code is None:
            continue
        if bc not in barcode_idx:
            barcode_idx[bc] = len(barcodes)
            barcodes.append(bc)
        if feat not in feature_idx:
            feature_idx[feat] = len(feature_idx)
        entries.append((feature_idx[feat] + 1, barcode_idx[bc] + 1, code))

    n_rows = len(feature_idx)
    n_cols = len(barcodes)
    n_nnz = len(entries)

    mtx_path = os.path.join(out_dir, "matrix.mtx.gz")
    with gzmod.open(mtx_path, "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write("%metadata_json: {\"format_version\": 1}\n")
        fh.write(f"{n_rows} {n_cols} {n_nnz}\n")
        for row, col, val in entries:
            fh.write(f"{row} {col} {val}\n")

    with gzmod.open(os.path.join(out_dir, "barcodes.tsv.gz"), "wt") as fh:
        for bc in barcodes:
            fh.write(bc + "\n")

    with gzmod.open(os.path.join(out_dir, "features.tsv.gz"), "wt") as fh:
        for feat in sorted(feature_idx, key=lambda x: feature_idx[x]):
            fh.write(feat + "\n")


def run_merge_matrix(
    sample_dir: str,
    sample: str,
    merge_threads: int,
    report_variants: str | None = None,
) -> None:
    round2 = _glob_files(
        os.path.join(sample_dir, "genotype2", "round2_vcfs", "*.vcf.gz")
    )
    if not round2:
        raise SystemExit("No round-2 VCFs found")

    merge_list = os.path.join(sample_dir, "final_merge_list.txt")
    with open(merge_list, "w") as fh:
        for v in round2:
            _run(["tabix", v])
            fh.write(v + "\n")

    final_vcf = os.path.join(sample_dir, f"{sample}.final_merged.vcf.gz")
    cmd = (
        f"bcftools merge --file-list {merge_list} "
        f"| bcftools norm --atomize --multiallelics - "
        f"| bcftools +fill-tags - -- -t all "
        f"| bgzip -@ {merge_threads} > {final_vcf}"
    )
    subprocess.run(cmd, shell=True, check=True)
    _run(["tabix", final_vcf])

    gmatrix = os.path.join(sample_dir, f"{sample}.genotype_matrix")
    variant_mex(final_vcf, gmatrix, report_variants)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_call1 = sub.add_parser("call1")
    p_call1.add_argument("--sample-dir", required=True)
    p_call1.add_argument("--ref", required=True)
    p_call1.add_argument("--threads", type=int, default=4)
    p_call1.add_argument("--min-alt-count", type=int, default=1)
    p_call1.add_argument("--min-cov", type=int, default=1)

    p_bulk = sub.add_parser("bulk")
    p_bulk.add_argument("--sample-dir", required=True)
    p_bulk.add_argument("--ref", required=True)
    p_bulk.add_argument("--min-alt-count", type=int, default=2)
    p_bulk.add_argument("--min-cov", type=int, default=2)

    p_merge = sub.add_parser("merge_candidates")
    p_merge.add_argument("--sample-dir", required=True)
    p_merge.add_argument("--merge-threads", type=int, default=8)

    p_g2 = sub.add_parser("genotype2")
    p_g2.add_argument("--sample-dir", required=True)
    p_g2.add_argument("--ref", required=True)
    p_g2.add_argument("--ref-fai", required=True)
    p_g2.add_argument("--threads", type=int, default=4)
    p_g2.add_argument("--depth-target", type=int, default=200)
    p_g2.add_argument("--min-alt-count", type=int, default=1)
    p_g2.add_argument("--min-cov", type=int, default=1)

    p_mtx = sub.add_parser("merge_matrix")
    p_mtx.add_argument("--sample-dir", required=True)
    p_mtx.add_argument("--sample", required=True)
    p_mtx.add_argument("--merge-threads", type=int, default=8)
    p_mtx.add_argument("--report-variants", default="")

    p_clip = sub.add_parser("clip_depth")
    p_clip.add_argument("--bed", required=True)
    p_clip.add_argument("--bam-in", required=True)
    p_clip.add_argument("--bam-out", required=True)
    p_clip.add_argument("--target-depth", type=int, default=200)

    args = parser.parse_args(argv)

    if args.command == "call1":
        run_call1(
            args.sample_dir,
            args.ref,
            args.threads,
            args.min_alt_count,
            args.min_cov,
        )
    elif args.command == "bulk":
        run_bulk(
            args.sample_dir,
            args.ref,
            args.min_alt_count,
            args.min_cov,
        )
    elif args.command == "merge_candidates":
        run_merge_candidates(args.sample_dir, args.merge_threads)
    elif args.command == "genotype2":
        run_genotype2(
            args.sample_dir,
            args.ref,
            args.ref_fai,
            args.threads,
            args.depth_target,
            args.min_alt_count,
            args.min_cov,
        )
    elif args.command == "merge_matrix":
        rv = args.report_variants or None
        run_merge_matrix(
            args.sample_dir,
            args.sample,
            args.merge_threads,
            rv,
        )
    elif args.command == "clip_depth":
        clip_depth(args.bed, args.bam_in, args.bam_out, args.target_depth)
    return 0


if __name__ == "__main__":
    sys.exit(main())
