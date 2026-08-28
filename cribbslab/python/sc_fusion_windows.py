#!/usr/bin/env python3
"""
sc_fusion_windows.py -- Single-cell targeted translocation scanner.

Adaptation of split-mm whitelist-based fusion window scanning for barcode-tagged
long-read BAMs from pipeline_sclong. Counts reads (and UMIs) whose alignments
land in both genomic windows of a chromosome-pair template, aggregated per
cell barcode.

Usage
-----
    python sc_fusion_windows.py \\
        --bam tagged_cells/sample.cells.bam \\
        --whitelist fusion_whitelist.tsv \\
        --sample sample \\
        --outdir targeted_fusions
"""

from __future__ import annotations

import argparse
import csv
import logging
import statistics
import subprocess
import shutil
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam


DEFAULT_WINDOW_SIZE = 1_000_000
NO_BARCODE = "__no_barcode__"


@dataclass
class ReadHit:
    cb: str = ""
    ub: str = ""
    is_split: bool = False


@dataclass
class Target:
    cluster: str
    chr_left: str
    chr_right: str
    window_size: int
    pad_left: float
    pad_right: float
    padding_units: str
    target_mode: str
    line_no: int
    bin_left: Optional[int] = None
    bin_right: Optional[int] = None
    mid_left_bp: Optional[int] = None
    mid_right_bp: Optional[int] = None
    exact_left_start: int = 0
    exact_left_end: int = 0
    exact_right_start: int = 0
    exact_right_end: int = 0
    broad_left_start: int = 0
    broad_left_end: int = 0
    broad_right_start: int = 0
    broad_right_end: int = 0


@dataclass
class CellSupport:
    reads: Set[str] = field(default_factory=set)
    split_reads: Set[str] = field(default_factory=set)
    umis: Set[str] = field(default_factory=set)


def bin_to_interval(bin_id: int, window_size: int) -> Tuple[int, int]:
    if bin_id < 1:
        raise ValueError(f"bin_id must be >= 1, got {bin_id}")
    start = (bin_id - 1) * window_size + 1
    end = bin_id * window_size
    return start, end


def midpoint_to_interval(mid_bp: int, window_size: int) -> Tuple[int, int]:
    if mid_bp < 1:
        raise ValueError(f"midpoint must be >= 1, got {mid_bp}")
    half = window_size // 2
    start = max(1, mid_bp - half)
    end = start + window_size - 1
    return start, end


def expand_interval(
    start: int,
    end: int,
    pad_value: float,
    window_size: int,
    padding_units: str,
) -> Tuple[int, int]:
    if pad_value < 0:
        raise ValueError(f"pad_value must be >= 0, got {pad_value}")
    if padding_units == "windows":
        pad_bp = pad_value * window_size
    elif padding_units == "bp":
        pad_bp = pad_value
    else:
        raise ValueError(f"Unknown padding_units: {padding_units}")
    pad_bp_int = int(round(pad_bp))
    return max(1, start - pad_bp_int), end + pad_bp_int


def build_contig_resolver(bam: pysam.AlignmentFile):
    refs = set(bam.references)
    alias: Dict[str, str] = {}
    for r in refs:
        alias[r] = r
        if r.startswith("chr"):
            alias[r[3:]] = r
        else:
            alias["chr" + r] = r

    def resolve(name: str) -> Optional[str]:
        name = str(name).strip()
        if name in alias:
            return alias[name]
        if name.startswith("chr") and name[3:] in alias:
            return alias[name[3:]]
        if not name.startswith("chr") and ("chr" + name) in alias:
            return alias["chr" + name]
        return None

    return resolve


def parse_optional_int(value: Optional[str], *, default: Optional[int] = None) -> Optional[int]:
    if value is None or value == "":
        return default
    return int(value)


def parse_optional_float(value: Optional[str], *, default: float = 1.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def load_whitelist(path: str, default_window_size: int) -> List[dict]:
    rows: List[dict] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"cluster", "chr_left", "chr_right"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Whitelist missing columns: " + ", ".join(sorted(missing)))
        for line_no, row in enumerate(reader, start=2):
            rows.append(
                {
                    "cluster": row["cluster"].strip(),
                    "chr_left": row["chr_left"].strip(),
                    "chr_right": row["chr_right"].strip(),
                    "bin_left": parse_optional_int(row.get("bin_left")),
                    "bin_right": parse_optional_int(row.get("bin_right")),
                    "mid_left_bp": parse_optional_int(row.get("mid_left_bp")),
                    "mid_right_bp": parse_optional_int(row.get("mid_right_bp")),
                    "window_size": parse_optional_int(
                        row.get("window_size"), default=default_window_size
                    ),
                    "pad_left": parse_optional_float(row.get("pad_left"), default=1.0),
                    "pad_right": parse_optional_float(row.get("pad_right"), default=1.0),
                    "line_no": line_no,
                }
            )
    return rows


def prepare_targets(
    raw_rows: List[dict],
    bam: pysam.AlignmentFile,
    default_padding_units: str,
    logger: logging.Logger,
) -> List[Target]:
    resolve = build_contig_resolver(bam)
    targets: List[Target] = []

    for row in raw_rows:
        chr_left = resolve(row["chr_left"])
        chr_right = resolve(row["chr_right"])
        if chr_left is None:
            raise ValueError(
                f"Could not resolve chr_left='{row['chr_left']}' "
                f"(whitelist line {row['line_no']})"
            )
        if chr_right is None:
            raise ValueError(
                f"Could not resolve chr_right='{row['chr_right']}' "
                f"(whitelist line {row['line_no']})"
            )

        window_size = row["window_size"]
        if window_size is None or window_size < 1:
            raise ValueError(
                f"Invalid window_size at whitelist line {row['line_no']}: {window_size}"
            )

        has_midpoints = row["mid_left_bp"] is not None or row["mid_right_bp"] is not None
        has_bins = row["bin_left"] is not None or row["bin_right"] is not None

        if has_midpoints:
            if row["mid_left_bp"] is None or row["mid_right_bp"] is None:
                raise ValueError(
                    f"Both mid_left_bp and mid_right_bp required "
                    f"(whitelist line {row['line_no']})"
                )
            exact_left_start, exact_left_end = midpoint_to_interval(
                row["mid_left_bp"], window_size
            )
            exact_right_start, exact_right_end = midpoint_to_interval(
                row["mid_right_bp"], window_size
            )
            target_mode = "midpoint"
            bin_left = bin_right = None
            mid_left_bp = row["mid_left_bp"]
            mid_right_bp = row["mid_right_bp"]
        elif has_bins:
            if row["bin_left"] is None or row["bin_right"] is None:
                raise ValueError(
                    f"Both bin_left and bin_right required "
                    f"(whitelist line {row['line_no']})"
                )
            exact_left_start, exact_left_end = bin_to_interval(
                row["bin_left"], window_size
            )
            exact_right_start, exact_right_end = bin_to_interval(
                row["bin_right"], window_size
            )
            target_mode = "bin"
            bin_left = row["bin_left"]
            bin_right = row["bin_right"]
            mid_left_bp = mid_right_bp = None
        else:
            raise ValueError(
                f"Whitelist line {row['line_no']} must provide "
                f"bin_left/bin_right or mid_left_bp/mid_right_bp"
            )

        broad_left_start, broad_left_end = expand_interval(
            exact_left_start, exact_left_end,
            row["pad_left"], window_size, default_padding_units,
        )
        broad_right_start, broad_right_end = expand_interval(
            exact_right_start, exact_right_end,
            row["pad_right"], window_size, default_padding_units,
        )

        targets.append(
            Target(
                cluster=row["cluster"],
                chr_left=chr_left,
                chr_right=chr_right,
                window_size=window_size,
                pad_left=row["pad_left"],
                pad_right=row["pad_right"],
                padding_units=default_padding_units,
                target_mode=target_mode,
                line_no=row["line_no"],
                bin_left=bin_left,
                bin_right=bin_right,
                mid_left_bp=mid_left_bp,
                mid_right_bp=mid_right_bp,
                exact_left_start=exact_left_start,
                exact_left_end=exact_left_end,
                exact_right_start=exact_right_start,
                exact_right_end=exact_right_end,
                broad_left_start=broad_left_start,
                broad_left_end=broad_left_end,
                broad_right_start=broad_right_start,
                broad_right_end=broad_right_end,
            )
        )

    logger.info("Loaded %d whitelist targets", len(targets))
    return targets


def configure_logging(log_path: Path, verbose: bool) -> logging.Logger:
    logger = logging.getLogger("sc_fusion_windows")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    for h in list(logger.handlers):
        logger.removeHandler(h)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, mode="w")
    fh.setLevel(logging.DEBUG if verbose else logging.INFO)
    fh.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.ERROR)
    ch.setFormatter(logging.Formatter("%(levelname)s | %(message)s"))
    logger.addHandler(ch)
    return logger


def is_split_alignment(aln) -> bool:
    return aln.is_supplementary or aln.has_tag("SA")


def _read_tags(aln) -> Tuple[str, str]:
    cb = aln.get_tag("CB") if aln.has_tag("CB") else ""
    ub = aln.get_tag("UB") if aln.has_tag("UB") else ""
    return str(cb), str(ub)


def collect_reads_in_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
    min_mapq: int,
    keep_secondary: bool,
) -> Tuple[Dict[str, ReadHit], List[int]]:
    """Return read_name -> ReadHit for alignments overlapping the interval."""
    hits: Dict[str, ReadHit] = {}
    mapqs: List[int] = []
    fetch_start0 = max(0, start_1based - 1)
    fetch_end0 = end_1based

    for aln in bam.fetch(chrom, fetch_start0, fetch_end0):
        if aln.is_unmapped:
            continue
        if not keep_secondary and aln.is_secondary:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        rid = aln.query_name
        cb, ub = _read_tags(aln)
        split = is_split_alignment(aln)
        mapqs.append(aln.mapping_quality)

        if rid not in hits:
            hits[rid] = ReadHit(cb=cb, ub=ub, is_split=split)
        else:
            h = hits[rid]
            if not h.cb and cb:
                h.cb = cb
            if not h.ub and ub:
                h.ub = ub
            if split:
                h.is_split = True

    return hits, mapqs


def _barcode_key(cb: str) -> str:
    return cb if cb else NO_BARCODE


def _umi_key(cb: str, ub: str) -> str:
    if not cb or not ub:
        return ""
    return f"{cb}:{ub}"


def aggregate_by_cell(
    supporting_reads: Dict[str, ReadHit],
) -> Dict[str, CellSupport]:
    """Aggregate supporting reads by cell barcode."""
    by_cell: Dict[str, CellSupport] = defaultdict(CellSupport)
    for rid, hit in supporting_reads.items():
        bc = _barcode_key(hit.cb)
        cell = by_cell[bc]
        cell.reads.add(rid)
        if hit.is_split:
            cell.split_reads.add(rid)
        umi = _umi_key(hit.cb, hit.ub)
        if umi:
            cell.umis.add(umi)
    return dict(by_cell)


def scan_target(
    bam: pysam.AlignmentFile,
    target: Target,
    min_mapq: int,
    keep_secondary: bool,
) -> dict:
    left_exact, left_exact_mapqs = collect_reads_in_interval(
        bam, target.chr_left,
        target.exact_left_start, target.exact_left_end,
        min_mapq, keep_secondary,
    )
    right_exact, right_exact_mapqs = collect_reads_in_interval(
        bam, target.chr_right,
        target.exact_right_start, target.exact_right_end,
        min_mapq, keep_secondary,
    )
    left_broad, left_broad_mapqs = collect_reads_in_interval(
        bam, target.chr_left,
        target.broad_left_start, target.broad_left_end,
        min_mapq, keep_secondary,
    )
    right_broad, right_broad_mapqs = collect_reads_in_interval(
        bam, target.chr_right,
        target.broad_right_start, target.broad_right_end,
        min_mapq, keep_secondary,
    )

    exact_support_names = set(left_exact) & set(right_exact)
    broad_support_names = set(left_broad) & set(right_broad)

    exact_support = {rid: left_exact[rid] for rid in exact_support_names}
    broad_support = {rid: left_broad[rid] for rid in broad_support_names}

    for rid in exact_support_names:
        h = exact_support[rid]
        if rid in right_exact:
            rh = right_exact[rid]
            if not h.cb and rh.cb:
                h.cb = rh.cb
            if not h.ub and rh.ub:
                h.ub = rh.ub
            if rh.is_split:
                h.is_split = True

    for rid in broad_support_names:
        h = broad_support[rid]
        if rid in right_broad:
            rh = right_broad[rid]
            if not h.cb and rh.cb:
                h.cb = rh.cb
            if not h.ub and rh.ub:
                h.ub = rh.ub
            if rh.is_split:
                h.is_split = True

    exact_split = {rid for rid in exact_support_names if exact_support[rid].is_split}
    broad_split = {rid for rid in broad_support_names if broad_support[rid].is_split}

    exact_by_cell = aggregate_by_cell(exact_support)
    broad_by_cell = aggregate_by_cell(broad_support)

    exact_all_mapqs = left_exact_mapqs + right_exact_mapqs
    broad_all_mapqs = left_broad_mapqs + right_broad_mapqs

    barcoded_cells = {
        bc for bc in exact_by_cell if bc != NO_BARCODE
    } | {bc for bc in broad_by_cell if bc != NO_BARCODE}

    all_umis: Set[str] = set()
    for bc, cs in exact_by_cell.items():
        if bc != NO_BARCODE:
            all_umis.update(cs.umis)

    return {
        "target": target,
        "exact_support": exact_support,
        "broad_support": broad_support,
        "exact_split": exact_split,
        "broad_split": broad_split,
        "exact_by_cell": exact_by_cell,
        "broad_by_cell": broad_by_cell,
        "exact_left_reads": len(left_exact),
        "exact_right_reads": len(right_exact),
        "exact_support_reads": len(exact_support_names),
        "exact_split_reads": len(exact_split),
        "exact_median_mapq": statistics.median(exact_all_mapqs) if exact_all_mapqs else "",
        "broad_left_reads": len(left_broad),
        "broad_right_reads": len(right_broad),
        "broad_support_reads": len(broad_support_names),
        "broad_split_reads": len(broad_split),
        "broad_median_mapq": statistics.median(broad_all_mapqs) if broad_all_mapqs else "",
        "n_cells": len(barcoded_cells),
        "n_umis": len(all_umis),
    }


def run_samtools_flagstat(bam_path: str, logger: logging.Logger) -> Dict[str, int]:
    if shutil.which("samtools") is None:
        raise SystemExit("samtools not found in PATH")
    logger.info("Running samtools flagstat")
    proc = subprocess.run(
        ["samtools", "flagstat", bam_path],
        capture_output=True,
        text=True,
        check=True,
    )
    total = primary_mapped = mapped = secondary = supplementary = 0
    line_re = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+(.+?)(?:\s+\(|$)")
    for raw_line in proc.stdout.splitlines():
        m = line_re.match(raw_line)
        if not m:
            continue
        count = int(m.group(1))
        label = m.group(2).strip()
        if label == "in total":
            total = count
        elif label == "primary mapped":
            primary_mapped = count
        elif label == "mapped":
            mapped = count
        elif label == "secondary":
            secondary = count
        elif label == "supplementary":
            supplementary = count
    if primary_mapped == 0 and mapped:
        primary_mapped = max(mapped - secondary - supplementary, 0)
    return {
        "total_primary_mapped_reads": primary_mapped,
        "total_chimeric_reads": supplementary,
        "samtools_total_reads": total,
        "samtools_mapped_reads": mapped,
        "samtools_secondary_alignments": secondary,
        "samtools_supplementary_alignments": supplementary,
    }


def write_per_target(path: Path, sample: str, results: List[dict], stats: dict) -> None:
    fieldnames = [
        "sample", "cluster", "chr_left", "chr_right", "target_mode", "line_no",
        "bin_left", "bin_right", "mid_left_bp", "mid_right_bp",
        "window_size", "pad_left", "pad_right", "padding_units",
        "total_primary_mapped_reads", "total_chimeric_reads",
        "samtools_total_reads", "samtools_mapped_reads",
        "samtools_secondary_alignments", "samtools_supplementary_alignments",
        "exact_left_interval", "exact_right_interval",
        "broad_left_interval", "broad_right_interval",
        "exact_left_reads", "exact_right_reads",
        "exact_support_reads", "exact_split_reads", "exact_median_mapq",
        "broad_left_reads", "broad_right_reads",
        "broad_support_reads", "broad_split_reads", "broad_median_mapq",
        "n_cells", "n_umis",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            t: Target = r["target"]
            writer.writerow({
                "sample": sample,
                "cluster": t.cluster,
                "chr_left": t.chr_left,
                "chr_right": t.chr_right,
                "target_mode": t.target_mode,
                "line_no": t.line_no,
                "bin_left": t.bin_left if t.bin_left is not None else "",
                "bin_right": t.bin_right if t.bin_right is not None else "",
                "mid_left_bp": t.mid_left_bp if t.mid_left_bp is not None else "",
                "mid_right_bp": t.mid_right_bp if t.mid_right_bp is not None else "",
                "window_size": t.window_size,
                "pad_left": t.pad_left,
                "pad_right": t.pad_right,
                "padding_units": t.padding_units,
                "total_primary_mapped_reads": stats["total_primary_mapped_reads"],
                "total_chimeric_reads": stats["total_chimeric_reads"],
                "samtools_total_reads": stats["samtools_total_reads"],
                "samtools_mapped_reads": stats["samtools_mapped_reads"],
                "samtools_secondary_alignments": stats["samtools_secondary_alignments"],
                "samtools_supplementary_alignments": stats["samtools_supplementary_alignments"],
                "exact_left_interval": f"{t.chr_left}:{t.exact_left_start}-{t.exact_left_end}",
                "exact_right_interval": f"{t.chr_right}:{t.exact_right_start}-{t.exact_right_end}",
                "broad_left_interval": f"{t.chr_left}:{t.broad_left_start}-{t.broad_left_end}",
                "broad_right_interval": f"{t.chr_right}:{t.broad_right_start}-{t.broad_right_end}",
                "exact_left_reads": r["exact_left_reads"],
                "exact_right_reads": r["exact_right_reads"],
                "exact_support_reads": r["exact_support_reads"],
                "exact_split_reads": r["exact_split_reads"],
                "exact_median_mapq": r["exact_median_mapq"],
                "broad_left_reads": r["broad_left_reads"],
                "broad_right_reads": r["broad_right_reads"],
                "broad_support_reads": r["broad_support_reads"],
                "broad_split_reads": r["broad_split_reads"],
                "broad_median_mapq": r["broad_median_mapq"],
                "n_cells": r["n_cells"],
                "n_umis": r["n_umis"],
            })


def write_per_cell(path: Path, sample: str, results: List[dict]) -> None:
    fieldnames = [
        "sample", "cluster", "barcode",
        "exact_reads", "exact_umis", "exact_split_reads",
        "broad_reads", "broad_umis", "broad_split_reads",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            cluster = r["target"].cluster
            all_barcodes = set(r["exact_by_cell"]) | set(r["broad_by_cell"])
            for bc in sorted(all_barcodes):
                exact = r["exact_by_cell"].get(bc, CellSupport())
                broad = r["broad_by_cell"].get(bc, CellSupport())
                if not exact.reads and not broad.reads:
                    continue
                writer.writerow({
                    "sample": sample,
                    "cluster": cluster,
                    "barcode": bc,
                    "exact_reads": len(exact.reads),
                    "exact_umis": len(exact.umis),
                    "exact_split_reads": len(exact.split_reads),
                    "broad_reads": len(broad.reads),
                    "broad_umis": len(broad.umis),
                    "broad_split_reads": len(broad.split_reads),
                })


def write_matrices(
    outdir: Path,
    prefix: str,
    sample: str,
    results: List[dict],
    matrix_window: str,
) -> None:
    clusters = [r["target"].cluster for r in results]
    all_barcodes: Set[str] = set()
    for r in results:
        by_cell = r["exact_by_cell"] if matrix_window == "exact" else r["broad_by_cell"]
        for bc in by_cell:
            if bc != NO_BARCODE:
                all_barcodes.add(bc)
    barcodes = sorted(all_barcodes)

    umi_path = outdir / f"{prefix}.targeted_fusion.matrix.umi.tsv"
    read_path = outdir / f"{prefix}.targeted_fusion.matrix.reads.tsv"

    for path, value_fn in ((umi_path, lambda cs: len(cs.umis)),
                           (read_path, lambda cs: len(cs.reads))):
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["barcode"] + clusters)
            for bc in barcodes:
                row = [bc]
                for r in results:
                    by_cell = (
                        r["exact_by_cell"] if matrix_window == "exact"
                        else r["broad_by_cell"]
                    )
                    cs = by_cell.get(bc, CellSupport())
                    row.append(value_fn(cs))
                writer.writerow(row)


def write_per_read(path: Path, sample: str, results: List[dict]) -> None:
    fieldnames = [
        "sample", "cluster", "read_id", "barcode", "umi",
        "window", "is_split",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            cluster = r["target"].cluster
            for window, support, split_set in (
                ("exact", r["exact_support"], r["exact_split"]),
                ("broad", r["broad_support"], r["broad_split"]),
            ):
                for rid, hit in sorted(support.items()):
                    writer.writerow({
                        "sample": sample,
                        "cluster": cluster,
                        "read_id": rid,
                        "barcode": hit.cb,
                        "umi": hit.ub,
                        "window": window,
                        "is_split": rid in split_set,
                    })


def scan_bam(
    bam_path: str,
    whitelist_path: str,
    sample: str,
    outdir: Path,
    out_prefix: str,
    window_size: int,
    padding_units: str,
    min_mapq: int,
    keep_secondary: bool,
    matrix_window: str,
    verbose: bool,
) -> None:
    log_path = outdir / f"{out_prefix}.targeted_fusion.log"
    logger = configure_logging(log_path, verbose)
    logger.info("Starting single-cell targeted fusion scan")
    logger.info("BAM: %s", bam_path)
    logger.info("Whitelist: %s", whitelist_path)
    logger.info("Sample: %s", sample)

    stats = run_samtools_flagstat(bam_path, logger)
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        if not bam.has_index():
            raise RuntimeError(
                f"BAM index not found for {bam_path}. "
                "Run samtools index before scanning."
            )
        raw_rows = load_whitelist(whitelist_path, window_size)
        targets = prepare_targets(raw_rows, bam, padding_units, logger)
        results = [
            scan_target(bam, t, min_mapq, keep_secondary) for t in targets
        ]
        for r in results:
            t = r["target"]
            logger.debug(
                "%s | exact_support=%d exact_split=%d n_cells=%d n_umis=%d",
                t.cluster,
                r["exact_support_reads"],
                r["exact_split_reads"],
                r["n_cells"],
                r["n_umis"],
            )
    finally:
        bam.close()

    outdir.mkdir(parents=True, exist_ok=True)
    write_per_target(
        outdir / f"{out_prefix}.targeted_fusion.per_target.tsv",
        sample, results, stats,
    )
    write_per_cell(
        outdir / f"{out_prefix}.targeted_fusion.per_cell.tsv",
        sample, results,
    )
    write_matrices(outdir, out_prefix, sample, results, matrix_window)
    write_per_read(
        outdir / f"{out_prefix}.targeted_fusion.per_read.tsv",
        sample, results,
    )
    logger.info("Wrote outputs with prefix %s to %s", out_prefix, outdir)
    logger.info("Done")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--whitelist", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", default="targeted_fusions")
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW_SIZE)
    parser.add_argument("--padding-units", choices=("windows", "bp"), default="windows")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--keep-secondary", action="store_true")
    parser.add_argument(
        "--matrix-window", choices=("exact", "broad"), default="exact",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_prefix = args.out_prefix or args.sample
    scan_bam(
        bam_path=args.bam,
        whitelist_path=args.whitelist,
        sample=args.sample,
        outdir=Path(args.outdir),
        out_prefix=out_prefix,
        window_size=args.window_size,
        padding_units=args.padding_units,
        min_mapq=args.min_mapq,
        keep_secondary=args.keep_secondary,
        matrix_window=args.matrix_window,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
