#!/usr/bin/env python3
"""
isoquant_matrices.py -- Build spliced/unspliced AnnData matrices from IsoQuant
single-cell output.

Reads IsoQuant's native barcode-grouped count matrices for molecule totals,
parses the per-read assignment file to label each molecule as spliced,
unspliced, or ambiguous, and writes:

  {outdir}/{sample}.transcript.spliced.h5ad  -- transcript x barcode (spliced)
  {outdir}/{sample}.gene.h5ad               -- gene x barcode with
                                               spliced/unspliced/total layers
  {outdir}/{sample}/{spliced_tx,...}/        -- MatrixMarket fallback exports
  {outdir}/{sample}.barcode_qc.tsv          -- per-barcode QC table

Sparse matrices are used throughout; genome-wide tables are never densified.

Strand note: minimap2 upstream uses -uf (forward-strand only). IsoQuant
inherits this from the BAM. Wrong-strand intronic reads would contaminate the
unspliced matrix; the -uf flag prevents this.
"""

import argparse
import glob
import gzip
import logging
import os
import re
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# Assignment types treated as cleanly spliced (no intron retention)
SPLICED_TYPES = {"unique", "unique_minor_difference"}

# Assignment event tokens that indicate intron retention
IR_TOKENS = {
    "intron_retention",
    "unspliced_intron_retention",
    "incomplete_intron_retention_5",
    "incomplete_intron_retention_3",
}


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------

def _open(path):
    """Open plain or gzip file for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _find_file(directory, patterns):
    """Return the first existing file matching any of the glob patterns."""
    for pattern in patterns:
        matches = glob.glob(os.path.join(directory, pattern))
        if matches:
            return sorted(matches)[0]
    return None


def find_grouped_counts(isoquant_dir, sample, feature):
    """
    Locate IsoQuant barcode-grouped count files for a given feature type
    ('gene' or 'transcript').

    IsoQuant names these files with the '_grouped_barcode_counts' infix when
    run in single-cell mode (--barcoded_bam / --barcoded_reads). The older
    '_grouped_counts' infix is also tried as a fallback for bulk-mode runs.

    Preference order:
      1. MTX trio: {sample}.{feature}_grouped_barcode_counts.matrix.mtx
                   + .features.tsv + .barcodes.tsv
      2. Linear TSV: {sample}.{feature}_grouped_barcode_counts.linear.tsv
      3. Fallback (bulk/older SC): {sample}.{feature}_grouped_counts.* variants

    Returns ('mtx', mtx_path, features_path, barcodes_path) or
            ('linear', linear_path) or ('tsv', tsv_path) or raises.
    """
    # SC mode uses _grouped_barcode_counts; try that first, then legacy names
    for infix in ("_grouped_barcode_counts", "_grouped_counts"):
        mtx_path = _find_file(
            isoquant_dir,
            [
                "{}.{}{}.matrix.mtx".format(sample, feature, infix),
                "*.{}{}.matrix.mtx".format(feature, infix),
            ],
        )
        if mtx_path:
            feat_path = mtx_path.replace(".matrix.mtx", ".features.tsv")
            bc_path = mtx_path.replace(".matrix.mtx", ".barcodes.tsv")
            if os.path.exists(feat_path) and os.path.exists(bc_path):
                return ("mtx", mtx_path, feat_path, bc_path)

        linear_path = _find_file(
            isoquant_dir,
            [
                "{}.{}{}.linear.tsv".format(sample, feature, infix),
                "*.{}{}.linear.tsv".format(feature, infix),
            ],
        )
        if linear_path:
            return ("linear", linear_path)

        tsv_path = _find_file(
            isoquant_dir,
            [
                "{}.{}{}.tsv".format(sample, feature, infix),
                "*.{}{}.tsv".format(feature, infix),
            ],
        )
        if tsv_path:
            return ("tsv", tsv_path)

    raise FileNotFoundError(
        "No {}_grouped_barcode_counts (or _grouped_counts) file found in "
        "{}".format(feature, isoquant_dir)
    )


def _scan_header(path):
    """
    Locate the tab-separated column header in an IsoQuant assignment file.

    Returns (col_names, skiprows) where skiprows is the number of lines before
    the first data row.
    """
    preview = []
    with _open(path) as fh:
        skip = 0
        for raw_line in fh:
            skip += 1
            if len(preview) < 5:
                preview.append(raw_line.rstrip("\n")[:200])
            parts = raw_line.lstrip("#").strip().split("\t")
            if len(parts) > 1:
                col_names = [p.lstrip("#").strip() for p in parts]
                return col_names, skip

    size = os.path.getsize(path)
    preview_txt = "\n  ".join(preview) if preview else "(file empty)"
    raise ValueError(
        "Could not find a tab-separated header line in: {} ({} bytes). "
        "First lines:\n  {}\n"
        "The file may be empty or from an interrupted IsoQuant run. "
        "Re-run IsoQuant with isoquant.large_output: read_info (or "
        "read_assignments) in pipeline.yml.".format(path, size, preview_txt)
    )


def _assignment_file_usable(path):
    """Return True if path looks like a non-empty IsoQuant assignment TSV."""
    try:
        _scan_header(path)
        return True
    except (ValueError, OSError):
        return False


def find_read_info(isoquant_dir, sample):
    """
    Locate a usable per-read assignment file.

    Newer IsoQuant produces read_info.tsv(.gz) with native barcode/umi cols.
    Older releases produce read_assignments.tsv(.gz); barcodes are recovered
    from the UMI_filtered allinfo file when needed.

    Skips empty or truncated files (common after OOM/interrupted runs).

    Returns (path, format_name) where format_name is 'read_info' or
    'read_assignments'.
    """
    candidates = []
    for fmt, patterns in [
        ("read_info", [
            "{}.read_info.tsv.gz".format(sample),
            "{}.read_info.tsv".format(sample),
            "*.read_info.tsv.gz",
            "*.read_info.tsv",
        ]),
        ("read_assignments", [
            "{}.read_assignments.tsv.gz".format(sample),
            "{}.read_assignments.tsv".format(sample),
            "*.read_assignments.tsv.gz",
            "*.read_assignments.tsv",
        ]),
    ]:
        path = _find_file(isoquant_dir, patterns)
        if path:
            candidates.append((path, fmt))

    for path, fmt in candidates:
        if _assignment_file_usable(path):
            log.info("Found usable per-read file (%s): %s", fmt, path)
            return path, fmt
        log.warning(
            "Skipping unusable per-read file (%s): %s "
            "(empty or missing header — likely incomplete IsoQuant output)",
            fmt, path,
        )

    if candidates:
        paths = ", ".join(p for p, _ in candidates)
        raise FileNotFoundError(
            "Per-read assignment files exist but none are readable in {}: {}. "
            "Re-run IsoQuant with isoquant.large_output: read_info.".format(
                isoquant_dir, paths
            )
        )
    raise FileNotFoundError(
        "No read_info or read_assignments file found in {}".format(isoquant_dir)
    )


def find_allinfo(isoquant_dir):
    """
    Locate the UMI-filtered allinfo file (post-dedup representative reads).

    This file carries barcode and umi columns and is the authoritative source
    for which reads survived UMI deduplication. It is always present when
    IsoQuant runs in single-cell mode (--barcoded_bam / --barcoded_reads).

    Returns path or None if absent.
    """
    path = _find_file(
        isoquant_dir,
        [
            "*.UMI_filtered.ED*.allinfo.gz",
            "*.UMI_filtered.ED*.allinfo",
        ],
    )
    if path:
        log.info("Found allinfo file: %s", path)
    else:
        log.warning("No UMI_filtered allinfo file found in %s. "
                    "Barcode/UMI join will not be available for "
                    "read_assignments format.", isoquant_dir)
    return path


# ---------------------------------------------------------------------------
# Loading native count matrices
# ---------------------------------------------------------------------------

def load_grouped_counts(isoquant_dir, sample, feature):
    """
    Load an IsoQuant barcode-grouped count matrix as a sparse DataFrame.

    Returns (csr_matrix, features_index, barcodes_index).
    """
    result = find_grouped_counts(isoquant_dir, sample, feature)
    fmt = result[0]

    if fmt == "mtx":
        _, mtx_path, feat_path, bc_path = result
        log.info("Loading %s counts from MTX: %s", feature, mtx_path)
        mat = sp.load_npz(mtx_path) if mtx_path.endswith(".npz") else \
            _load_mtx(mtx_path)
        features = pd.read_csv(feat_path, sep="\t", header=None)[0].tolist()
        barcodes = pd.read_csv(bc_path, sep="\t", header=None)[0].tolist()
        return mat.tocsr(), features, barcodes

    if fmt == "linear":
        _, path = result
        log.info("Loading %s counts from linear TSV: %s", feature, path)
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["feature_id", "group_id", "count"],
                         comment="#")
        features = sorted(df["feature_id"].unique())
        barcodes = sorted(df["group_id"].unique())
        feat_idx = {f: i for i, f in enumerate(features)}
        bc_idx = {b: i for i, b in enumerate(barcodes)}
        rows = df["feature_id"].map(feat_idx).values
        cols = df["group_id"].map(bc_idx).values
        data = df["count"].values.astype(np.float32)
        mat = sp.csr_matrix(
            (data, (rows, cols)), shape=(len(features), len(barcodes))
        )
        return mat, features, barcodes

    if fmt == "tsv":
        _, path = result
        log.info("Loading %s counts from wide TSV: %s", feature, path)
        df = pd.read_csv(path, sep="\t", index_col=0)
        features = df.index.tolist()
        barcodes = df.columns.tolist()
        mat = sp.csr_matrix(df.values.astype(np.float32))
        return mat, features, barcodes

    raise RuntimeError("Unknown format: {}".format(fmt))


def _load_mtx(path):
    """Load a Matrix Market file (.mtx) as a scipy CSR matrix."""
    from scipy.io import mmread
    mat = mmread(path)
    return mat.tocsr()


# ---------------------------------------------------------------------------
# Per-read labelling
# ---------------------------------------------------------------------------

def _has_ir(events_str):
    """Return True if any intron-retention token appears in assignment_events."""
    if not events_str or pd.isna(events_str):
        return False
    for token in IR_TOKENS:
        if token in events_str:
            return True
    return False


def load_read_labels(path, fmt, allinfo_path=None):
    """
    Parse the per-read assignment file and return a DataFrame with columns:
      barcode, umi, gene_id, isoform_id, assignment_type, assignment_events,
      label  ('spliced', 'unspliced', 'ambiguous')

    Parameters
    ----------
    path : str
        Path to read_info.tsv(.gz) or read_assignments.tsv(.gz).
    fmt : str
        'read_info' or 'read_assignments'.
    allinfo_path : str or None
        Path to *.UMI_filtered.ED*.allinfo(.gz). When provided, restrict to
        the surviving post-dedup read set before collapsing.
    """
    log.info("Loading per-read assignments (%s): %s", fmt, path)

    col_names, skip = _scan_header(path)
    log.info("Detected %d columns (skipped %d header/comment lines): %s",
             len(col_names), skip, col_names[:6])

    df = pd.read_csv(
        path, sep="\t",
        header=None,
        names=col_names,
        skiprows=skip,
        dtype=str,
        low_memory=False,
    )

    if df.empty:
        raise ValueError(
            "Per-read assignment file has a header but no data rows: {}".format(
                path
            )
        )

    # Normalise column names (strip any residual leading #/whitespace)
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    log.info("Read assignment columns: %s", list(df.columns))

    if fmt == "read_info":
        df = _normalise_read_info(df)
    else:
        # _normalise_read_assignments joins allinfo to add barcode/umi.
        # The allinfo file only contains post-dedup representative reads,
        # so non-surviving reads receive NaN barcode/umi and are dropped
        # naturally by collapse_to_molecules. No second filter needed.
        df = _normalise_read_assignments(df, allinfo_path)

    # Label each read
    df["label"] = "ambiguous"
    spliced_mask = (
        df["assignment_type"].isin(SPLICED_TYPES)
        & ~df["assignment_events"].apply(_has_ir)
    )
    ir_mask = df["assignment_events"].apply(_has_ir)
    df.loc[spliced_mask, "label"] = "spliced"
    df.loc[ir_mask & ~spliced_mask, "label"] = "unspliced"

    n_spliced = spliced_mask.sum()
    n_unspliced = (ir_mask & ~spliced_mask).sum()
    n_ambiguous = len(df) - n_spliced - n_unspliced
    log.info(
        "Label counts (all reads): spliced=%d unspliced=%d ambiguous=%d",
        n_spliced, n_unspliced, n_ambiguous,
    )

    return df


def _normalise_read_info(df):
    """
    Normalise a read_info.tsv DataFrame to a common schema.

    read_info.tsv has native barcode and umi columns.
    Isoform assignment is in isoform_assignment_type or assignment_type.
    """
    rename = {}
    # Barcode / UMI
    for src, dst in [("barcode", "barcode"), ("umi", "umi"),
                     ("gene_id", "gene_id"), ("isoform_id", "isoform_id")]:
        if src in df.columns:
            rename[src] = dst
    # Assignment type: prefer isoform_assignment_type, fall back to gene_assignment_type
    if "isoform_assignment_type" in df.columns:
        rename["isoform_assignment_type"] = "assignment_type"
    elif "assignment_type" in df.columns:
        pass
    if "assignment_events" not in df.columns and "isoform_assignment_events" in df.columns:
        rename["isoform_assignment_events"] = "assignment_events"
    df = df.rename(columns=rename)
    _require_columns(df, ["barcode", "umi", "gene_id", "assignment_type"])
    if "assignment_events" not in df.columns:
        df["assignment_events"] = ""
    return df


def _normalise_read_assignments(df, allinfo_path):
    """
    Normalise a read_assignments.tsv DataFrame to a common schema.

    read_assignments.tsv lacks native barcode/umi columns; extract from
    the allinfo file or the 'additional' column (--bam_tags CB,UB).
    """
    if "assignment_events" not in df.columns:
        df["assignment_events"] = ""

    # Try to get barcode/umi from allinfo
    if allinfo_path:
        allinfo = _load_allinfo_barcodes(allinfo_path)
        if "read_id" in df.columns and not allinfo.empty:
            # Cast both sides to str so the merge key types always match
            df["read_id"] = df["read_id"].astype(str)
            allinfo["read_id"] = allinfo["read_id"].astype(str)
            df = df.merge(allinfo[["read_id", "barcode", "umi"]],
                          on="read_id", how="left")
            log.info("Joined allinfo barcodes: %d/%d reads matched",
                     df["barcode"].notna().sum(), len(df))
            return df

    # Try 'additional' column for CB:Z:/UB:Z: tags
    if "additional" in df.columns:
        df["barcode"] = df["additional"].str.extract(r"CB:Z:([^\t,;]+)")
        df["umi"] = df["additional"].str.extract(r"UB:Z:([^\t,;]+)")
    else:
        log.warning(
            "read_assignments.tsv has no barcode/umi columns and no allinfo "
            "file was found. Barcode/UMI will be missing."
        )
        df["barcode"] = np.nan
        df["umi"] = np.nan

    _require_columns(df, ["gene_id", "assignment_type"])
    return df


def _require_columns(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            "Required columns missing from assignment file: {}. "
            "Available: {}".format(missing, list(df.columns))
        )


# Canonical column order for IsoQuant *.UMI_filtered.ED*.allinfo files.
# These files often have no header row; the first data line looks like a
# multi-field TSV and must not be mistaken for a header.
ALLINFO_COLUMNS = [
    "read_id",
    "gene_id",
    "cell_type",
    "barcode",
    "umi",
    "introns",
    "TSS",
    "polyA",
    "exons",
    "read_type",
    "intron_count",
    "transcript_id",
    "transcript_type",
]


def _scan_tsv_header(fh):
    """
    Advance fh past comment lines, return (col_names, n_skipped).

    Reads lines until a multi-field tab-separated line is found; that line
    is treated as the header (stripping a leading '#' if present).
    """
    for n, raw_line in enumerate(fh, start=1):
        parts = raw_line.lstrip("#").strip().split("\t")
        if len(parts) > 1:
            return parts, n
    raise ValueError("No tab-separated header found in file.")


def _looks_like_allinfo_header(parts):
    """True if the multi-field line is a real header, not a data row."""
    lowered = {p.strip().lower() for p in parts}
    return "read_id" in lowered or "barcode" in lowered


def _load_allinfo_df(path):
    """
    Load an IsoQuant allinfo file as a DataFrame.

    The allinfo file often has no header: after optional comment lines the
    first multi-field line is already data. Detect that case and apply the
    known ALLINFO_COLUMNS schema instead of treating the data row as names.
    """
    with _open(path) as fh:
        first_fields, line_no = _scan_tsv_header(fh)

    if _looks_like_allinfo_header(first_fields):
        col_names = [c.lstrip("#").strip() for c in first_fields]
        skip = line_no  # skip comments + header row
    else:
        # No header: first multi-field line is data. Skip only preceding
        # comment lines so that row is included in the DataFrame.
        n_cols = len(first_fields)
        if n_cols == len(ALLINFO_COLUMNS):
            col_names = list(ALLINFO_COLUMNS)
        elif n_cols >= 5:
            # Minimal SC fields we need; pad any extras generically
            col_names = list(ALLINFO_COLUMNS[:5]) + [
                "col{}".format(i) for i in range(5, n_cols)
            ]
        else:
            raise ValueError(
                "allinfo file has unexpected column count {}: {}".format(
                    n_cols, first_fields[:5]
                )
            )
        skip = line_no - 1  # keep the first data line
        log.info(
            "allinfo has no header row; using schema for %d columns",
            len(col_names),
        )

    log.info("allinfo columns (%d, skiprows=%d): %s",
             len(col_names), skip, col_names[:6])
    df = pd.read_csv(
        path, sep="\t",
        header=None,
        names=col_names,
        skiprows=skip,
        dtype=str,
        low_memory=False,
        comment="#",
    )
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def _load_allinfo_read_ids(path):
    """Return set of surviving read_ids from an allinfo file."""
    log.info("Loading allinfo read IDs from: %s", path)
    df = _load_allinfo_df(path)
    if "read_id" not in df.columns:
        log.warning("No 'read_id' column in allinfo; post-dedup filter skipped.")
        return set()
    return set(df["read_id"].dropna().astype(str))


def _load_allinfo_barcodes(path):
    """Return DataFrame with read_id, barcode, umi from allinfo file."""
    log.info("Loading allinfo barcodes from: %s", path)
    df = _load_allinfo_df(path)
    rename = {}
    for candidate, target in [("read_id", "read_id"),
                               ("barcode", "barcode"),
                               ("umi", "umi"),
                               ("UMI", "umi"),
                               ("cell_barcode", "barcode"),
                               ("CB", "barcode"),
                               ("UB", "umi")]:
        if candidate in df.columns and target not in rename.values():
            rename[candidate] = target
    df = df.rename(columns=rename)
    for col in ["read_id", "barcode", "umi"]:
        if col not in df.columns:
            log.warning("allinfo missing expected column '%s'", col)
            df[col] = np.nan
        else:
            df[col] = df[col].astype(str)
    n_ok = df["barcode"].notna().sum()
    log.info("allinfo barcode table: %d rows, %d with barcode",
             len(df), n_ok)
    return df[["read_id", "barcode", "umi"]]


# ---------------------------------------------------------------------------
# Molecule-level collapse
# ---------------------------------------------------------------------------

def collapse_to_molecules(df, feature_col="gene_id"):
    """
    Collapse read-level labels to molecule-level by (barcode, umi, feature).

    feature_col is typically 'gene_id' (gene-level) or 'isoform_id'
    (transcript-level). Strategy: if any read in the group is labelled
    'spliced', the molecule is spliced; if any is 'unspliced' and none are
    spliced, it is unspliced; otherwise ambiguous.
    """
    log.info(
        "Collapsing %d reads to molecules by (barcode, umi, %s)...",
        len(df), feature_col,
    )
    needed = ["barcode", "umi", feature_col, "label"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(
            "collapse_to_molecules missing columns: {}. Available: {}".format(
                missing, list(df.columns)
            )
        )

    df = df.dropna(subset=["barcode", "umi", feature_col]).copy()
    df = df[df["barcode"] != "."]
    # Drop rows where barcode/umi were never joined (literal 'nan' from cast)
    df = df[~df["barcode"].isin(["nan", "None", ""])]
    df = df[~df["umi"].isin(["nan", "None", ""])]

    def _mol_label(labels):
        s = set(labels)
        if "spliced" in s:
            return "spliced"
        if "unspliced" in s:
            return "unspliced"
        return "ambiguous"

    mol = (
        df.groupby(["barcode", "umi", feature_col], sort=False)["label"]
        .apply(_mol_label)
        .reset_index()
    )
    # Normalise feature column name to gene_id for build_molecule_matrix
    if feature_col != "gene_id":
        mol = mol.rename(columns={feature_col: "gene_id"})
    log.info("Collapsed to %d molecules.", len(mol))
    return mol


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------

def build_molecule_matrix(mol_df, label, features, barcodes, feature_col="gene_id"):
    """
    Build a sparse (features x barcodes) count matrix for reads with a given
    label ('spliced' or 'unspliced').

    Returns scipy.sparse.csr_matrix with shape (len(features), len(barcodes)).
    """
    subset = mol_df[mol_df["label"] == label]
    feat_idx = {f: i for i, f in enumerate(features)}
    bc_idx = {b: i for i, b in enumerate(barcodes)}

    rows, cols, data = [], [], []
    for _, row in subset.iterrows():
        feat = row[feature_col]
        bc = row["barcode"]
        if feat in feat_idx and bc in bc_idx:
            rows.append(feat_idx[feat])
            cols.append(bc_idx[bc])
            data.append(1)

    if not rows:
        return sp.csr_matrix((len(features), len(barcodes)), dtype=np.float32)

    # Use coo_matrix then sum duplicates
    mat = sp.coo_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(len(features), len(barcodes)),
    ).tocsr()
    return mat


# ---------------------------------------------------------------------------
# Reconciliation
# ---------------------------------------------------------------------------

def reconcile(native_total, derived_spliced, derived_unspliced,
              derived_ambiguous, barcodes, tolerance):
    """
    Compare spliced+unspliced+ambiguous totals per barcode against the native
    IsoQuant count matrix. Emit a warning for barcodes that diverge beyond
    tolerance (relative difference).
    """
    n_sp = np.asarray(derived_spliced.sum(axis=0)).flatten()
    n_un = np.asarray(derived_unspliced.sum(axis=0)).flatten()
    n_am = np.asarray(derived_ambiguous.sum(axis=0)).flatten()
    derived_total = n_sp + n_un + n_am

    native = np.asarray(native_total.sum(axis=0)).flatten()

    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = np.where(
            native > 0,
            np.abs(derived_total - native) / native,
            np.where(derived_total > 0, 1.0, 0.0),
        )

    flagged = rel_diff > tolerance
    n_flagged = flagged.sum()
    if n_flagged > 0:
        worst_idx = int(np.argmax(rel_diff))
        log.warning(
            "Reconciliation: %d/%d barcodes exceed tolerance %.4f. "
            "Worst barcode '%s': native=%g derived=%g (rel_diff=%.4f). "
            "Spliced/unspliced labels may not fully account for IsoQuant totals.",
            n_flagged, len(barcodes), tolerance,
            barcodes[worst_idx], native[worst_idx], derived_total[worst_idx],
            rel_diff[worst_idx],
        )
    else:
        log.info(
            "Reconciliation passed: all %d barcodes within tolerance %.4f.",
            len(barcodes), tolerance,
        )


# ---------------------------------------------------------------------------
# AnnData construction and export
# ---------------------------------------------------------------------------

def _common_barcodes(bc_gene, bc_tx):
    """Union of barcodes from gene and transcript matrices."""
    return sorted(set(bc_gene) | set(bc_tx))


def _reindex_matrix(mat, src_features, src_barcodes,
                    tgt_features, tgt_barcodes):
    """Reindex a sparse matrix to new feature/barcode orderings."""
    feat_map = {f: i for i, f in enumerate(src_features)}
    bc_map = {b: i for i, b in enumerate(src_barcodes)}

    rows, cols, data = [], [], []
    mat_coo = mat.tocoo()
    for r, c, v in zip(mat_coo.row, mat_coo.col, mat_coo.data):
        feat = src_features[r]
        bc = src_barcodes[c]
        if feat in feat_map and bc in bc_map:
            ti = tgt_features.index(feat) if feat in tgt_features else None
            bci = tgt_barcodes.index(bc) if bc in tgt_barcodes else None
            if ti is not None and bci is not None:
                rows.append(ti)
                cols.append(bci)
                data.append(v)

    return sp.csr_matrix(
        (data, (rows, cols)) if data else ([], ([], [])),
        shape=(len(tgt_features), len(tgt_barcodes)),
        dtype=np.float32,
    )


def build_gene_anndata(gene_total_mat, gene_features, gene_barcodes,
                       mol_df, tolerance):
    """
    Build a gene-level AnnData with layers: total, spliced, unspliced.

    total   = IsoQuant native gene matrix (molecules, post dedup).
    spliced = molecules assigned spliced (projected to gene).
    unspliced = molecules assigned unspliced.
    """
    log.info("Building gene-level AnnData (%d genes, %d barcodes)...",
             len(gene_features), len(gene_barcodes))

    spliced_mat = build_molecule_matrix(
        mol_df, "spliced", gene_features, gene_barcodes)
    unspliced_mat = build_molecule_matrix(
        mol_df, "unspliced", gene_features, gene_barcodes)
    ambiguous_mat = build_molecule_matrix(
        mol_df, "ambiguous", gene_features, gene_barcodes)

    reconcile(gene_total_mat, spliced_mat, unspliced_mat, ambiguous_mat,
              gene_barcodes, tolerance)

    adata = ad.AnnData(
        X=gene_total_mat.T.tocsr(),
        obs=pd.DataFrame(index=gene_barcodes),
        var=pd.DataFrame(index=gene_features),
    )
    adata.layers["total"] = gene_total_mat.T.tocsr()
    adata.layers["spliced"] = spliced_mat.T.tocsr()
    adata.layers["unspliced"] = unspliced_mat.T.tocsr()
    return adata


def build_gene_anndata_total_only(gene_total_mat, gene_features, gene_barcodes):
    """Build gene AnnData from native counts when per-read labels are unavailable."""
    log.warning(
        "Building total-only gene AnnData (no spliced/unspliced split). "
        "Set isoquant.large_output: read_info and re-run IsoQuant for "
        "splice-stratified layers."
    )
    zero = sp.csr_matrix(
        gene_total_mat.shape, dtype=np.float32,
    )
    adata = ad.AnnData(
        X=gene_total_mat.T.tocsr(),
        obs=pd.DataFrame(index=gene_barcodes),
        var=pd.DataFrame(index=gene_features),
    )
    adata.layers["total"] = gene_total_mat.T.tocsr()
    adata.layers["spliced"] = zero.T.tocsr()
    adata.layers["unspliced"] = zero.T.tocsr()
    return adata


def build_tx_anndata_total_only(tx_total_mat, tx_features, tx_barcodes):
    """Build transcript AnnData from native counts when labels are unavailable."""
    return ad.AnnData(
        X=tx_total_mat.T.tocsr(),
        obs=pd.DataFrame(index=tx_barcodes),
        var=pd.DataFrame(index=tx_features),
    )


def build_barcode_qc_total_only(gene_total_mat, gene_features, gene_barcodes):
    """Per-barcode QC from native gene matrix only (no splice labels)."""
    mat_csc = gene_total_mat.tocsc()
    rows = []
    for i, bc in enumerate(gene_barcodes):
        col_vec = mat_csc.getcol(i)
        total_umis = float(col_vec.sum())
        n_genes = int((col_vec > 0).sum())
        rows.append({
            "barcode": bc,
            "total_umis": total_umis,
            "spliced_umis": 0.0,
            "unspliced_umis": 0.0,
            "ambiguous_umis": 0.0,
            "unspliced_fraction": 0.0,
            "n_genes": n_genes,
        })
    return pd.DataFrame(rows)


def project_tx_to_gene(tx_spliced_mat, tx_features, gene_features,
                       isoquant_dir, sample):
    """
    Project transcript-level spliced counts to gene level using IsoQuant's
    transcript->gene mapping from the GTF/genedb or from read_info.

    Attempts to read a transcript-to-gene map from read_info; if unavailable,
    uses a naive prefix strip (ENST -> ENSG) as a last resort.
    """
    tx_to_gene = _load_tx_gene_map(isoquant_dir, sample)

    gene_list = sorted(gene_features)
    gene_idx = {g: i for i, g in enumerate(gene_list)}

    coo = tx_spliced_mat.tocoo()
    rows, cols, data = [], [], []
    for r, c, v in zip(coo.row, coo.col, coo.data):
        tx = tx_features[r]
        gene = tx_to_gene.get(tx)
        if gene and gene in gene_idx:
            rows.append(gene_idx[gene])
            cols.append(c)
            data.append(v)

    if not rows:
        log.warning("No transcript->gene mappings found; gene-projected "
                    "spliced layer will be empty.")
        return sp.csr_matrix(
            (len(gene_list), tx_spliced_mat.shape[1]), dtype=np.float32
        ), gene_list

    mat = sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(gene_list), tx_spliced_mat.shape[1]),
        dtype=np.float32,
    ).tocsr()
    return mat, gene_list


def _load_tx_gene_map(isoquant_dir, sample):
    """
    Build a transcript -> gene mapping from read_info.tsv or
    read_assignments.tsv (gene_id + isoform_id columns).
    Returns dict {transcript_id: gene_id}.
    """
    path = _find_file(
        isoquant_dir,
        [
            "{}.read_info.tsv.gz".format(sample),
            "{}.read_info.tsv".format(sample),
            "*.read_info.tsv.gz",
            "*.read_info.tsv",
            "{}.read_assignments.tsv.gz".format(sample),
            "{}.read_assignments.tsv".format(sample),
            "*.read_assignments.tsv.gz",
            "*.read_assignments.tsv",
        ],
    )
    if path:
        log.info("Building tx->gene map from: %s", path)
    if not path:
        return {}

    tx_to_gene = {}
    with _open(path) as fh:
        header, _ = _scan_tsv_header(fh)
        gene_col = next((i for i, c in enumerate(header)
                         if c in ("gene_id", "gene")), None)
        iso_col = next((i for i, c in enumerate(header)
                        if c in ("isoform_id", "transcript_id")), None)
        if gene_col is None or iso_col is None:
            return {}
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) > max(gene_col, iso_col):
                tx_to_gene[parts[iso_col]] = parts[gene_col]
    return tx_to_gene


def build_tx_anndata(tx_spliced_mat, tx_features, barcodes):
    """Build transcript-level AnnData for spliced molecules."""
    log.info("Building transcript-level AnnData (%d transcripts, %d barcodes)...",
             len(tx_features), len(barcodes))
    adata = ad.AnnData(
        X=tx_spliced_mat.T.tocsr(),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=tx_features),
    )
    return adata


# ---------------------------------------------------------------------------
# MTX export
# ---------------------------------------------------------------------------

def export_mtx(mat, features, barcodes, directory):
    """Write matrix.mtx, features.tsv, barcodes.tsv to directory."""
    os.makedirs(directory, exist_ok=True)
    from scipy.io import mmwrite
    mmwrite(os.path.join(directory, "matrix.mtx"), mat.T.tocoo())
    pd.Series(features).to_csv(
        os.path.join(directory, "features.tsv"), index=False, header=False)
    pd.Series(barcodes).to_csv(
        os.path.join(directory, "barcodes.tsv"), index=False, header=False)
    log.info("MTX export -> %s", directory)


# ---------------------------------------------------------------------------
# Per-barcode QC table
# ---------------------------------------------------------------------------

def build_barcode_qc(mol_df, gene_total_mat, gene_features, gene_barcodes):
    """
    Build a per-barcode QC table:
      total_umis, spliced_umis, unspliced_umis, ambiguous_umis,
      unspliced_fraction, n_genes.
    """
    bc_counts = (
        mol_df.groupby(["barcode", "label"])
        .size()
        .unstack(fill_value=0)
    )
    for col in ["spliced", "unspliced", "ambiguous"]:
        if col not in bc_counts.columns:
            bc_counts[col] = 0

    bc_counts["total_umis"] = (
        bc_counts["spliced"] + bc_counts["unspliced"] + bc_counts["ambiguous"]
    )
    bc_counts["unspliced_fraction"] = np.where(
        bc_counts["total_umis"] > 0,
        bc_counts["unspliced"] / bc_counts["total_umis"],
        0.0,
    )

    # Count expressed genes per barcode from native total matrix
    bc_idx = {b: i for i, b in enumerate(gene_barcodes)}
    n_genes = {}
    mat_csc = gene_total_mat.tocsc()
    for bc in bc_counts.index:
        if bc in bc_idx:
            col_vec = mat_csc.getcol(bc_idx[bc])
            n_genes[bc] = int((col_vec > 0).sum())
        else:
            n_genes[bc] = 0
    bc_counts["n_genes"] = pd.Series(n_genes)

    qc = bc_counts[
        ["total_umis", "spliced", "unspliced", "ambiguous",
         "unspliced_fraction", "n_genes"]
    ].rename(columns={"spliced": "spliced_umis",
                      "unspliced": "unspliced_umis",
                      "ambiguous": "ambiguous_umis"})
    qc.index.name = "barcode"
    return qc.reset_index()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build spliced/unspliced matrices from IsoQuant SC output."
    )
    parser.add_argument("--isoquant-dir", required=True,
                        help="IsoQuant output directory for the sample.")
    parser.add_argument("--sample", required=True,
                        help="Sample prefix (IsoQuant -p value).")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for matrices.")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Reconciliation tolerance (default 0.01).")
    args = parser.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    # ---- Load native count matrices ----------------------------------------
    gene_result = load_grouped_counts(args.isoquant_dir, args.sample, "gene")
    gene_total_mat, gene_features, gene_barcodes = gene_result

    tx_result = load_grouped_counts(
        args.isoquant_dir, args.sample, "transcript")
    tx_total_mat, tx_features, tx_barcodes = tx_result

    all_barcodes = sorted(set(gene_barcodes) | set(tx_barcodes))

    # ---- Load per-read assignments (optional for splice labelling) ----------
    allinfo_path = find_allinfo(args.isoquant_dir)
    reads_df = None
    try:
        read_path, read_fmt = find_read_info(args.isoquant_dir, args.sample)
        reads_df = load_read_labels(read_path, read_fmt, allinfo_path)
    except (FileNotFoundError, ValueError) as exc:
        log.warning(
            "Per-read assignments unavailable (%s). "
            "Writing total-only matrices from native IsoQuant counts.",
            exc,
        )

    if reads_df is not None:
        mol_df = collapse_to_molecules(reads_df)
        gene_adata = build_gene_anndata(
            gene_total_mat, gene_features, gene_barcodes,
            mol_df, args.tolerance,
        )
    else:
        mol_df = None
        gene_adata = build_gene_anndata_total_only(
            gene_total_mat, gene_features, gene_barcodes,
        )
    gene_h5ad = os.path.join(args.outdir, "{}.gene.h5ad".format(args.sample))
    gene_adata.write_h5ad(gene_h5ad, compression="gzip")
    log.info("Written: %s", gene_h5ad)

    # MTX fallback for gene matrices
    for layer_name, mat in [
        ("total_gene",    gene_total_mat),
        ("spliced_gene",  gene_adata.layers["spliced"].T.tocsr()),
        ("unspliced_gene", gene_adata.layers["unspliced"].T.tocsr()),
    ]:
        export_mtx(
            mat, gene_features, gene_barcodes,
            os.path.join(args.outdir, args.sample, layer_name),
        )

    # ---- Transcript-level AnnData (spliced only) ----------------------------
    if reads_df is not None and "isoform_id" in reads_df.columns:
        mol_tx_df = collapse_to_molecules(reads_df, feature_col="isoform_id")
        tx_spliced_mat = build_molecule_matrix(
            mol_tx_df, "spliced", tx_features, tx_barcodes)
    else:
        if reads_df is not None:
            log.warning(
                "isoform_id column not found in per-read file; "
                "transcript matrix uses native TX totals."
            )
        tx_spliced_mat = tx_total_mat

    if reads_df is not None:
        tx_adata = build_tx_anndata(tx_spliced_mat, tx_features, tx_barcodes)
    else:
        tx_adata = build_tx_anndata_total_only(
            tx_total_mat, tx_features, tx_barcodes,
        )
    tx_h5ad = os.path.join(
        args.outdir, "{}.transcript.spliced.h5ad".format(args.sample))
    tx_adata.write_h5ad(tx_h5ad, compression="gzip")
    log.info("Written: %s", tx_h5ad)

    export_mtx(
        tx_spliced_mat, tx_features, tx_barcodes,
        os.path.join(args.outdir, args.sample, "spliced_tx"),
    )

    # ---- Per-barcode QC table -----------------------------------------------
    if mol_df is not None:
        qc = build_barcode_qc(mol_df, gene_total_mat, gene_features, gene_barcodes)
    else:
        qc = build_barcode_qc_total_only(
            gene_total_mat, gene_features, gene_barcodes,
        )
    qc_path = os.path.join(
        args.outdir, "{}.barcode_qc.tsv".format(args.sample))
    qc.to_csv(qc_path, sep="\t", index=False)
    log.info("Written: %s", qc_path)

    log.info("isoquant_matrices.py complete.")


if __name__ == "__main__":
    sys.exit(main())
