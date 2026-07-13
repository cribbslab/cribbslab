"""Discover IsoQuant output files in a sample directory."""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def _find_file(directory: str, patterns: list[str]) -> Optional[str]:
    """Return first existing path matching any glob pattern."""
    for pattern in patterns:
        matches = glob.glob(os.path.join(directory, pattern))
        if matches:
            return sorted(matches)[0]
    return None


@dataclass
class IsoQuantPaths:
    """Resolved paths to IsoQuant artefacts for one sample."""

    isoquant_dir: str
    sample: str
    gene_grouped_mtx: Optional[str] = None
    gene_grouped_features: Optional[str] = None
    gene_grouped_barcodes: Optional[str] = None
    gene_grouped_linear: Optional[str] = None
    gene_grouped_tpm_mtx: Optional[str] = None
    transcript_grouped_mtx: Optional[str] = None
    transcript_grouped_features: Optional[str] = None
    transcript_grouped_barcodes: Optional[str] = None
    transcript_grouped_linear: Optional[str] = None
    transcript_grouped_tpm_mtx: Optional[str] = None
    gene_counts: Optional[str] = None
    transcript_counts: Optional[str] = None
    read_assignments: Optional[str] = None
    read_info: Optional[str] = None
    allinfo: Optional[str] = None
    transcript_models_gtf: Optional[str] = None
    extended_annotation_gtf: Optional[str] = None
    missing: list[str] = field(default_factory=list)

    def completeness(self) -> dict[str, bool]:
        """Return which expected artefacts were found."""
        return {
            "gene_grouped_counts": bool(
                self.gene_grouped_mtx or self.gene_grouped_linear
            ),
            "transcript_grouped_counts": bool(
                self.transcript_grouped_mtx or self.transcript_grouped_linear
            ),
            "read_assignments": bool(
                self.read_assignments or self.read_info
            ),
            "allinfo": bool(self.allinfo),
            "transcript_models_gtf": bool(self.transcript_models_gtf),
            "gene_counts": bool(self.gene_counts),
            "transcript_counts": bool(self.transcript_counts),
        }


def _grouped_paths(
    isoquant_dir: str,
    sample: str,
    feature: str,
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Locate grouped count MTX trio or linear TSV for gene/transcript."""
    for infix in ("_grouped_barcode_counts", "_grouped_counts"):
        mtx = _find_file(
            isoquant_dir,
            [
                f"{sample}.{feature}{infix}.matrix.mtx",
                f"*.{feature}{infix}.matrix.mtx",
            ],
        )
        if mtx:
            feat = mtx.replace(".matrix.mtx", ".features.tsv")
            bc = mtx.replace(".matrix.mtx", ".barcodes.tsv")
            if os.path.exists(feat) and os.path.exists(bc):
                return mtx, feat, bc, None
        linear = _find_file(
            isoquant_dir,
            [
                f"{sample}.{feature}{infix}.linear.tsv",
                f"*.{feature}{infix}.linear.tsv",
            ],
        )
        if linear:
            return None, None, None, linear
    return None, None, None, None


def discover_isoquant(isoquant_dir: str, sample: str) -> IsoQuantPaths:
    """
    Discover IsoQuant SC output files adaptively.

    Parameters
    ----------
    isoquant_dir
        Directory containing IsoQuant output (often isoquant/{sample}/{sample}).
    sample
        Sample prefix used in filenames.
    """
    paths = IsoQuantPaths(isoquant_dir=isoquant_dir, sample=sample)

    g_mtx, g_feat, g_bc, g_lin = _grouped_paths(isoquant_dir, sample, "gene")
    paths.gene_grouped_mtx = g_mtx
    paths.gene_grouped_features = g_feat
    paths.gene_grouped_barcodes = g_bc
    paths.gene_grouped_linear = g_lin

    t_mtx, t_feat, t_bc, t_lin = _grouped_paths(
        isoquant_dir, sample, "transcript"
    )
    paths.transcript_grouped_mtx = t_mtx
    paths.transcript_grouped_features = t_feat
    paths.transcript_grouped_barcodes = t_bc
    paths.transcript_grouped_linear = t_lin

    for infix in ("_grouped_barcode_tpm", "_grouped_tpm"):
        paths.gene_grouped_tpm_mtx = _find_file(
            isoquant_dir,
            [
                f"{sample}.gene{infix}.matrix.mtx",
                f"*.gene{infix}.matrix.mtx",
            ],
        )
        if paths.gene_grouped_tpm_mtx:
            break

    for infix in ("_grouped_barcode_tpm", "_grouped_tpm"):
        paths.transcript_grouped_tpm_mtx = _find_file(
            isoquant_dir,
            [
                f"{sample}.transcript{infix}.matrix.mtx",
                f"*.transcript{infix}.matrix.mtx",
            ],
        )
        if paths.transcript_grouped_tpm_mtx:
            break

    paths.gene_counts = _find_file(
        isoquant_dir,
        [f"{sample}.gene_counts.tsv", "*.gene_counts.tsv"],
    )
    paths.transcript_counts = _find_file(
        isoquant_dir,
        [f"{sample}.transcript_counts.tsv", "*.transcript_counts.tsv"],
    )

    paths.read_info = _find_file(
        isoquant_dir,
        [
            f"{sample}.read_info.tsv.gz",
            f"{sample}.read_info.tsv",
            "*.read_info.tsv.gz",
            "*.read_info.tsv",
        ],
    )
    paths.read_assignments = _find_file(
        isoquant_dir,
        [
            f"{sample}.read_assignments.tsv.gz",
            f"{sample}.read_assignments.tsv",
            "*.read_assignments.tsv.gz",
            "*.read_assignments.tsv",
        ],
    )
    paths.allinfo = _find_file(
        isoquant_dir,
        ["*.UMI_filtered.ED*.allinfo.gz", "*.UMI_filtered.ED*.allinfo"],
    )
    paths.transcript_models_gtf = _find_file(
        isoquant_dir,
        [f"{sample}.transcript_models.gtf", "*.transcript_models.gtf"],
    )
    paths.extended_annotation_gtf = _find_file(
        isoquant_dir,
        [f"{sample}.extended_annotation.gtf", "*.extended_annotation.gtf"],
    )

    comp = paths.completeness()
    for key, found in comp.items():
        if not found:
            paths.missing.append(key)

    return paths
