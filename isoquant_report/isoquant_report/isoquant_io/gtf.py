"""Parse IsoQuant transcript model GTF files."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from typing import Optional

import pandas as pd


@dataclass
class TranscriptModel:
    """One transcript model from IsoQuant GTF."""

    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    exon_count: int
    length: int
    transcript_type: str = "unknown"


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_attrs(field: str) -> dict[str, str]:
    out = {}
    for part in re.findall(r'(\w+)\s+"([^"]*)"', field):
        out[part[0]] = part[1]
    return out


def parse_transcript_models(gtf_path: str) -> pd.DataFrame:
    """
    Parse transcript exon lines from an IsoQuant transcript_models.gtf.

    Returns a DataFrame with one row per transcript (aggregated exons).
    """
    if not gtf_path:
        return pd.DataFrame()

    exons: dict[str, list[tuple]] = {}
    meta: dict[str, dict] = {}

    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom, start, end, strand = (
                parts[0],
                int(parts[3]),
                int(parts[4]),
                parts[6],
            )
            attrs = _parse_attrs(parts[8])
            tx_id = attrs.get("transcript_id", "")
            gene_id = attrs.get("gene_id", attrs.get("gene", ""))
            if not tx_id:
                continue
            exons.setdefault(tx_id, []).append((start, end))
            meta[tx_id] = {
                "transcript_id": tx_id,
                "gene_id": gene_id,
                "chrom": chrom,
                "strand": strand,
                "transcript_type": attrs.get("transcript_type", "unknown"),
            }

    rows = []
    for tx_id, coords in exons.items():
        m = meta[tx_id]
        length = sum(e - s + 1 for s, e in coords)
        rows.append(
            {
                **m,
                "exon_count": len(coords),
                "length": length,
                "start": min(s for s, _ in coords),
                "end": max(e for _, e in coords),
            }
        )
    return pd.DataFrame(rows)
