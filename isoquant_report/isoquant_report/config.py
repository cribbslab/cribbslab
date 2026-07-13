"""Load and validate report configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class ReportConfig:
    """Paths and QC defaults for one sample."""

    sample: str = ""
    isoquant_dir: str = ""
    gene_h5ad: str = ""
    transcript_h5ad: str = ""
    barcode_qc: str = ""
    mito_prefix: str = "MT-"
    ribo_prefixes: list[str] = field(default_factory=lambda: ["RPS", "RPL"])
    min_counts: int = 500
    min_genes: int = 200
    max_mito: float = 0.2
    max_unspliced_fraction: float = 0.95

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReportConfig":
        """Build config from a plain dict (YAML or sidebar)."""
        known = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in known}
        return cls(**filtered)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ReportConfig":
        """Load configuration from a YAML file."""
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        return cls.from_dict(data)

    def to_dict(self) -> dict[str, Any]:
        """Serialise config for writing config.yaml."""
        return {
            "sample": self.sample,
            "isoquant_dir": self.isoquant_dir,
            "gene_h5ad": self.gene_h5ad,
            "transcript_h5ad": self.transcript_h5ad,
            "barcode_qc": self.barcode_qc,
            "mito_prefix": self.mito_prefix,
            "ribo_prefixes": self.ribo_prefixes,
            "min_counts": self.min_counts,
            "min_genes": self.min_genes,
            "max_mito": self.max_mito,
            "max_unspliced_fraction": self.max_unspliced_fraction,
        }

    def write_yaml(self, path: str | Path) -> None:
        """Write configuration to YAML."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            yaml.safe_dump(self.to_dict(), fh, default_flow_style=False)
