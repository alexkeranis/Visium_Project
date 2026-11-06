"""Library-size normalization and housekeeping scaling."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scanpy as sc


def log_normalize(adata: ad.AnnData, target_sum: float) -> ad.AnnData:
    """Keep raw counts in layers['counts']; set X to log1p(counts per `target_sum`)."""
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def housekeeping_genes(var_names, configured: list[str]) -> list[str]:
    """Configured housekeeping genes that exist in the probe set, in config order."""
    present = set(var_names)
    return [g for g in configured if g in present]


def housekeeping_scale(counts: np.ndarray) -> np.ndarray:
    """Per-spot mean count of the housekeeping genes (rows = spots). Zeros become NaN."""
    scale = np.asarray(counts, dtype=float).mean(axis=1)
    return np.where(scale > 0, scale, np.nan)

