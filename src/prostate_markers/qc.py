"""Spot-level quality control."""

from __future__ import annotations

import anndata as ad
import pandas as pd
import scanpy as sc


def annotate(adata: ad.AnnData, qc: dict) -> ad.AnnData:
    """Add total_counts, n_genes_by_counts and a boolean `qc_pass` column to .obs."""
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    obs = adata.obs
    obs["qc_pass"] = (
        (obs["total_counts"] >= qc["min_counts"])
        & (obs["total_counts"] <= qc["max_counts"])
        & (obs["n_genes_by_counts"] >= qc["min_genes"])
    )
    return adata


def spot_table(adata: ad.AnnData) -> pd.DataFrame:
    """Per-spot QC values with hires-image coordinates, before filtering."""
    scale = hires_scale(adata)
    return pd.DataFrame({
        "slide": adata.obs["slide"].to_numpy(),
        "barcode": adata.obs_names.to_numpy(),
        "total_counts": adata.obs["total_counts"].to_numpy(),
        "n_genes": adata.obs["n_genes_by_counts"].to_numpy(),
        "qc_pass": adata.obs["qc_pass"].to_numpy(),
        "x_hires": adata.obsm["spatial"][:, 0] * scale,
        "y_hires": adata.obsm["spatial"][:, 1] * scale,
    })


def summary(spots: pd.DataFrame) -> pd.DataFrame:
    """One row per slide: spot counts before/after filtering and median depth."""
    rows = []
    for slide, df in spots.groupby("slide", sort=False):
        kept = df[df["qc_pass"]]
        rows.append({
            "slide": slide,
            "spots_in_tissue": len(df),
            "spots_kept": len(kept),
            "median_counts": kept["total_counts"].median(),
            "median_genes": kept["n_genes"].median(),
        })
    return pd.DataFrame(rows)


def hires_scale(adata: ad.AnnData) -> float:
    library = next(iter(adata.uns["spatial"]))
    return adata.uns["spatial"][library]["scalefactors"]["tissue_hires_scalef"]
