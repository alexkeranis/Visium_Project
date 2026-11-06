"""Shared plotting helpers. Functions draw on a given Axes; callers save figures."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

SLIDE_COLORS = {"normal": "#4C78A8", "cancer": "#E45756"}


def set_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "font.size": 10,
    })


def save(fig: plt.Figure, folder: Path, name: str) -> None:
    fig.savefig(folder / name, metadata={"Software": None})


def hires_gray(adata: ad.AnnData) -> np.ndarray:
    library = next(iter(adata.uns["spatial"]))
    image = np.asarray(adata.uns["spatial"][library]["images"]["hires"])[..., :3]
    return image @ np.array([0.299, 0.587, 0.114])


def hires_coords(adata: ad.AnnData) -> np.ndarray:
    library = next(iter(adata.uns["spatial"]))
    scale = adata.uns["spatial"][library]["scalefactors"]["tissue_hires_scalef"]
    return adata.obsm["spatial"] * scale


def tissue_overlay(
    ax: plt.Axes, adata: ad.AnnData, values: np.ndarray, *,
    cmap: str = "viridis", vmin: float | None = None, vmax: float | None = None,
    size: float = 6, title: str = "",
):
    """Per-spot values over the grayscale H&E image, cropped to the tissue. Returns the scatter."""
    coords = hires_coords(adata)
    ax.imshow(hires_gray(adata), cmap="gray")
    order = np.argsort(values, kind="stable")  # draw high values last
    sc = ax.scatter(coords[order, 0], coords[order, 1], c=np.asarray(values)[order],
                    s=size, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
    pad = 40
    ax.set_xlim(coords[:, 0].min() - pad, coords[:, 0].max() + pad)
    ax.set_ylim(coords[:, 1].max() + pad, coords[:, 1].min() - pad)
    ax.set_title(title)
    ax.axis("off")
    return sc


def paired_overlay(slides: dict, values: dict, label: str, cmap: str = "viridis",
                   shared_scale: bool = True) -> plt.Figure:
    """Normal and cancer side by side; one color scale unless `shared_scale` is False."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    lo = min(np.nanpercentile(v, 1) for v in values.values())
    hi = max(np.nanpercentile(v, 99) for v in values.values())
    for ax, slide in zip(axes, ("normal", "cancer")):
        v = values[slide]
        vmin, vmax = (lo, hi) if shared_scale else np.nanpercentile(v, [1, 99])
        sc = tissue_overlay(ax, slides[slide], v, cmap=cmap, vmin=vmin, vmax=vmax, title=slide)
        if not shared_scale:
            fig.colorbar(sc, ax=ax, fraction=0.035, label=label)
    if shared_scale:
        fig.colorbar(sc, ax=axes, fraction=0.025, label=label)
    return fig
