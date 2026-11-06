"""
Per-spot H&E image features from the hires image.

Two families, both computed on a square crop centered on each spot:
- stain density: mean hematoxylin and eosin optical density after color deconvolution
  (hematoxylin ~ nuclei, eosin ~ cytoplasm and stroma);
- gray-level co-occurrence (Haralick-type) statistics, averaged over distances and angles.

These values depend on staining, scanning and quantization, so they are compared
within a slide, not between slides.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from skimage.color import rgb2gray, rgb2hed
from skimage.feature import graycomatrix, graycoprops

GLCM_PROPS = ("contrast", "dissimilarity", "homogeneity", "energy", "correlation")


def hires_image(adata: ad.AnnData) -> tuple[np.ndarray, float, float]:
    """Return (RGB float image in [0, 1], hires scale factor, spot diameter in hires px)."""
    library = next(iter(adata.uns["spatial"]))
    spatial = adata.uns["spatial"][library]
    image = np.asarray(spatial["images"]["hires"])[..., :3]
    if image.dtype == np.uint8:
        image = image / 255.0
    scale = spatial["scalefactors"]["tissue_hires_scalef"]
    diameter = spatial["scalefactors"]["spot_diameter_fullres"] * scale
    return image, scale, diameter


def crop(image: np.ndarray, x: float, y: float, side: int) -> np.ndarray:
    half = side // 2
    r0, c0 = max(int(round(y)) - half, 0), max(int(round(x)) - half, 0)
    return image[r0:r0 + side, c0:c0 + side]


def glcm_features(gray: np.ndarray, levels: int, distances, angles_rad) -> dict[str, float]:
    q = np.clip((gray * levels).astype(np.uint8), 0, levels - 1)
    glcm = graycomatrix(q, distances, angles_rad, levels=levels, symmetric=True, normed=True)
    out = {}
    for prop in GLCM_PROPS:
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = graycoprops(glcm, prop)
        out[f"glcm_{prop}"] = float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan
    return out


def spot_features(adata: ad.AnnData, params: dict) -> pd.DataFrame:
    """One row per spot (index = obs_names) with stain-density and GLCM features."""
    image, scale, diameter = hires_image(adata)
    side = max(int(round(params["spot_scale"] * diameter)), 5)
    angles = np.deg2rad(params["angles_deg"])
    hed = rgb2hed(image)
    gray = rgb2gray(image)
    coords = adata.obsm["spatial"] * scale
    rows = []
    for x, y in coords:
        g = crop(gray, x, y, side)
        h = crop(hed, x, y, side)
        row = {
            "hematoxylin_od": float(h[..., 0].mean()),
            "eosin_od": float(h[..., 1].mean()),
            "gray_mean": float(g.mean()),
            "gray_sd": float(g.std()),
        }
        row.update(glcm_features(g, params["gray_levels"], params["distances"], angles))
        rows.append(row)
    return pd.DataFrame(rows, index=adata.obs_names)
