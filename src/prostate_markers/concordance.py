"""
Within-slide association between molecular features and H&E image features.

Correlations are computed inside each slide separately, so slide identity cannot drive
them. A link is called concordant when it has the same sign and at least a minimum
strength in both slides. Spots are spatially autocorrelated, so no p-values are given.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def spearman_matrix(A: pd.DataFrame, B: pd.DataFrame) -> pd.DataFrame:
    """Spearman correlation of every column of A with every column of B (same rows)."""
    ra = rankdata(A.to_numpy(float), axis=0)
    rb = rankdata(B.to_numpy(float), axis=0)
    ra = (ra - ra.mean(0)) / ra.std(0)
    rb = (rb - rb.mean(0)) / rb.std(0)
    with np.errstate(invalid="ignore"):
        rho = ra.T @ rb / len(ra)
    return pd.DataFrame(rho, index=A.columns, columns=B.columns)


def within_slide(
    molecular: pd.DataFrame, image: pd.DataFrame, slide: pd.Series, min_abs_rho: float,
) -> pd.DataFrame:
    """Long table: feature x image feature with rho per slide and a concordance flag."""
    ok = image.notna().all(axis=1).to_numpy()
    per_slide = {}
    for name in ("normal", "cancer"):
        rows = ok & (slide == name).to_numpy()
        per_slide[name] = spearman_matrix(molecular[rows], image[rows]).stack()
    out = pd.DataFrame({f"rho_{k}": v for k, v in per_slide.items()})
    out.index.names = ["feature", "image_feature"]
    out = out.reset_index()
    same_sign = np.sign(out["rho_normal"]) == np.sign(out["rho_cancer"])
    out["min_abs_rho"] = out[["rho_normal", "rho_cancer"]].abs().min(axis=1)
    out["concordant"] = same_sign & (out["min_abs_rho"] >= min_abs_rho)
    return out


def strongest_link(table: pd.DataFrame) -> pd.DataFrame:
    """
    Per feature, the concordant image feature with the largest min |rho| (if any).
    `rho` is that weaker of the two correlations, with its shared sign.
    """
    hits = table[table["concordant"]].sort_values("min_abs_rho", ascending=False)
    best = hits.drop_duplicates("feature").copy()
    best["rho"] = np.sign(best["rho_normal"]) * best["min_abs_rho"]
    return best[["feature", "image_feature", "rho_normal", "rho_cancer", "rho"]]
