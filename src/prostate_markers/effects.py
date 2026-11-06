"""Between-slide effect sizes for genes and spot-level features."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu


def _column_mean_var(X) -> tuple[np.ndarray, np.ndarray]:
    """Column means and unbiased variances for dense or sparse X."""
    n = X.shape[0]
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = X.multiply(X) if sp.issparse(X) else np.square(X)
    mean_sq = np.asarray(sq.mean(axis=0)).ravel()
    var = (mean_sq - mean**2) * n / (n - 1)
    return mean, np.clip(var, 0, None)


def _detection(counts) -> np.ndarray:
    return np.asarray((counts > 0).mean(axis=0)).ravel()


def _pseudobulk_cp10k(counts) -> np.ndarray:
    """Summed counts per gene divided by all counts in the slide, per 10k."""
    gene_sums = np.asarray(counts.sum(axis=0)).ravel()
    return gene_sums / gene_sums.sum() * 1e4


def cohens_d(mean1, var1, n1, mean0, var0, n0) -> np.ndarray:
    pooled = np.sqrt(((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(pooled > 0, (mean1 - mean0) / pooled, 0.0)


def gene_effects(normal: ad.AnnData, cancer: ad.AnnData, pseudocount: float) -> pd.DataFrame:
    """
    Per-gene comparison of the two slides on the full transcriptome.

    log2fc: log2 ratio of pseudo-bulk CP10k (cancer / normal).
    cohens_d: standardized difference of spot-level log1p(CP10k).
    """
    if not normal.var_names.equals(cancer.var_names):
        raise ValueError("Both slides must share the same gene order.")
    m0, v0 = _column_mean_var(normal.X)
    m1, v1 = _column_mean_var(cancer.X)
    pb0 = _pseudobulk_cp10k(normal.layers["counts"])
    pb1 = _pseudobulk_cp10k(cancer.layers["counts"])
    n0, n1 = normal.n_obs, cancer.n_obs
    return pd.DataFrame({
        "gene": normal.var_names.to_numpy(),
        "mean_log_expr": (m0 * n0 + m1 * n1) / (n0 + n1),
        "detect_normal": _detection(normal.layers["counts"]),
        "detect_cancer": _detection(cancer.layers["counts"]),
        "cp10k_normal": pb0,
        "cp10k_cancer": pb1,
        "log2fc": np.log2((pb1 + pseudocount) / (pb0 + pseudocount)),
        "cohens_d": cohens_d(m1, v1, n1, m0, v0, n0),
    })


def prob_superiority(x_cancer: np.ndarray, x_normal: np.ndarray) -> float:
    """P(X_cancer > X_normal) + 0.5 P(tie): the Mann-Whitney U scaled to [0, 1]."""
    u = mannwhitneyu(x_cancer, x_normal, alternative="two-sided").statistic
    return float(u / (len(x_cancer) * len(x_normal)))


def feature_effects(values: pd.DataFrame, slide: pd.Series) -> pd.DataFrame:
    """Cohen's d and probability of superiority for each column of a spot x feature table."""
    is_cancer = (slide == "cancer").to_numpy()
    rows = []
    for name in values.columns:
        x = values[name].to_numpy(dtype=float)
        x1, x0 = x[is_cancer], x[~is_cancer]
        d = cohens_d(x1.mean(), x1.var(ddof=1), len(x1), x0.mean(), x0.var(ddof=1), len(x0))
        rows.append({
            "feature": name,
            "mean_normal": x0.mean(),
            "mean_cancer": x1.mean(),
            "cohens_d": float(d),
            "prob_superiority": prob_superiority(x1, x0),
        })
    return pd.DataFrame(rows)
