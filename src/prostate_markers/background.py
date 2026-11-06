"""
Slide-background reference.

With one slide per class, nearly every gene differs between the slides because of
donor, processing and sequencing depth. A panel gene is only interesting if it shifts
more than genes of similar expression level do.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def detected(effects: pd.DataFrame, min_detection: float) -> pd.Series:
    return effects[["detect_normal", "detect_cancer"]].max(axis=1) >= min_detection


def expression_bins(mean_expr: pd.Series, n_bins: int) -> pd.Series:
    """Quantile bins of mean expression; labels are 0..n_bins-1."""
    return pd.qcut(mean_expr.rank(method="first"), n_bins, labels=False)


def background_percentiles(
    effects: pd.DataFrame, genes: list[str], min_detection: float, n_bins: int,
) -> pd.DataFrame:
    """
    For each gene, the percentile of its |log2fc| and |d| among all other detected genes
    in the same expression bin. Deterministic: every matched gene is used.
    """
    universe = effects[detected(effects, min_detection)].copy()
    universe["bin"] = expression_bins(universe["mean_log_expr"], n_bins)
    universe = universe.set_index("gene")
    rows = []
    for gene in genes:
        if gene not in universe.index:
            rows.append({"gene": gene, "n_matched": 0,
                         "bg_pct_log2fc": np.nan, "bg_pct_d": np.nan})
            continue
        own = universe.loc[gene]
        matched = universe[(universe["bin"] == own["bin"]) & (universe.index != gene)]
        rows.append({
            "gene": gene,
            "n_matched": len(matched),
            "bg_pct_log2fc": 100 * np.mean(matched["log2fc"].abs() < abs(own["log2fc"])),
            "bg_pct_d": 100 * np.mean(matched["cohens_d"].abs() < abs(own["cohens_d"])),
        })
    return pd.DataFrame(rows)


def separability(y: np.ndarray, score: np.ndarray) -> float:
    """ROC AUC folded around 0.5: how well a score separates the two slides, either way."""
    auc = roc_auc_score(y, score)
    return float(max(auc, 1 - auc))


def random_set_separability(
    X, y: np.ndarray, gene_pool: np.ndarray, set_size: int, n_sets: int, seed: int,
) -> np.ndarray:
    """Separability of the mean expression of `n_sets` random gene sets of `set_size` genes."""
    rng = np.random.default_rng(seed)
    out = np.empty(n_sets)
    for i in range(n_sets):
        idx = np.sort(rng.choice(gene_pool, size=set_size, replace=False))
        out[i] = separability(y, np.asarray(X[:, idx].mean(axis=1)).ravel())
    return out
