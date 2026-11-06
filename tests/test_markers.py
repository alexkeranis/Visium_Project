"""Unit tests for the statistical helpers (no raw data needed)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from prostate_markers import background, concordance, effects, normalize, panel


def _slide(counts: np.ndarray, genes: list[str]) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(counts.astype(np.float32)),
                       var=pd.DataFrame(index=genes))
    adata.obs_names = [f"s{i}" for i in range(adata.n_obs)]
    adata.obs["total_counts"] = counts.sum(axis=1)
    return normalize.log_normalize(adata, target_sum=1e4)


# --- normalization ------------------------------------------------------------

def test_log_normalize_keeps_counts_and_scales_rows():
    counts = np.array([[1, 3], [10, 30]])
    adata = _slide(counts, ["a", "b"])
    np.testing.assert_array_equal(adata.layers["counts"].toarray(), counts)
    # Same composition, different depth -> identical normalized rows
    np.testing.assert_allclose(adata.X.toarray()[0], adata.X.toarray()[1], rtol=1e-6)
    np.testing.assert_allclose(np.expm1(adata.X.toarray()).sum(axis=1), 1e4, rtol=1e-5)


def test_housekeeping_helpers():
    assert normalize.housekeeping_genes(["ACTB", "B2M", "X"], ["ACTB", "GAPDH", "B2M"]) == [
        "ACTB", "B2M"]
    scale = normalize.housekeeping_scale(np.array([[2, 4], [0, 0]]))
    assert scale[0] == 3 and np.isnan(scale[1])


# --- effect sizes -------------------------------------------------------------

def test_gene_effects_log2fc_uses_pseudobulk():
    normal = _slide(np.array([[10, 90], [10, 90]]), ["g1", "g2"])
    cancer = _slide(np.array([[40, 60], [40, 60]]), ["g1", "g2"])
    fx = effects.gene_effects(normal, cancer, pseudocount=0).set_index("gene")
    assert fx.loc["g1", "log2fc"] == pytest.approx(2.0)  # 10% -> 40% of counts
    assert fx.loc["g2", "log2fc"] == pytest.approx(np.log2(60 / 90))


def test_prob_superiority_bounds():
    assert effects.prob_superiority(np.array([3, 4, 5]), np.array([0, 1, 2])) == 1.0
    assert effects.prob_superiority(np.array([1, 1]), np.array([1, 1])) == 0.5


# --- slide background ---------------------------------------------------------

def _effects_table(n=100, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "gene": [f"g{i}" for i in range(n)],
        "mean_log_expr": np.linspace(0.1, 3, n),
        "detect_normal": 0.5, "detect_cancer": 0.5,
        "log2fc": rng.normal(0, 0.5, n),
        "cohens_d": rng.normal(0, 0.5, n),
    })


def test_background_percentile_extreme_gene_is_top():
    fx = _effects_table()
    fx.loc[50, ["log2fc", "cohens_d"]] = 10.0
    out = background.background_percentiles(fx, ["g50"], min_detection=0.1, n_bins=5)
    assert out.loc[0, "bg_pct_log2fc"] == 100.0
    assert out.loc[0, "n_matched"] == 19  # 20 genes per bin, minus the gene itself


def test_background_ignores_undetected_genes():
    fx = _effects_table()
    fx.loc[3, ["detect_normal", "detect_cancer"]] = 0.0
    out = background.background_percentiles(fx, ["g3"], min_detection=0.1, n_bins=5)
    assert out.loc[0, "n_matched"] == 0 and np.isnan(out.loc[0, "bg_pct_log2fc"])


def test_random_sets_are_seed_deterministic():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(40, 30))
    y = np.repeat([0, 1], 20)
    pool = np.arange(30)
    a = background.random_set_separability(X, y, pool, 3, 20, seed=7)
    b = background.random_set_separability(X, y, pool, 3, 20, seed=7)
    c = background.random_set_separability(X, y, pool, 3, 20, seed=8)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)
    assert ((a >= 0.5) & (a <= 1)).all()


# --- literature and tiers -----------------------------------------------------

@pytest.mark.parametrize("expected, log2fc, result", [
    ("up", 1.0, "concordant"), ("up", -1.0, "discordant"),
    ("down", -0.5, "concordant"), ("down", 0.5, "discordant"),
    ("context", 2.0, "none"), ("unknown", -2.0, "none"), ("up", np.nan, "none"),
])
def test_literature_agreement(expected, log2fc, result):
    assert panel.literature_agreement(expected, log2fc) == result


def test_assign_tier():
    assert panel.assign_tier(True, True, True, True) == "A"
    assert panel.assign_tier(True, False, True, True) == "B"
    assert panel.assign_tier(True, False, False, True) == "-"
    assert panel.assign_tier(True, True, True, False) == "-"


def test_directional_sets_respect_probe_set():
    pt = pd.DataFrame({"gene": ["A", "B", "C", "D"],
                       "expected": ["up", "up", "down", "context"]})
    assert panel.directional_sets(pt, ["A", "C", "D"]) == {"up": ["A"], "down": ["C"]}


# --- within-slide concordance ---------------------------------------------------

def test_spearman_matrix_matches_pandas():
    rng = np.random.default_rng(3)
    A = pd.DataFrame(rng.normal(size=(50, 2)), columns=["a1", "a2"])
    B = pd.DataFrame(rng.normal(size=(50, 3)), columns=["b1", "b2", "b3"])
    rho = concordance.spearman_matrix(A, B)
    assert rho.loc["a2", "b3"] == pytest.approx(A["a2"].corr(B["b3"], method="spearman"))


def test_within_slide_requires_same_sign_in_both_slides():
    n = 60
    x = np.arange(2 * n, dtype=float)
    slide = pd.Series(["normal"] * n + ["cancer"] * n)
    molecular = pd.DataFrame({"m": x})
    image = pd.DataFrame({
        "same": x,                                          # +1 in both slides
        "flip": np.concatenate([x[:n], -x[n:]]),            # +1 normal, -1 cancer
    })
    out = concordance.within_slide(molecular, image, slide, min_abs_rho=0.2).set_index(
        "image_feature")
    assert out.loc["same", "concordant"]
    assert not out.loc["flip", "concordant"]
    best = concordance.strongest_link(out.reset_index())
    assert best["image_feature"].tolist() == ["same"] and best["rho"].iloc[0] == pytest.approx(1)
