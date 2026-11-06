"""Literature panel, gene-set scores, log-ratios, normalization sensitivity and tiers."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .effects import cohens_d
from .normalize import housekeeping_genes, housekeeping_scale


def panel_table(config: dict) -> pd.DataFrame:
    df = pd.DataFrame(config["panel"])
    df["pmid"] = df.get("pmid", pd.Series(dtype=float)).astype("Int64")
    return df[["gene", "program", "expected", "pmid"]]


def directional_sets(panel: pd.DataFrame, var_names) -> dict[str, list[str]]:
    """Genes reported up or down in cancer, restricted to the probe set."""
    present = set(var_names)
    return {
        direction: [g for g in panel.loc[panel["expected"] == direction, "gene"] if g in present]
        for direction in ("up", "down")
    }


def literature_agreement(expected: str, log2fc: float) -> str:
    """'concordant' / 'discordant' for up/down expectations, otherwise 'none'."""
    if expected not in ("up", "down") or np.isnan(log2fc) or log2fc == 0:
        return "none"
    observed = "up" if log2fc > 0 else "down"
    return "concordant" if observed == expected else "discordant"


def add_set_scores(combined: ad.AnnData, sets: dict[str, list[str]], seed: int) -> list[str]:
    """
    scanpy score_genes (mean minus expression-matched control genes) for each set.
    Scores are computed on both slides together so they share one control set.
    """
    names = []
    for direction, genes in sets.items():
        if not genes:
            continue
        name = f"score_{direction}"
        sc.tl.score_genes(combined, genes, score_name=name, random_state=seed)
        names.append(name)
    return names


def add_log_ratios(combined: ad.AnnData, pairs: list[list[str]]) -> list[str]:
    """log1p(CP10k a) - log1p(CP10k b) per spot; X already holds log1p(CP10k)."""
    names = []
    for a, b in pairs:
        if a in combined.var_names and b in combined.var_names:
            name = f"logratio_{a}_{b}"
            combined.obs[name] = _gene(combined, a) - _gene(combined, b)
            names.append(name)
    return names


def _gene(adata: ad.AnnData, gene: str, layer: str | None = None) -> np.ndarray:
    col = adata[:, gene].layers[layer] if layer else adata[:, gene].X
    return np.asarray(col.todense() if hasattr(col, "todense") else col).ravel()


def housekeeping_check(
    normal: ad.AnnData, cancer: ad.AnnData, configured: list[str],
) -> pd.DataFrame:
    """
    How suitable the available housekeeping genes are as a reference: their own shift
    between slides, their spot-to-spot variability, and how much they track depth.
    """
    rows = []
    for gene in configured:
        row = {"gene": gene, "in_probe_set": gene in normal.var_names}
        if row["in_probe_set"]:
            for slide, adata in (("normal", normal), ("cancer", cancer)):
                counts = _gene(adata, gene, "counts")
                cp10k = counts / adata.obs["total_counts"].to_numpy() * 1e4
                row[f"cv_cp10k_{slide}"] = cp10k.std(ddof=1) / cp10k.mean()
                row[f"rho_depth_{slide}"] = pd.Series(counts).corr(
                    adata.obs["total_counts"].reset_index(drop=True), method="spearman")
            n_cp = _gene(normal, gene, "counts").sum() / normal.obs["total_counts"].sum()
            c_cp = _gene(cancer, gene, "counts").sum() / cancer.obs["total_counts"].sum()
            row["log2fc"] = np.log2(c_cp / n_cp)
        rows.append(row)
    return pd.DataFrame(rows)


def housekeeping_effects(
    normal: ad.AnnData, cancer: ad.AnnData, genes: list[str], configured_hk: list[str],
) -> pd.DataFrame:
    """Effect sizes when counts are scaled by the housekeeping mean instead of library size."""
    hk = housekeeping_genes(normal.var_names, configured_hk)
    per_slide = {}
    for slide, adata in (("normal", normal), ("cancer", cancer)):
        counts = adata[:, genes].layers["counts"]
        counts = np.asarray(counts.todense() if hasattr(counts, "todense") else counts, float)
        hk_counts = adata[:, hk].layers["counts"]
        hk_counts = np.asarray(hk_counts.todense() if hasattr(hk_counts, "todense") else hk_counts)
        scale = housekeeping_scale(hk_counts)
        ok = ~np.isnan(scale)
        per_slide[slide] = (counts[ok], scale[ok])
    ref = np.median(np.concatenate([s for _, s in per_slide.values()]))
    (cn, sn), (cc, sc_) = per_slide["normal"], per_slide["cancer"]
    log_n = np.log1p(cn / sn[:, None] * ref)
    log_c = np.log1p(cc / sc_[:, None] * ref)
    rate_n = cn.sum(axis=0) / sn.sum()
    rate_c = cc.sum(axis=0) / sc_.sum()
    eps = 1e-6
    return pd.DataFrame({
        "gene": genes,
        "log2fc_hk": np.log2((rate_c + eps) / (rate_n + eps)),
        "cohens_d_hk": cohens_d(
            log_c.mean(0), log_c.var(0, ddof=1), len(log_c),
            log_n.mean(0), log_n.var(0, ddof=1), len(log_n),
        ),
    })


def epithelial_mask(
    adata: ad.AnnData, epithelial_genes: list[str], top_fraction: float,
) -> np.ndarray:
    """Spots in the top `top_fraction` of this slide by mean log-expression of the genes."""
    genes = [g for g in epithelial_genes if g in adata.var_names]
    score = np.asarray(adata[:, genes].X.mean(axis=1)).ravel()
    return score >= np.quantile(score, 1 - top_fraction)


def epithelial_effects(
    normal: ad.AnnData, cancer: ad.AnnData, genes: list[str], composition: dict, pseudocount: float,
) -> pd.DataFrame:
    """Pseudo-bulk log2 fold change using only the epithelium-rich spots of each slide."""
    rates = {}
    for slide, adata in (("normal", normal), ("cancer", cancer)):
        sub = adata[epithelial_mask(adata, composition["epithelial_genes"],
                                    composition["top_fraction"])]
        counts = sub[:, genes].layers["counts"]
        counts = np.asarray(counts.sum(axis=0)).ravel()
        rates[slide] = counts / sub.obs["total_counts"].sum() * 1e4
    ratio = (rates["cancer"] + pseudocount) / (rates["normal"] + pseudocount)
    return pd.DataFrame({"gene": genes, "log2fc_epithelial": np.log2(ratio)})


def assign_tier(beyond_background: bool, concordant: bool, robust: bool, detected: bool) -> str:
    """A: all three evidence criteria. B: two of three. '-': fewer, or too sparsely detected."""
    if not detected:
        return "-"
    hits = int(beyond_background) + int(concordant) + int(robust)
    return {3: "A", 2: "B"}.get(hits, "-")
