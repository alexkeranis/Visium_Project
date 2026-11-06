"""Pipeline stages: prepare -> image -> markers -> report."""

from __future__ import annotations

import os
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd

from . import background, concordance, effects, io, normalize, panel, qc, repro, texture

NOTEBOOKS = ("01_data_qc.ipynb", "02_gene_markers.ipynb", "03_morphology_concordance.ipynb")


def _child_seeds(seed: int, n: int) -> list[int]:
    return [int(s.generate_state(1)[0]) for s in np.random.SeedSequence(seed).spawn(n)]


def prepare(config: dict) -> None:
    """Verify inputs, apply spot QC, normalize, write data/processed/<slide>.h5ad."""
    repro.verify_inputs(config)
    spots = []
    for slide in io.SLIDES:
        adata = qc.annotate(io.read_visium(config, slide), config["qc"])
        spots.append(qc.spot_table(adata))
        adata = adata[adata.obs["qc_pass"]].copy()
        normalize.log_normalize(adata, config["normalization"]["target_sum"])
        io.write_slide(adata, config, slide)
        print(f"{slide}: {adata.n_obs} spots kept")
    spots = pd.concat(spots, ignore_index=True)
    io.write_table(spots, config, "qc_spots.csv.gz", key="processed")
    io.write_table(qc.summary(spots), config, "qc_summary.csv")


def image(config: dict) -> None:
    """H&E stain-density and co-occurrence features for every kept spot."""
    tables = []
    for slide in io.SLIDES:
        adata = io.read_slide(config, slide)
        feats = texture.spot_features(adata, config["texture"])
        feats.insert(0, "barcode", feats.index)
        feats.insert(0, "slide", slide)
        tables.append(feats.reset_index(drop=True))
        print(f"{slide}: image features for {len(feats)} spots")
    io.write_table(pd.concat(tables, ignore_index=True), config, "image_features.csv.gz",
                   key="processed")


def markers(config: dict) -> None:
    """Effect sizes, slide background, literature agreement, concordance, candidate table."""
    seed_scores, seed_sets = _child_seeds(config["seed"], 2)
    norm_cfg, mk = config["normalization"], config["markers"]
    slides = io.read_slides(config)
    normal, cancer = slides["normal"], slides["cancer"]

    # Full-transcriptome effects and the matched-expression background
    genes_all = effects.gene_effects(normal, cancer, norm_cfg["pseudocount"])
    io.write_table(genes_all, config, "gene_effects_all.csv.gz", key="processed")
    pt = panel.panel_table(config)
    pt = pt[pt["gene"].isin(normal.var_names)].reset_index(drop=True)
    genes = pt["gene"].tolist()
    bg = background.background_percentiles(
        genes_all, genes, mk["min_detection"], mk["n_expression_bins"])

    # Housekeeping normalization as a sensitivity analysis
    hk_cfg = norm_cfg["housekeeping"]
    hk_check = panel.housekeeping_check(normal, cancer, hk_cfg)
    io.write_table(hk_check, config, "housekeeping_check.csv")
    hk = panel.housekeeping_effects(normal, cancer, genes, hk_cfg)
    epi = panel.epithelial_effects(normal, cancer, genes, config["composition"],
                                   norm_cfg["pseudocount"])

    # Spot-level features on both slides: panel genes, set scores, log-ratios
    combined = ad.concat(slides, index_unique="|")  # obs_names: "<barcode>|<slide>"
    sets = panel.directional_sets(pt, combined.var_names)
    score_names = panel.add_set_scores(combined, sets, seed_scores)
    ratio_names = panel.add_log_ratios(combined, config["ratios"])
    gene_values = pd.DataFrame(
        _dense(combined[:, genes].X), columns=genes, index=combined.obs_names)
    spot = pd.concat([gene_values, combined.obs[score_names + ratio_names]], axis=1)
    keyed = spot[score_names + ratio_names].copy()
    keyed.insert(0, "barcode", combined.obs_names.str.rsplit("|", n=1).str[0])
    keyed.insert(0, "slide", combined.obs["slide"].to_numpy())
    io.write_table(keyed, config, "spot_scores.csv.gz", key="processed")
    spot_fx = effects.feature_effects(spot, combined.obs["slide"])

    # Reference: how well random gene sets separate the two slides
    y = (combined.obs["slide"] == "cancer").to_numpy().astype(int)
    pool = np.flatnonzero(background.detected(genes_all, mk["min_detection"]).to_numpy())
    random_rows, score_rows = [], []
    for direction, members in sets.items():
        null = background.random_set_separability(
            combined.X, y, pool, len(members), mk["n_random_sets"], seed_sets)
        random_rows.append(pd.DataFrame({"set": direction, "set_size": len(members),
                                         "separability": null}))
        set_mean = np.asarray(combined[:, members].X.mean(axis=1)).ravel()
        observed = background.separability(y, set_mean)
        fx = spot_fx.set_index("feature").loc[f"score_{direction}"]
        score_rows.append({
            "feature": f"score_{direction}", "genes": " ".join(members),
            "cohens_d": fx["cohens_d"], "prob_superiority": fx["prob_superiority"],
            "set_mean_separability": observed,
            "random_set_pct": 100 * np.mean(null < observed),
        })
    for name in ratio_names:
        fx = spot_fx.set_index("feature").loc[name]
        score_rows.append({"feature": name, "cohens_d": fx["cohens_d"],
                           "prob_superiority": fx["prob_superiority"]})
    io.write_table(pd.concat(random_rows, ignore_index=True), config, "random_set_separability.csv")

    # Within-slide link to H&E image features
    img = io.read_table(config, "image_features.csv.gz", key="processed")
    img.index = img["barcode"] + "|" + img["slide"]
    img = img.drop(columns=["slide", "barcode"]).loc[combined.obs_names]
    both = genes_all.set_index("gene").loc[genes, ["detect_normal", "detect_cancer"]].min(axis=1)
    mol_cols = [g for g in genes if both[g] >= mk["min_detection"]] + score_names + ratio_names
    conc = concordance.within_slide(
        spot[mol_cols], img, combined.obs["slide"], config["concordance"]["min_abs_rho"])
    io.write_table(conc, config, "image_concordance.csv")
    link = concordance.strongest_link(conc)
    io.write_table(effects.feature_effects(img, combined.obs["slide"]), config,
                   "image_between_slides.csv")

    scores = pd.DataFrame(score_rows).merge(
        link.rename(columns={"image_feature": "image_link", "rho": "image_link_rho"})
        [["feature", "image_link", "image_link_rho"]], on="feature", how="left")
    io.write_table(scores, config, "score_effects.csv")

    # Candidate table
    cand = (
        pt.merge(genes_all, on="gene")
        .merge(bg, on="gene")
        .merge(hk, on="gene")
        .merge(epi, on="gene")
        .merge(spot_fx.rename(columns={"feature": "gene"})[["gene", "prob_superiority"]],
               on="gene")
        .merge(link.rename(columns={"feature": "gene", "image_feature": "image_link",
                                    "rho": "image_link_rho"})
               [["gene", "image_link", "image_link_rho"]], on="gene", how="left")
    )
    cand["observed"] = np.where(cand["log2fc"] > 0, "up", "down")
    cand["literature"] = [panel.literature_agreement(e, f)
                          for e, f in zip(cand["expected"], cand["log2fc"])]
    cand["detected"] = cand[["detect_normal", "detect_cancer"]].max(axis=1) >= mk["min_detection"]
    cand["beyond_background"] = cand["bg_pct_log2fc"] >= mk["background_percentile"]
    signs = np.sign(cand[["log2fc", "log2fc_hk", "log2fc_epithelial"]])
    cand["direction_robust"] = signs.nunique(axis=1) == 1
    cand["tier"] = [
        panel.assign_tier(b, lit == "concordant", s, d)
        for b, lit, s, d in zip(cand["beyond_background"], cand["literature"],
                                cand["direction_robust"], cand["detected"])
    ]
    cand["tier_rank"] = cand["tier"].map({"A": 0, "B": 1, "-": 2})
    cand = cand.sort_values(["tier_rank", "bg_pct_log2fc"], ascending=[True, False])
    cols = ["gene", "program", "tier", "observed", "expected", "literature", "pmid",
            "log2fc", "cohens_d", "prob_superiority", "detect_normal", "detect_cancer",
            "bg_pct_log2fc", "bg_pct_d", "n_matched", "log2fc_hk", "cohens_d_hk",
            "log2fc_epithelial", "beyond_background", "direction_robust", "detected",
            "image_link", "image_link_rho"]
    io.write_table(cand[cols], config, "candidate_markers.csv")
    print(cand[["gene", "tier", "observed", "expected", "log2fc", "bg_pct_log2fc"]]
          .head(12).to_string(index=False))


def report(config: dict) -> None:
    """Execute the notebooks in place; they read the tables above and write the figures."""
    nb_dir = io.project_root() / "notebooks"
    for name in NOTEBOOKS:
        print(f"executing {name}")
        subprocess.run(
            [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
             "--inplace", "--ExecutePreprocessor.timeout=1200",
             "--ExecutePreprocessor.kernel_name=python3", str(nb_dir / name)],
            check=True, env={**os.environ, "PYTHONHASHSEED": str(config["seed"])},
        )
        _strip_execution_timing(nb_dir / name)


def _strip_execution_timing(path) -> None:
    """Remove per-cell execution timestamps so reruns leave the notebook byte-identical."""
    import nbformat

    nb = nbformat.read(path, as_version=4)
    for cell in nb.cells:
        cell.metadata.pop("execution", None)
    nbformat.write(nb, path)


def _dense(X) -> np.ndarray:
    return np.asarray(X.todense() if hasattr(X, "todense") else X)


STAGES = {"prepare": prepare, "image": image, "markers": markers, "report": report}
