"""Configuration, paths, raw Visium input and processed artifacts."""

from __future__ import annotations

import tarfile
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

SLIDES = ("normal", "cancer")


def project_root() -> Path:
    """Folder that contains config.yaml (works for editable installs and notebooks)."""
    for start in (Path.cwd(), Path(__file__).resolve()):
        for parent in (start, *start.parents):
            if (parent / "config.yaml").is_file():
                return parent
    raise FileNotFoundError("config.yaml not found above the working directory or the package.")


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else project_root() / "config.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def path(config: dict, key: str) -> Path:
    """Absolute path for config['paths'][key]; directories are created on demand."""
    p = project_root() / config["paths"][key]
    (p.parent if p.suffix else p).mkdir(parents=True, exist_ok=True)
    return p


def sample_dir(config: dict, slide: str) -> Path:
    return project_root() / config["paths"]["data"] / config["samples"][slide]


def ensure_spatial_folder(folder: Path) -> None:
    """Unpack <sample>_spatial.tar.gz into <sample>/spatial if that folder is missing."""
    if (folder / "spatial").is_dir():
        return
    archive = folder / f"{folder.name}_spatial.tar.gz"
    with tarfile.open(archive) as tar:
        tar.extractall(folder, filter="data")


def read_visium(config: dict, slide: str) -> ad.AnnData:
    """Read one Space Ranger output folder with its hires H&E image."""
    import squidpy as sq

    folder = sample_dir(config, slide)
    ensure_spatial_folder(folder)
    adata = sq.read.visium(
        path=folder,
        counts_file="filtered_feature_bc_matrix.h5",
        library_id=slide,
        load_images=True,
    )
    adata.var_names_make_unique()
    adata.obs["slide"] = slide
    return adata


def processed_file(config: dict, name: str) -> Path:
    return path(config, "processed") / name


def write_slide(adata: ad.AnnData, config: dict, slide: str) -> Path:
    out = processed_file(config, f"{slide}.h5ad")
    adata.write_h5ad(out, compression="gzip")
    return out


def read_slide(config: dict, slide: str) -> ad.AnnData:
    file = processed_file(config, f"{slide}.h5ad")
    if not file.exists():
        raise FileNotFoundError(f"{file} is missing. Run: python scripts/run_pipeline.py")
    return ad.read_h5ad(file)


def read_slides(config: dict) -> dict[str, ad.AnnData]:
    return {slide: read_slide(config, slide) for slide in SLIDES}


def write_table(df: pd.DataFrame, config: dict, name: str, key: str = "tables") -> Path:
    """Write a CSV with fixed float formatting so reruns are byte-identical."""
    out = path(config, key) / name
    df.to_csv(out, index=False, float_format="%.6g", lineterminator="\n")
    return out


def read_table(config: dict, name: str, key: str = "tables") -> pd.DataFrame:
    return pd.read_csv(path(config, key) / name)
