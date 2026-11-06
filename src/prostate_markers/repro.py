"""Seeding and the run manifest."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import random
import subprocess
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

from . import io

TRACKED_PACKAGES = (
    "anndata", "matplotlib", "numpy", "pandas", "scanpy",
    "scikit-image", "scipy", "seaborn", "squidpy",
)


def seed_everything(seed: int) -> None:
    """Seed Python and NumPy global state. Stochastic calls also get explicit seeds."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def sha256(file: Path) -> str:
    h = hashlib.sha256()
    with open(file, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_checksums(file: Path) -> dict[str, str]:
    """Parse a `sha256sum`-style file into {relative path: digest}."""
    out = {}
    for line in file.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            out[name.lstrip("*").strip()] = digest
    return out


def verify_inputs(config: dict) -> dict[str, str]:
    """Check every raw input against data/CHECKSUMS.sha256 and return the digests."""
    data = io.project_root() / config["paths"]["data"]
    expected = read_checksums(data / "CHECKSUMS.sha256")
    actual = {}
    for name, digest in expected.items():
        actual[name] = sha256(data / name)
        if actual[name] != digest:
            raise ValueError(f"Checksum mismatch for data/{name}")
    return actual


def _git(*args: str) -> str | None:
    try:
        return subprocess.run(
            ["git", *args], cwd=io.project_root(), capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _version(pkg: str) -> str | None:
    try:
        return version(pkg)
    except PackageNotFoundError:
        return None


def write_manifest(config: dict, stages: list[str]) -> Path:
    config_file = io.project_root() / "config.yaml"
    manifest = {
        "finished_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "stages": stages,
        "seed": config["seed"],
        "config_sha256": sha256(config_file),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain", "--untracked-files=no")),
        "python": platform.python_version(),
        "packages": {pkg: _version(pkg) for pkg in TRACKED_PACKAGES},
        "inputs_sha256": verify_inputs(config),
    }
    out = io.path(config, "manifest")
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    return out
