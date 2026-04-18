"""Reproducibility package for the Empty Quarter amplicon ecology paper."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
CACHE_DIR = REPO_ROOT / "cache"
FIGURES_DIR = REPO_ROOT / "figures"

CACHE_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

__all__ = ["REPO_ROOT", "DATA_DIR", "CACHE_DIR", "FIGURES_DIR"]
