from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = (
    ROOT / "analysis/v3/make_submission_figures.py"
).read_text(encoding="utf-8")


def test_reviewed_figure_runtime_is_explicit_and_matches_recipe() -> None:
    runtime = json.loads(
        (
            ROOT
            / "empty-quarter-amplicon/figures/figure_runtime.json"
        ).read_text(encoding="utf-8")
    )
    dependencies = yaml.safe_load(
        (ROOT / "environment/environment.yml").read_text(encoding="utf-8")
    )["dependencies"]
    assert runtime == {
        "schema_version": "figure-runtime-v1",
        "python": "3.11.14",
        "matplotlib": "3.9.4",
        "freetype": "2.14.3",
    }
    for component in ("python", "matplotlib", "freetype"):
        assert f"{component}={runtime[component]}" in dependencies


def test_generator_fails_closed_before_rendering_with_another_runtime() -> None:
    assert "EXPECTED_FIGURE_RUNTIME" in GENERATOR
    assert "if observed != EXPECTED_FIGURE_RUNTIME" in GENERATOR
    assert "environment/conda-linux-64.lock" in GENERATOR
