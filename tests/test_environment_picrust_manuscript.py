"""Cross-check the new environmental and predicted-function Results sections."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from manuscript_paths import PAPER


ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENT = ROOT / "analysis/v3/environment_associations"
PICRUST = ROOT / "analysis/v3/picrust2_ecology"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def _without_value_math(text: str) -> str:
    """Normalize numeric LaTeX typesetting for prose-level claim checks."""
    return text.replace("$", "").replace("{,}", ",")


def _verify_checksums(directory: Path) -> None:
    for line in (directory / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest


def test_environment_associations_match_the_main_results() -> None:
    decision = _json(ENVIRONMENT / "analysis_decision.json")
    results = pd.read_csv(
        ENVIRONMENT / "climate_alpha_correlations.tsv", sep="\t"
    )
    main = _without_value_math(_flat(PAPER / "main.tex"))

    assert decision["status"] == "observational_climate_associations_supported"
    assert decision["analysis_unit"] == "core site"
    assert decision["climate_coverage"] == {
        "core_sites": 60,
        "first_month": "2022-01",
        "last_month": "2026-01",
        "monthly_records": 2940,
        "months_per_site": 49,
    }
    assert len(results) == 9
    assert results["supported_q_lt_0_05"].all()
    assert results["q_global_9"].max() == pytest.approx(0.02059311071)
    assert decision["genus_tests"] == 600
    assert decision["genus_tests_q_lt_0_05_by_climate"] == {
        "mean_air_temperature_c": 112,
        "mean_monthly_rain_mm": 112,
        "mean_relative_humidity_pct": 111,
    }
    for text in (
        "Climate gradients track bacterial diversity and composition",
        "corrected the 9 climate--diversity tests together",
        "Among 200 genera detected in at least 20\\,\\% of profiles, 112, 112 and 111",
        "kept all comparisons involving the same site together",
    ):
        assert text in main
    _verify_checksums(ENVIRONMENT)


def test_environmental_limits_and_negative_diagnostics_are_in_the_right_places() -> None:
    main = _flat(PAPER / "main.tex")
    supplement = _flat(PAPER / "supplement.tex")

    assert "covaried with route position" in supplement
    assert "shared environmental and geographic gradient" in main
    assert "statistical support depended on the assumed spatial covariance" in main
    assert r"adjusted $q\geq0.527$" not in main
    assert r"adjusted $q\geq0.527$" in supplement
    assert "Rainfall models identify an early sampling window" in main
    assert "Fitted short-term rainfall pulse" in supplement
    assert "conditional on exchangeability of the annual rainfall fields" in main


def test_picrust_ecology_matches_the_main_results_and_bounded_claim() -> None:
    decision = _json(PICRUST / "analysis_decision.json")
    main = _without_value_math(_flat(PAPER / "main.tex"))
    supplement = _without_value_math(_flat(PAPER / "supplement.tex"))

    assert decision["status"] == "predicted_functional_structure_supported"
    assert decision["cohort"]["ecology_profiles"] == 1227
    assert decision["cohort"]["grouped_profiles"] == 633
    assert decision["cohort"]["pathways"] == 462
    assert decision["primary_geography"]["quadratic_transect_r2"] == pytest.approx(
        0.23326512170467784
    )
    assert decision["primary_geography"]["permutation_p"] == 0.0001
    assert decision["primary_geography"]["all_sensitivities_p_lt_0_05"] is True
    assert decision["pathway_level"] == {
        "geographic_tests": 200,
        "geographic_tests_q_lt_0_05": 92,
        "position_tests": 600,
        "position_tests_q_lt_0_05": 270,
    }
    assert decision["position_contrast_sensitivity"]["Rhizosphere-Deep"][
        "all_p_lt_0_05"
    ] is True
    assert decision["position_contrast_sensitivity"]["Rhizosphere-Surface"][
        "all_p_lt_0_05"
    ] is True
    assert decision["position_contrast_sensitivity"]["Deep-Surface"][
        "all_p_lt_0_05"
    ] is False

    for text in (
        "Predicted metabolic pathways follow geography and compartment",
        "23.3\\,\\% of the differences among site-level pathway profiles",
        "Of 600 tests (200 pathways \\times 3 compartment pairs), 270",
        "geographic changes in fatty-acid pathways",
        "glucose, xylose and fructuronate degradation",
        "substrate-uptake and metabolite-exchange experiments",
        "median rank correlation was 0.743",
        "KO annotations of metagenome-assembled genomes (MAGs) with their CoverM abundance estimates",
    ):
        assert text in main
    assert "did not pass when Trip~3 was omitted" not in main
    assert "did not pass when Trip~3 was omitted" in supplement
    assert "dedicated UV-repair pathway" not in main
    assert "UV protection" in supplement
    _verify_checksums(PICRUST)
