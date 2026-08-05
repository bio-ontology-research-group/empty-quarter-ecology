"""Regression tests for the fitted antecedent-rain pulse analysis."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis/v3/rain_pulse_response.py"
PRIMARY = ROOT / "analysis/v3/rain_pulse_response"
OPEN_METEO = ROOT / "analysis/v3/rain_pulse_response_open_meteo"
SENSITIVITIES = ROOT / "analysis/v3/rain_pulse_sensitivities"


def load_module():
    specification = importlib.util.spec_from_file_location(
        "rain_pulse_response_for_tests", SCRIPT
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_pulse_kernel_rises_to_peak_and_then_decays():
    module = load_module()
    peak = 3.0
    days = np.asarray([0.5, 1.0, 2.0, peak, 4.0, 8.0, 20.0])
    values = module.pulse_kernel(days, peak)
    assert values[0] < values[1] < values[2] < values[3]
    assert values[3] == pytest.approx(1.0)
    assert values[3] > values[4] > values[5] > values[6]


def test_analysis_cohorts_use_all_available_groups_and_complete_pairs():
    module = load_module()
    base = module.load_base_module(ROOT)
    geodata = base.load_geodata(ROOT)
    alpha = ROOT / "analysis/v2/review/cache/alpha.tsv"
    grouped = module.load_group_means(alpha, geodata)
    paired = module.load_paired_means(alpha, geodata)
    assert len(grouped) == 617
    assert grouped["Type"].value_counts().to_dict() == {
        "Rhizosphere": 218,
        "Deep": 203,
        "Surface": 196,
    }
    assert len(paired) == 179
    assert paired["Site"].nunique() == 60


def test_primary_result_is_borderline_and_not_promoted_by_filtering():
    primary = read_json(PRIMARY / "analysis_decision.json")
    filtered = read_json(
        SENSITIVITIES / "control_filtered" / "analysis_decision.json"
    )
    assert primary["community_table_role"] == "primary_unfiltered"
    assert primary["weather_product"] == "nasa_power"
    assert primary["analysis_status"] == "temporally_localized_association_borderline"
    assert primary["selected_endpoint"] == "richness_hurlbert_25000"
    assert primary["selected_peak_complete_days"] == 2.0
    assert primary["selected_estimate_per_mm_at_peak"] == pytest.approx(
        179.5673223157866
    )
    family = primary["familywise_inference"]
    assert family["n_lag_rotation_draws"] == 19_999
    assert family["conditional_lag_rotation_one_sided_p"] == pytest.approx(0.056)
    assert family["conditional_lag_rotation_two_sided_p"] == pytest.approx(
        0.07865
    )

    assert filtered["community_table_role"] == "control_filtered_sensitivity"
    assert filtered["selected_peak_complete_days"] == 2.0
    assert abs(
        filtered["selected_estimate_per_mm_at_peak"]
        - primary["selected_estimate_per_mm_at_peak"]
    ) < 0.2
    assert filtered["familywise_inference"][
        "conditional_lag_rotation_one_sided_p"
    ] == pytest.approx(0.0485)
    assert "threshold crossing" in " ".join(primary["prohibited_wording"])


def test_independent_weather_product_places_peak_in_same_early_period():
    decision = read_json(OPEN_METEO / "analysis_decision.json")
    assert decision["community_table_role"] == "primary_unfiltered"
    assert decision["weather_product"] == "open_meteo"
    assert decision["selected_endpoint"] == "richness_hurlbert_25000"
    assert decision["selected_peak_complete_days"] == 1.0
    assert decision["familywise_inference"][
        "conditional_lag_rotation_one_sided_p"
    ] == pytest.approx(0.04315)
    assert decision["familywise_inference"][
        "conditional_lag_rotation_two_sided_p"
    ] == pytest.approx(0.04865)


def test_shape_position_and_future_rain_diagnostics_are_bounded():
    positions = pd.read_csv(PRIMARY / "soil_position_sensitivity.tsv", sep="\t")
    richness = positions[
        (positions["endpoint"] == "richness_hurlbert_25000")
        & (positions["cluster_definition"] == "site")
    ]
    assert set(richness["soil_position"]) == {
        "Surface",
        "Deep",
        "Rhizosphere",
    }
    assert (richness["estimate_per_mm_at_kernel_peak"] > 0).all()

    bins = pd.read_csv(PRIMARY / "disjoint_lag_bin_sensitivity.tsv", sep="\t")
    early = bins[
        (bins["lag_start_complete_days"] == 3)
        & (bins["lag_end_complete_days"] == 4)
    ].iloc[0]
    assert early["estimate_per_mm"] > 0
    assert bins.iloc[-1]["estimate_per_mm"] < 0

    decision = read_json(PRIMARY / "analysis_decision.json")
    placebo = decision["future_rain_placebo"]
    assert placebo["richness_maximum_t"] < 0.1
    assert placebo["inferential_role"] == "descriptive negative control only"


def test_suite_does_not_require_rainfall_stability_across_campaigns():
    manifest = read_json(SENSITIVITIES / "suite_manifest.json")
    summary = pd.read_csv(SENSITIVITIES / "summary.tsv", sep="\t")
    assert manifest["campaign_stability_requirement"] is False
    assert "rare and uneven" in manifest["reason"]
    assert manifest["consolidated_interpretation"]["reporting_verdict"] == (
        "bounded_observational_association"
    )
    assert len(summary) == 8
    assert (summary["estimate_per_mm_at_peak"] > 0).all()
    assert summary["selected_peak_complete_days"].max() <= 5.0


def test_all_run_checksums_are_current():
    outputs = [PRIMARY, OPEN_METEO]
    manifest = read_json(SENSITIVITIES / "suite_manifest.json")
    outputs.extend(ROOT / run["output"] for run in manifest["runs"])
    for output in dict.fromkeys(outputs):
        for line in (output / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, filename = line.split("  ", 1)
            assert sha256_file(output / filename) == digest


def test_fixed_event_exposures_cover_five_products_with_positive_agreement():
    table = pd.read_csv(
        ROOT / "data/processed/climate/rain_event_product_exposures.tsv", sep="\t"
    )
    assert len(table) == 300
    assert table["product_id"].nunique() == 5
    assert set(table.groupby("product_id").size()) == {60}
    wide = table.pivot(index="site", columns="product_id", values="precipitation_mm")
    correlations = wide.corr(method="spearman")
    off_diagonal = correlations.to_numpy()[
        ~np.eye(len(correlations), dtype=bool)
    ]
    assert (off_diagonal > 0).all()
