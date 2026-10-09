"""Cross-check the new environmental and predicted-function Results sections."""

from __future__ import annotations

import hashlib
import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from manuscript_paths import PAPER
import manuscript_text as mt
from manuscript_text import fmt, pct


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
    assert (results["spearman_rho"] < 0).all()
    assert results["q_global_9"].max() == pytest.approx(0.02059311071)
    assert decision["genus_tests"] == 600
    by_climate = decision["genus_tests_q_lt_0_05_by_climate"]
    assert by_climate == {
        "mean_air_temperature_c": 112,
        "mean_monthly_rain_mm": 112,
        "mean_relative_humidity_pct": 111,
    }
    main = mt.section("main", "Climate covaries with the geographic pattern", "Soil pH tracks")
    assert "All correlations between the climate variables and the three alpha diversity measures were negative and significant" in main
    assert f"q≤{fmt(results['q_global_9'].max(), 3)}" in main
    counts = sorted(by_climate.values())
    assert f"({counts[0]}–{counts[-1]} genera per variable" in main
    assert "200 genera tested" in main
    # Climate summary ranges: one decimal in the main text, two in the supplement.
    site = pd.read_csv(ENVIRONMENT / "climate_site_summary.tsv", sep="\t")
    temperature, rain, humidity = (site[c] for c in ("mean_air_temperature_c", "mean_monthly_rain_mm", "mean_relative_humidity_pct"))
    assert re.search(
        rf"temperatures ranged from {fmt(temperature.min(), 1)} to {fmt(temperature.max(), 1)}\^\\circC, "
        rf"monthly rainfall from {fmt(rain.min(), 1)} to {fmt(rain.max(), 1)} mm, and relative humidity from "
        rf"{fmt(humidity.min(), 1)} to {fmt(humidity.max(), 1)}%",
        main,
    )
    supplement = mt.text("supplement")
    assert (
        f"{fmt(temperature.min(), 2)}-{fmt(temperature.max(), 2)}^\\circC for temperature, "
        f"{fmt(rain.min(), 2)}-{fmt(rain.max(), 2)} mm per month for rain and "
        f"{fmt(humidity.min(), 2)}-{fmt(humidity.max(), 2)}% for relative humidity"
    ) in supplement
    # Field-weather models: no association after correction.
    weather = pd.read_csv(ENVIRONMENT / "field_weather_adjusted_models.tsv", sep="\t")
    floor = f"{np.floor(weather['q_global_9'].min() * 1000) / 1000:.3f}"
    assert f"(all q≥{floor})" in main
    # Supplementary climate-diversity table.
    variables = {"Temperature": "mean_air_temperature_c", "Rain": "mean_monthly_rain_mm",
                 "Relative humidity": "mean_relative_humidity_pct"}
    measures = {"Shannon": "shannon", "Rarefied richness": "expected_richness_25k",
                "Normalized Shannon": "normalized_evenness"}
    current = None
    checked = 0
    for row in mt.table_rows("supplement", "tab:climate-alpha")[1:]:
        if row[0]:
            current = variables[row[0]]
        source = results[(results.climate_variable == current) & (results.response == measures[row[1]])].iloc[0]
        assert row[2] == (
            f"{fmt(source['spearman_rho'], 3)} [{fmt(source['bootstrap_ci_low'], 3)},{fmt(source['bootstrap_ci_high'], 3)}]"
        ), row
        q = source["q_global_9"]
        assert row[3] == (fmt(q, 4) if q >= 0.01 else _sci3(q)), row
        checked += 1
    assert checked == 9
    # Table 4: transect rho and genus counts per environmental variable.
    gradients = pd.read_csv(
        ROOT / "analysis/v3/taxon_context_corrected_20260909/site_gradients.tsv", sep="\t"
    ).set_index("variable")
    for label, variable in (("Temperature", "mean_air_temperature_c"), ("Relative humidity", "mean_relative_humidity_pct"),
                            ("Rainfall", "mean_monthly_rain_mm")):
        row = mt.table_row("main", "tab:environment", label)
        assert row[1:3] == [fmt(gradients.loc[variable, "spearman_rho_route_position"], 2), str(by_climate[variable])], row
    assert "kept all comparisons involving the same site together" in mt.text("main")
    _verify_checksums(ENVIRONMENT)


def test_environmental_limits_and_negative_diagnostics_are_in_the_right_places() -> None:
    main = mt.text("main")
    supplement = mt.text("supplement")
    # Climate, pH and elemental gradients are confounded with geography.
    assert "describe a shared climatic and geographic gradient rather than independent effects of climate" in main
    assert "our design cannot separate their individual effects" in main
    assert "Geographic position is therefore confounded with collection order and time of day" in main
    assert "remained spatially autocorrelated after accounting for the transect" in main
    assert re.search(r"transect trend (?:was supported only when|depended on the assumed spatial covariance)", main)
    # Field-weather null diagnostics stay in the supplement with their reason.
    weather = pd.read_csv(ENVIRONMENT / "field_weather_adjusted_models.tsv", sep="\t")
    assert (weather["q_global_9"] > 0.05).all()
    assert re.search(r"none[^.]*remained after correction", supplement)
    assert "capture transient conditions rather than the climate" in supplement
    # Measurement scope.
    assert "Climate summaries describe air, rather than soil, condition" in main
    assert "XRF measures total elemental composition rather than soluble ions or salinity" in main
    # Rainfall: exploratory, with the exchangeability assumption in the supplement.
    assert "Exploratory rainfall models identify the first days after rain as a window for event-based sampling" in main
    assert "exchangeability of the three observed annual rainfall fields" in supplement
    assert "q≥0.527" not in main and "q≥0.53" not in main


def test_picrust_ecology_matches_the_main_results_and_bounded_claim() -> None:
    decision = _json(PICRUST / "analysis_decision.json")
    assert decision["status"] == "predicted_functional_structure_supported"
    assert decision["cohort"]["ecology_profiles"] == 1227
    assert decision["cohort"]["grouped_profiles"] == 633
    assert decision["cohort"]["pathways"] == 462
    geography = decision["primary_geography"]
    assert geography["quadratic_transect_r2"] == pytest.approx(0.23326512170467784)
    assert geography["permutation_p"] == 0.0001
    assert geography["all_sensitivities_p_lt_0_05"] is True
    assert decision["pathway_level"] == {
        "geographic_tests": 200,
        "geographic_tests_q_lt_0_05": 92,
        "position_tests": 600,
        "position_tests_q_lt_0_05": 270,
    }
    sensitivity = decision["position_contrast_sensitivity"]
    assert sensitivity["Rhizosphere-Deep"]["all_p_lt_0_05"] is True
    assert sensitivity["Rhizosphere-Surface"]["all_p_lt_0_05"] is True
    assert sensitivity["Deep-Surface"]["all_p_lt_0_05"] is False

    main = mt.section("main", "Predicted functional potential follows geography and compartment", "We examined mechanisms")
    low, high = geography["quadratic_transect_r2_95_jackknife_ci"]
    tests = pd.read_csv(PICRUST / "geographic_profile_tests.tsv", sep="\t")
    assert (
        f"Transect position explained {pct(geography['quadratic_transect_r2'])}% of the variation in pathway composition "
        f"among sites (redundancy analysis; 95% CI {pct(low)}–{pct(high)}%; "
        f"{pct(tests.quadratic_transect_r2.min())}–{pct(tests.quadratic_transect_r2.max())}% across sensitivity analyses)"
    ) in main
    assert "92 of the 200 pathways were correlated with transect position (q<0.05)" in main
    position = pd.read_csv(PICRUST / "position_profile_tests.tsv", sep="\t")
    primary = position[position.analysis == "primary"].set_index("contrast")
    omnibus = primary.loc["omnibus_three_positions"]
    assert f"pseudo-F={fmt(omnibus['pseudo_f'], 2)}, p={omnibus['permutation_p']:.4f}; {int(omnibus['n_complete_blocks'])} site visits" in main
    root = [primary.loc[c] for c in ("Rhizosphere-Deep", "Rhizosphere-Surface")]
    assert (
        f"standardized displacement {fmt(root[0]['standardized_displacement'], 2)} and "
        f"{fmt(root[1]['standardized_displacement'], 2)}; q={fmt(max(r['q_primary_three'] for r in root), 4)}"
    ) in main
    deep = primary.loc["Deep-Surface"]
    assert f"({fmt(deep['standardized_displacement'], 2)}, q={fmt(deep['q_primary_three'], 3)}) was not significant when Trip 3 was omitted" in main
    omit3 = position[(position.analysis == "leave_one_campaign_out") & (position.omitted_campaign == 3)
                     & (position.contrast == "Deep-Surface")].iloc[0]
    assert omit3["permutation_p"] > 0.05
    assert "Of 600 pathway–compartment contrasts, 270 were significant (q<0.05)" in main
    # Prediction quality and genome-derived agreement.
    quality = decision["prediction_quality"]
    assert f"(NSTI) was {fmt(quality['weighted_nsti']['median'], 3)}" in main
    metrics = pd.read_csv(ROOT / "analysis/v3/measured_function_summary_results/summary_metrics.tsv", sep="\t").set_index("metric")
    median = metrics.loc["per_sample_ko_profile_spearman_median"]
    assert (
        f"median Spearman ρ={fmt(median['estimate'], 2)}, 95% CI {fmt(median['interval_low'], 2)}–{fmt(median['interval_high'], 2)}"
    ) in main
    assert f"(ρ={fmt(quality['community_mean_ko_spearman'], 2)})" in main
    assert f"Across {quality['shared_kos']:,} shared KOs" in main
    assert f"using the {quality['shotgun_matched_samples']} libraries with matched amplicon profiles" in main
    # Pathway-geography table: the printed rhos are the strongest supported ones.
    correlations = pd.read_csv(PICRUST / "pathway_geographic_correlations.tsv", sep="\t")
    supported = correlations[correlations.supported_q_lt_0_05]
    rows = mt.table_rows("supplement", "tab:pathway_geography")[1:]
    printed_up = [row[1] for row in rows]
    printed_down = [row[3] for row in rows]
    expected_up = [fmt(v, 2) for v in supported.spearman_rho.sort_values(ascending=False).head(len(rows))]
    expected_down = [fmt(v, 2) for v in supported.spearman_rho.sort_values().head(len(rows))]
    assert printed_up == expected_up and printed_down == expected_down
    assert (supported.q_global_200.loc[supported.spearman_rho.abs() >= abs(float(printed_down[-1]))] < 0.001).all()
    ranking = pd.read_csv(PICRUST / "pathway_ranking.tsv", sep="\t")
    assert f"Of the 462 predicted pathways, {len(ranking)} were present in at least 20% of the 633 grouped sample profiles" in mt.text("supplement")
    # Bounded interpretation: potential, not activity.
    assert "We did not measure function" in mt.text("main")
    assert "dedicated UV-repair pathway" not in mt.text("main")
    _verify_checksums(PICRUST)


def _sci3(value: float) -> str:
    exponent = int(f"{value:.2e}".split("e")[1])
    return f"{value / 10 ** exponent:.2f}×10^{exponent}"
