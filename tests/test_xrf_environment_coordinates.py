"""Coordinate provenance of the XRF and climate-association results.

The Trip 1 and Trip 3 field sheets carried the Site 53 position for Site 52.
These checks pin the regenerated XRF and environment-association outputs to
the corrected position and to the locked data revision.
"""

from __future__ import annotations

import hashlib
import json
import subprocess

from manuscript_paths import DATA_REPO
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
XRF = ROOT / "analysis/v3/xrf_community_rescue"
ENVIRONMENT = ROOT / "analysis/v3/environment_associations"
SITE_52 = (20.82784, 53.57835)
OLD_SITE_52 = (20.851514166666668, 53.75788444444444)


def _lock() -> dict[str, str]:
    frame = pd.read_csv(ROOT / "DATA_REPOSITORY.lock", sep="\t", dtype=str)
    return dict(zip(frame["field"], frame["value"]))


def _verify_data_revision(recorded_commit: str) -> None:
    # Keep the run's original provenance; a later paper/evidence release may
    # advance the lock only while all original scientific input bytes agree.
    current = _lock()["commit"]
    subprocess.run(["git", "-C", str(DATA_REPO), "merge-base", "--is-ancestor",
                    recorded_commit, current], check=True)
    def scientific_blobs(commit):
        rows = subprocess.check_output(
            ["git", "-C", str(DATA_REPO), "ls-tree", "-r", commit,
             "metadata", "processed"], text=True).splitlines()
        return {row.split("\t", 1)[1]: row.split("\t", 1)[0] for row in rows}
    previous, updated = scientific_blobs(recorded_commit), scientific_blobs(current)
    for relative, blob in previous.items():
        assert updated.get(relative) == blob, relative



def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_checksums(directory: Path) -> None:
    for line in (directory / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        assert _sha256(directory / name) == digest, name


def test_xrf_alpha_table_uses_corrected_site_52_coordinates() -> None:
    table = pd.read_csv(XRF / "xrf_alpha_analysis_table.tsv", sep="\t")
    site = table[table["Site"] == 52]
    assert set(site["Trip"]) == {1, 3, 4}
    assert site["Latitude"].to_numpy() == pytest.approx(SITE_52[0], abs=1e-6)
    assert site["Longitude"].to_numpy() == pytest.approx(SITE_52[1], abs=1e-6)
    assert not (
        (table["Site"] == 52)
        & (table["Latitude"].sub(OLD_SITE_52[0]).abs() < 1e-6)
    ).any()


def test_xrf_results_and_provenance() -> None:
    summary = json.loads((XRF / "xrf_claim_summary.json").read_text())
    alpha = pd.read_csv(XRF / "xrf_alpha_models.tsv", sep="\t").set_index("model")
    community = pd.read_csv(XRF / "xrf_community_models.tsv", sep="\t")
    manifest = json.loads((XRF / "run_manifest.json").read_text())

    assert summary["status"] == "not_supported"
    assert summary["counts"] == {
        "alpha_joined_observations": 622,
        "community_joined_observations": 621,
        "core_sites": 60,
        "laboratory_xrf_records": 725,
        "trip5_records": 178,
        "trips1_4_records": 547,
    }
    assert summary["axis"]["pc1_variance_explained"] == pytest.approx(0.392, abs=5e-4)
    assert summary["axis"]["pc1_evaporite_reference_correlation"] == pytest.approx(
        0.825, abs=5e-4
    )
    assert summary["axis"]["within_vs_global_pc1_correlation"] == pytest.approx(
        0.976, abs=5e-4
    )
    # Coordinate-free primary model; the geographic-trend sensitivity is the
    # only alpha model that reads the site coordinates.
    assert alpha.loc["site_fixed_primary", "p"] == pytest.approx(0.9896, abs=5e-5)
    assert alpha.loc["geographic_trend_sensitivity", "p"] > 0.05
    by_fraction = community.set_index("pcoa_target_fraction")
    assert by_fraction.loc[0.80, "pcoa_axes"] == 38
    assert by_fraction.loc[0.80, "partial_r2"] == pytest.approx(0.003455, abs=5e-7)
    assert by_fraction.loc[0.80, "p"] == pytest.approx(0.056)
    assert by_fraction.loc[0.95, "partial_r2"] == pytest.approx(0.003304, abs=5e-7)
    assert by_fraction.loc[0.95, "p"] == pytest.approx(0.044)

    _verify_data_revision(manifest["data_repository"]["commit"])
    for relative, entry in manifest["inputs"].items():
        assert _sha256(ROOT / relative) == entry["sha256"], relative
    for trip in range(1, 6):
        relative = f"data/metadata/geodata/trip{trip}_geodata.tsv"
        assert manifest["inputs"][relative]["sha256"] == _sha256(ROOT / relative)
    _verify_checksums(XRF)


def test_environment_site_summary_uses_corrected_site_52_coordinates() -> None:
    site = pd.read_csv(ENVIRONMENT / "climate_site_summary.tsv", sep="\t")
    row = site.set_index("site").loc[52]
    assert row["latitude"] == pytest.approx(SITE_52[0], abs=1e-6)
    assert row["longitude"] == pytest.approx(SITE_52[1], abs=1e-6)
    # The corrected position lies between Sites 51 and 53 on the route.
    order = site.sort_values("transect_km")["site"].tolist()
    assert order.index(51) < order.index(52) < order.index(53)


def test_environment_results_after_coordinate_correction() -> None:
    decision = json.loads((ENVIRONMENT / "analysis_decision.json").read_text())
    covariation = pd.read_csv(ENVIRONMENT / "climate_covariation.tsv", sep="\t")
    adjusted = pd.read_csv(ENVIRONMENT / "field_weather_adjusted_models.tsv", sep="\t")
    raw = pd.read_csv(ENVIRONMENT / "field_weather_raw_correlations.tsv", sep="\t")

    transect = covariation[covariation["variable_b"] == "transect_km"].set_index(
        "variable_a"
    )["spearman_rho"]
    assert transect["mean_air_temperature_c"] == pytest.approx(0.98208, abs=5e-6)
    assert transect["mean_relative_humidity_pct"] == pytest.approx(0.91889, abs=5e-6)
    assert transect["mean_monthly_rain_mm"] == pytest.approx(0.55093, abs=5e-6)

    # Collection weather: raw pressure and humidity correlations survive BH,
    # none survive campaign and quadratic transect adjustment.
    assert (raw.loc[raw["weather_variable"] == "pressure_mbar", "q_global_9"] < 0.05).all()
    assert adjusted["q_global_9"].min() == pytest.approx(0.3192, abs=5e-5)
    assert adjusted["q_global_9"].max() == pytest.approx(0.7966, abs=5e-5)

    _verify_data_revision(decision["data_repository"]["commit"])
    for trip in range(1, 6):
        path = ROOT / f"data/metadata/geodata/trip{trip}_geodata.tsv"
        assert decision["coordinate_input_sha256"][f"trip{trip}_geodata"] == _sha256(path)
    # Field-weather sheet at data commit 5a17782 (Site 52 coordinates restored
    # in a3f1c21); older data checkouts carry 1b8cb5c7...
    assert decision["input_sha256"]["field_weather"] == (
        "51177fb6b21e4b10712daaa731cdf3783d0fc351c536f562851efd233386573f"
    )
    _verify_checksums(ENVIRONMENT)


def test_hypothesis_registry_uses_corrected_field_weather_and_calendar_refit() -> None:
    registry = json.loads(
        (ROOT / "analysis/v3/hypothesis_families_20260909/hypothesis_families.json").read_text()
    )
    members = pd.read_csv(
        ROOT / "analysis/v3/hypothesis_families_20260909/hypothesis_members.tsv", sep="\t"
    )
    family = {item["family"]: item for item in registry["families"]}["field_weather"]
    assert family["source_sha256"] == _sha256(ENVIRONMENT / "field_weather_adjusted_models.tsv")
    adjusted = pd.read_csv(ENVIRONMENT / "field_weather_adjusted_models.tsv", sep="\t")
    rows = members[members["family"] == "field_weather"]
    assert rows["q"].to_numpy() == pytest.approx(adjusted["q_global_9"].to_numpy(), rel=1e-12)
    assert rows["q"].min() == pytest.approx(0.3192, abs=5e-5)
    rainfall = registry["other_tests"]
    assert rainfall["rainfall_search"].startswith("Primary: calendar-year refit in rain_calendar_refit_20260909")
    assert "6 permutations" in rainfall["rainfall_search"]
    assert "circular" not in rainfall["rainfall_search"]
    assert rainfall["rainfall_search_historical"].startswith("Superseded diagnostic only")


@pytest.mark.parametrize(
    "directory",
    ["biology_context_corrected_20260909", "taxon_context_corrected_20260909"],
)
def test_corrected_context_runs_record_climate_input_reconciliation(directory: str) -> None:
    folder = ROOT / "analysis/v3" / directory
    manifest = json.loads((folder / "run_manifest.json").read_text())
    record = manifest["input_reconciliation"]["climate_site_summary"]
    # The digest of the table actually read by the run is retained.
    assert manifest["inputs"]["climate_site_summary"]["sha256"] == record["sha256_read_by_run"]
    assert record["sha256_read_by_run"] == (
        "f1e5c44b47d603d8f9f199500db95ed033ef37e17822049d28b24eebad0327c3"
    )
    assert record["sha256_current"] == _sha256(ENVIRONMENT / "climate_site_summary.tsv")
    assert record["changed_columns"] == ["latitude", "longitude", "transect_km"]
    assert not set(record["consumed_columns"]) & set(record["changed_columns"])
    assert record["consumed_columns_identical"] is True
    _verify_checksums(folder)
