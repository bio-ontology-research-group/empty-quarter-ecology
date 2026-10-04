"""The rainfall analyses consume the Open-Meteo input with the corrected Site 52 series."""
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "analysis/v3/open_meteo_site52_corrected_20261004"
CORRECTED = PACKAGE / "daily_weather_canonical_site52_corrected.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_package_records_corrected_request_and_output_hash():
    manifest = json.loads((PACKAGE / "manifest.json").read_text())
    assert manifest["data_repository_commit"].startswith("5a17782")
    assert manifest["source"]["sha256"].startswith("6cb570e7")
    request = manifest["requests"]["corrected"]
    assert (request["latitude"], request["longitude"]) == (20.82784, 53.57835)
    assert request["start_date"] == "2022-01-01" and request["end_date"] == "2026-02-01"
    response = PACKAGE / "open_meteo_site52_corrected_response.json"
    assert sha256(response) == request["response_sha256"]
    check = manifest["verification"]
    assert check["averaged_position_reproduces_pinned_precipitation"] is True
    assert check["rows_replaced"] == 1493
    assert manifest["output"]["sha256"] == sha256(CORRECTED)


def test_site52_series_is_the_corrected_position_response():
    table = pd.read_csv(CORRECTED, sep="\t", dtype={"Site": str})
    site = table[table.Site == "52"].set_index("Date")["Precip_mm"]
    daily = json.loads((PACKAGE / "open_meteo_site52_corrected_response.json").read_text())["daily"]
    response = pd.Series(daily["precipitation_sum"], index=daily["time"], dtype=float)
    assert len(site) == len(response) == 1493
    assert (site.reindex(response.index).to_numpy() == response.to_numpy()).all()
    assert round(site.sum(), 1) == 183.0


def test_rain_analyses_consume_the_corrected_input():
    refit = json.loads(
        (ROOT / "analysis/v3/rain_calendar_refit_20260909/manifest.json").read_text()
    )
    inputs = {Path(item["path"]).name: item["sha256"] for item in refit["inputs"]}
    assert inputs[CORRECTED.name] == sha256(CORRECTED)
    assert "daily_weather_canonical.tsv" not in inputs
    for script in ("run_rain_pulse_suite.py", "rain_response_window.py"):
        text = (ROOT / "analysis/v3" / script).read_text()
        assert "open_meteo_site52_corrected_20261004" in text
        assert '"data/processed/climate/daily_weather_canonical.tsv"' not in text


PRODUCTS = ROOT / "analysis/v3/rain_event_product_exposures_site52_corrected_20261004"


def test_five_product_exposures_use_corrected_site52_position_and_inputs():
    table = pd.read_csv(PRODUCTS / "rain_event_product_exposures.tsv", sep="\t")
    site = table[table.site == 52]
    assert len(site) == 5
    assert set(site.latitude) == {20.82784} and set(site.longitude) == {53.57835}
    grids = dict(zip(site.product_id, site.grid_id))
    assert grids["CHIRPS_v2.0"] == "CHIRPS_20.82500_53.57500"
    assert grids["CMORPH_V1.0_ADJ"] == "CMORPH_20.875_53.625"
    assert grids["GPM_3IMERGDF_07B"] == "IMERG_20.85_53.55"
    sources = pd.read_csv(PRODUCTS / "rain_event_product_sources.tsv", sep="\t")
    open_meteo = sources[sources.product_id == "Open_Meteo_precipitation_sum"].iloc[0]
    assert open_meteo["sha256"] == sha256(CORRECTED)


def test_five_product_agreement_matches_reported_ranges():
    table = pd.read_csv(PRODUCTS / "product_correlations.tsv", sep="\t")
    assert len(table) == 10
    assert round(table.pearson_r.min(), 3) == 0.501 and round(table.pearson_r.max(), 3) == 0.917
    assert round(table.spearman_rho.min(), 3) == 0.554 and round(table.spearman_rho.max(), 3) == 0.881
