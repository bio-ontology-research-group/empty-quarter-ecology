import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.v3.claim_rescue import (
    TYPE_LABELS,
    TYPE_ORDER,
    WINDOWS,
    bh_fdr,
    downstream_decisions,
)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis/v3/results"
RAIN = RESULTS / "rain_window_models.tsv"
LOCO = RESULTS / "rain_leave_one_campaign.tsv"
LEDGER = RESULTS / "claim_ledger.tsv"

# Corrected site-campaign denominators after the Trip-3 "19.5" parser fix.
EXPECTED_DENOMINATORS = {"Surface": 198, "Deep": 211, "Rhizosphere": 224}


def test_bh_fdr_is_monotone_in_rank_and_bounded():
    raw = np.array([0.04, 0.001, 0.02, 0.8])
    adjusted = bh_fdr(raw)
    assert np.all((0 <= adjusted) & (adjusted <= 1))
    order = np.argsort(raw)
    assert np.all(np.diff(adjusted[order]) >= 0)


def test_bh_fdr_known_values():
    adjusted = bh_fdr([0.01, 0.04, 0.03])
    assert np.allclose(adjusted, [0.03, 0.04, 0.04])


@pytest.fixture(scope="module")
def rain() -> pd.DataFrame:
    if not RAIN.exists():
        pytest.skip(f"{RAIN} not staged")
    return pd.read_csv(RAIN, sep="\t")


def test_rain_family_covers_every_compartment_and_window(rain):
    assert len(rain) == len(TYPE_ORDER) * len(WINDOWS) == 18
    observed = {
        (row.compartment, row.window_days) for row in rain.itertuples()
    }
    assert observed == {
        (compartment, window)
        for compartment in TYPE_ORDER
        for window in WINDOWS
    }
    assert set(rain["compartment_label"]) == set(TYPE_LABELS.values())


def test_rain_denominators_match_the_corrected_geodata_parse(rain):
    counts = (
        rain.groupby("compartment")["n_site_campaigns"].agg(["nunique", "max"])
    )
    for compartment, expected in EXPECTED_DENOMINATORS.items():
        # One denominator per compartment across all six windows.
        assert counts.loc[compartment, "nunique"] == 1
        assert counts.loc[compartment, "max"] == expected
    assert set(rain["n_sites"]) == {60}


def test_rain_family_declares_one_global_correction(rain):
    assert set(rain["n_tests_in_global_family"]) == {18}
    assert (rain["q_global"] <= 1).all()
    assert (rain["q_within_compartment"] <= 1).all()
    # q_global is one Benjamini-Hochberg family over all 18 tests.
    assert np.allclose(rain["q_global"], bh_fdr(rain["p"]))
    # q_within_compartment is a separate six-test family per compartment.
    for compartment, part in rain.groupby("compartment"):
        assert len(part) == len(WINDOWS)
        assert np.allclose(part["q_within_compartment"], bh_fdr(part["p"]))


def test_rain_family_extreme_values_match_the_reported_pair(rain):
    # The manuscript reports the smallest within-compartment and the
    # smallest global q for this family.
    assert round(float(rain["q_within_compartment"].min()), 3) == 0.066
    assert round(float(rain["q_global"].min()), 3) == 0.185
    strongest = rain.loc[rain["p"].idxmin()]
    assert strongest["compartment"] == "Deep"
    assert strongest["window_days"] == 3
    assert round(float(strongest["estimate_per_mm"]), 3) == 0.208
    assert round(float(strongest["p"]), 4) == 0.0110


def test_surface_fourteen_day_model_matches_the_reported_values(rain):
    row = rain[
        (rain["compartment"] == "Surface") & (rain["window_days"] == 14)
    ].iloc[0]
    assert int(row["n_site_campaigns"]) == 198
    assert round(float(row["estimate_per_mm"]), 5) == -0.02268
    assert round(float(row["ci_low"]), 5) == -0.05607
    assert round(float(row["ci_high"]), 5) == 0.01071
    assert round(float(row["q_within_compartment"]), 5) == 0.21961


def test_no_rain_association_survives_the_declared_family(rain):
    assert (rain["q_global"] >= 0.05).all()
    assert (rain["q_within_compartment"] >= 0.05).all()


def test_rain_leave_one_campaign_covers_the_whole_family():
    if not LOCO.exists():
        pytest.skip(f"{LOCO} not staged")
    loco = pd.read_csv(LOCO, sep="\t")
    # Five campaigns per compartment-by-window cell.
    assert len(loco) == 18 * 5
    grouped = loco.groupby(["compartment", "window_days"])["omitted_trip"]
    assert set(grouped.nunique()) == {5}


def test_claim_ledger_has_no_pending_advanced_rows():
    if not LEDGER.exists():
        pytest.skip(f"{LEDGER} not staged")
    ledger = pd.read_csv(LEDGER, sep="\t")
    assert "pending_advanced_rescue" not in set(ledger["status"])
    assert "antecedent-rainfall response in any compartment" in set(
        ledger["claim"]
    )


def test_legacy_rain_artifacts_are_explicitly_superseded(rain):
    assert set(rain["artifact_status"]) == {"SUPERSEDED"}
    assert set(rain["superseded_by"]) == {
        "analysis/v3/rain_response_window/"
    }

    loco = pd.read_csv(LOCO, sep="\t")
    assert set(loco["artifact_status"]) == {"SUPERSEDED"}
    assert set(loco["superseded_by"]) == {
        "analysis/v3/rain_response_window/"
    }

    ledger = pd.read_csv(LEDGER, sep="\t")
    rain_rows = ledger[ledger["claim"].str.contains("rain", case=False)]
    assert len(rain_rows) == 2
    assert set(rain_rows["status"]) == {"superseded"}
    assert set(rain_rows["evidence_file"]) == {
        "../rain_response_window/analysis_decision.json"
    }

    readme = (RESULTS / "README.md").read_text(encoding="utf-8")
    sidecar = (RESULTS / "SUPERSEDED_RAIN_ANALYSIS.md").read_text(
        encoding="utf-8"
    )
    assert "Supersession warning" in readme
    assert "Status: **SUPERSEDED**" in sidecar


def test_result_checksums_cover_every_non_checksum_file():
    checksum_path = RESULTS / "SHA256SUMS"
    assert checksum_path.exists(), f"{checksum_path} must be staged"
    listed = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        listed[name] = digest

    expected = {
        path.name for path in RESULTS.iterdir()
        if path.is_file() and path.name != "SHA256SUMS"
    }
    assert set(listed) == expected
    for name, digest in listed.items():
        observed = hashlib.sha256((RESULTS / name).read_bytes()).hexdigest()
        assert observed == digest


def test_downstream_decisions_import_the_canonical_verdicts():
    decisions = downstream_decisions(ROOT)
    by_claim = {item.claim: item for item in decisions}
    assert set(by_claim) == {
        "network restructuring",
        "dispersal limitation from betaNTI/RCBray",
        "functional redundancy",
        "trace-gas-powered carbon fixation",
    }
    verdict = json.loads(
        (
            ROOT
            / "analysis/v3/network_rescue/results/claim_verdict.json"
        ).read_text()
    )
    assert (
        by_claim["network restructuring"].permitted_wording
        == verdict["permitted_wording"]
    )
    assert all(item.status != "pending_advanced_rescue" for item in decisions)


def test_downstream_decisions_fail_closed_on_a_missing_verdict(tmp_path):
    with pytest.raises(FileNotFoundError):
        downstream_decisions(tmp_path)
