import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis" / "v3"))

from distance_decay_turnover import (  # noqa: E402
    bray_partition,
    euclidean_matrix,
    permutation_slopes,
    sorensen_partition,
    two_sided_p,
    upper_triangle,
)
from spatial_resolution_sensitivity import (  # noqa: E402
    DEFAULT_COORDINATES,
    NEIGHBOUR_COUNTS,
    asv_counts_on_reference_groups,
    load_site_coordinates,
    reference_sample_ids,
)

ASV_DIR = ROOT / "analysis/v3/spatial_resolution_sensitivity"
DECAY_DIR = ROOT / "analysis/v3/distance_decay_turnover"


def test_sorensen_partition_sums_to_sorensen():
    presence = np.array(
        [
            [1, 1, 1, 0, 0],
            [1, 1, 0, 0, 0],
            [0, 0, 1, 1, 1],
        ],
        dtype=bool,
    )
    sorensen, simpson, nestedness = sorensen_partition(presence)
    assert np.allclose(sorensen, simpson + nestedness)
    assert np.allclose(np.diag(sorensen), 0.0)
    # Row 1 is a strict subset of row 0: pure nestedness, zero replacement.
    assert simpson[0, 1] == pytest.approx(0.0)
    assert nestedness[0, 1] == pytest.approx(sorensen[0, 1])
    assert sorensen[0, 1] > 0


def test_bray_partition_sums_and_is_symmetric():
    abundance = np.array(
        [[10.0, 5.0, 0.0], [2.0, 8.0, 1.0], [0.0, 0.0, 20.0]]
    )
    bray, balanced, gradient = bray_partition(abundance)
    assert np.allclose(bray, bray.T)
    assert np.allclose(bray, balanced + gradient)
    assert np.all(bray >= -1e-12)
    assert np.all(bray <= 1 + 1e-12)


def test_bray_abundance_gradient_vanishes_at_equal_totals():
    # This is why the presence-absence partition is the reportable arm:
    # coverage standardisation drives the abundance-gradient term to zero.
    rng = np.random.default_rng(0)
    abundance = rng.multinomial(1000, [0.2, 0.3, 0.5], size=4).astype(float)
    _, _, gradient = bray_partition(abundance)
    assert np.allclose(gradient, 0.0, atol=1e-12)


def test_permutation_slope_recovers_a_planted_decay():
    rng = np.random.default_rng(11)
    coordinates = rng.normal(size=(25, 2)) * 50
    geographic = euclidean_matrix(pd.DataFrame(coordinates))
    response = {"planted": upper_triangle(0.3 * geographic)}
    observed, null = permutation_slopes(geographic, response, 199, seed=3)
    assert observed["planted"] == pytest.approx(0.3, abs=1e-9)
    assert two_sided_p(observed["planted"], null["planted"]) <= 0.01


def test_two_sided_p_is_centred_on_the_null_mean():
    null = np.full(99, 5.0)
    assert two_sided_p(5.0, null) == pytest.approx(1.0)
    assert two_sided_p(9.0, null) == pytest.approx(0.01)


def test_asv_counts_follow_the_reference_cohort_exactly(tmp_path):
    """Groups come from the genus cohort; low-depth profiles are excluded."""
    # e0001_1Sr1 and e0002_1Sr2 form campaign-1 site-1 Surface; e0003_1Sr3
    # is a profile the genus cache dropped (below the cleaning threshold);
    # e0004_T2Dr1 is a group the genus cohort did not retain.
    asv = pd.DataFrame(
        {
            "e0001_1Sr1": [10, 0, 5],
            "e0002_1Sr2": [1, 2, 3],
            "e0003_1Sr3": [100, 100, 100],
            "e0004_T2Dr1": [7, 7, 7],
            "e0005_5Dr1": [0, 1, 0],
        },
        index=["asv_a", "asv_b", "asv_c"],
    )
    asv_path = tmp_path / "asv.tsv"
    asv.to_csv(asv_path, sep="\t")
    reference = pd.DataFrame(
        {
            "campaign": [1, 1],
            "site": [1, 5],
            "compartment": ["Surface", "Deep"],
        }
    )
    reference_samples = {"e0001_1Sr1", "e0002_1Sr2", "e0005_5Dr1", "e0004_T2Dr1"}
    counts, metadata, info = asv_counts_on_reference_groups(
        asv_path, reference_samples, reference
    )
    assert counts.shape == (3, 2)
    assert list(counts.index) == ["asv_a", "asv_b", "asv_c"]
    # The excluded profile's reads never enter the group sum.
    assert counts.iloc[:, 0].tolist() == [11, 2, 8]
    assert counts.iloc[:, 1].tolist() == [0, 1, 0]
    assert metadata[["campaign", "site", "compartment"]].equals(reference)
    assert info["reference_groups"] == 2
    assert info["aligned_groups"] == 2
    assert info["profiles_excluded_as_absent_from_genus_cache"] == ["e0003_1Sr3"]
    assert info["profiles_used"] == 4
    assert info["asv_cache_groups_outside_reference_cohort"] == [
        "campaign 2, site 2, Deep"
    ]
    assert info["filtered_asv_reads_per_group_min"] == 1
    assert set(info["reference_groups_below_2000_filtered_asv_reads"]) == {
        "campaign 1, site 1, Surface",
        "campaign 1, site 5, Deep",
    }


def test_asv_counts_refuse_a_reference_group_without_profiles(tmp_path):
    asv = pd.DataFrame({"e0001_1Sr1": [1, 2]}, index=["asv_a", "asv_b"])
    asv_path = tmp_path / "asv.tsv"
    asv.to_csv(asv_path, sep="\t")
    reference = pd.DataFrame(
        {"campaign": [1, 1], "site": [1, 2], "compartment": ["Surface", "Surface"]}
    )
    with pytest.raises(ValueError, match="no ASV profile"):
        asv_counts_on_reference_groups(asv_path, {"e0001_1Sr1"}, reference)


def test_reference_sample_ids_are_the_genus_cache_columns(tmp_path):
    genus_path = tmp_path / "genus.tsv"
    pd.DataFrame({"s1": [1], "s2": [2]}, index=["g"]).to_csv(genus_path, sep="\t")
    assert reference_sample_ids(genus_path) == {"s1", "s2"}


def test_site_coordinate_table_must_hold_the_sixty_core_sites(tmp_path):
    columns = ["site", "latitude", "longitude", "x_km", "y_km", "transect_km"]
    good = pd.DataFrame(
        {
            "site": range(60, 0, -1),
            **{column: np.linspace(0, 1, 60) for column in columns[1:]},
        }
    )
    path = tmp_path / "coordinates.tsv"
    good.to_csv(path, sep="\t", index=False)
    loaded = load_site_coordinates(path)
    assert loaded["site"].tolist() == list(range(1, 61))
    good.iloc[:-1].to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="60 core sites"):
        load_site_coordinates(path)
    good.drop(columns=["transect_km"]).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="transect_km"):
        load_site_coordinates(path)


@pytest.mark.skipif(
    not (ROOT / DEFAULT_COORDINATES).exists(),
    reason="corrected coordinate table not staged",
)
def test_default_coordinates_carry_the_site_52_correction():
    coordinates = load_site_coordinates(ROOT / DEFAULT_COORDINATES)
    site_52 = coordinates.set_index("site").loc[52]
    assert site_52["latitude"] == pytest.approx(20.82784, abs=1e-5)
    assert site_52["longitude"] == pytest.approx(53.57835, abs=1e-5)
    # Site 52 must no longer share the uncorrected campaign-1/3 position of
    # site 53.
    site_53 = coordinates.set_index("site").loc[53]
    assert abs(site_52["longitude"] - site_53["longitude"]) > 0.1


@pytest.mark.skipif(
    not (ASV_DIR / "asv_resolution_sensitivity.tsv").exists(),
    reason="ASV sensitivity not staged",
)
def test_asv_resolution_supports_the_genus_primary_result():
    frame = pd.read_csv(ASV_DIR / "asv_resolution_sensitivity.tsv", sep="\t")
    genus = frame[frame["resolution"] == "genus"].iloc[0]
    asv = frame[frame["resolution"] == "asv"]
    assert len(asv) >= 2
    assert genus["n_groups"] == 630 and genus["n_sites"] == 60
    # Both arms run on the primary cohort: the same 630 groups at 60 sites.
    assert (asv["n_groups"] == 630).all() and (asv["n_sites"] == 60).all()
    assert (asv["partial_r2"] > 0.5 * genus["partial_r2"]).all()
    assert (asv["permutation_p"] < 0.05).all()
    assert (asv["partial_r2_ci_low"] < asv["partial_r2"]).all()
    assert (asv["partial_r2_ci_high"] > asv["partial_r2"]).all()


@pytest.mark.skipif(
    not (ASV_DIR / "claim_verdict.json").exists(),
    reason="ASV sensitivity not staged",
)
def test_genus_arm_reproduces_the_corrected_primary_fit():
    """The genus comparison value is the corrected-coordinate primary fit."""
    frame = pd.read_csv(ASV_DIR / "asv_resolution_sensitivity.tsv", sep="\t")
    verdict = json.loads((ASV_DIR / "claim_verdict.json").read_text())
    primary = json.loads(
        (ROOT / "analysis/v3/spatial_turnover_rescue/results/claim_verdict.json")
        .read_text()
    )
    genus = frame[frame["resolution"] == "genus"].iloc[0]
    assert genus["partial_r2"] == pytest.approx(
        primary["primary_partial_r2"], abs=1e-9
    )
    assert verdict["genus_primary_partial_r2"] == pytest.approx(
        primary["primary_partial_r2"], abs=1e-9
    )
    assert verdict["genus_primary_partial_r2_95_jackknife_ci"] == pytest.approx(
        primary["primary_partial_r2_95_jackknife_ci"], abs=1e-9
    )
    assert genus["residual_moran_i"] == pytest.approx(
        primary["primary_residual_moran_i"], abs=1e-9
    )
    assert verdict["input"]["coordinates_path"].endswith(
        str(DEFAULT_COORDINATES)
    )


@pytest.mark.skipif(
    not (ASV_DIR / "claim_verdict.json").exists(),
    reason="ASV sensitivity not staged",
)
def test_asv_cohort_is_the_full_primary_cohort_and_fully_described():
    frame = pd.read_csv(ASV_DIR / "asv_resolution_sensitivity.tsv", sep="\t")
    verdict = json.loads((ASV_DIR / "claim_verdict.json").read_text())
    cohort = verdict["asv_cohort"]
    assert cohort["reference_groups"] == 630
    assert cohort["aligned_groups"] == 630
    asv_groups = set(frame[frame["resolution"] == "asv"]["n_groups"])
    assert asv_groups == {630}
    # The ASV cache's extra material is named, not silently absorbed: nine
    # core-site profiles below the cleaning threshold, and three groups whose
    # retained profiles fall short of 2,000 genus-assigned reads.
    assert len(cohort["profiles_excluded_as_absent_from_genus_cache"]) == 9
    assert cohort["asv_cache_groups_outside_reference_cohort"] == [
        "campaign 1, site 46, Rhizosphere",
        "campaign 2, site 1, Deep",
        "campaign 4, site 53, Deep",
    ]
    assert cohort["reference_groups_below_2000_filtered_asv_reads"] == {
        "campaign 3, site 54, Rhizosphere": 894
    }

    wording = verdict["permitted_wording"]
    assert "same 630 site-campaign-compartment groups" in wording
    assert "corrected site coordinates" in wording
    for stale in ("archived", "intersection", "earlier coordinate", "629"):
        assert stale not in wording


@pytest.mark.skipif(
    not (ASV_DIR / "moran_k_sensitivity.tsv").exists(),
    reason="Moran k sensitivity not staged",
)
def test_moran_k_sensitivity_is_reported_and_bounded():
    frame = pd.read_csv(ASV_DIR / "moran_k_sensitivity.tsv", sep="\t")
    assert list(frame["neighbours_k"]) == list(NEIGHBOUR_COUNTS)
    assert frame["residual_moran_i"].is_monotonic_decreasing
    verdict = json.loads((ASV_DIR / "claim_verdict.json").read_text())
    detected = verdict["neighbour_counts_with_detected_autocorrelation"]
    undetected = verdict["neighbour_counts_without_detected_autocorrelation"]
    assert sorted(detected + undetected) == sorted(NEIGHBOUR_COUNTS)
    if undetected:
        assert verdict["moran_k_status"] == (
            "residual_autocorrelation_depends_on_neighbour_count"
        )
        assert "neighbour count" in verdict["prohibited_wording"]


@pytest.mark.skipif(
    not (DECAY_DIR / "distance_decay_slopes.tsv").exists(),
    reason="distance decay not staged",
)
def test_distance_decay_uses_paired_contrasts_and_whole_site_permutations():
    frame = pd.read_csv(DECAY_DIR / "distance_decay_slopes.tsv", sep="\t")
    verdict = json.loads((DECAY_DIR / "claim_verdict.json").read_text())
    aitchison = frame[frame["family"] == "aitchison"]
    assert set(aitchison["response"]) == {"Surface", "Deep", "Rhizosphere"}
    assert (aitchison["slope_per_100km"] > 0).all()
    contrasts = frame[frame["family"] == "contrast"]
    assert set(contrasts["response"]) == {
        "Surface-Rhizosphere",
        "Deep-Rhizosphere",
    }
    # Family-wise control must never be more permissive than the raw test.
    assert (
        contrasts["max_t_adjusted_p"] >= contrasts["two_sided_p"] - 1e-12
    ).all()
    assert verdict["matched_sites"] == 60
    assert verdict["site_pairs"] == 60 * 59 // 2
    assert verdict["permutations"] >= 999
    assert "independent observations" in verdict["prohibited_wording"]


@pytest.mark.skipif(
    not (DECAY_DIR / "distance_decay_pairs.tsv").exists(),
    reason="distance-decay display data not staged",
)
def test_distance_decay_display_data_preserve_the_full_matched_design():
    pairs = pd.read_csv(DECAY_DIR / "distance_decay_pairs.tsv", sep="\t")
    slopes = pd.read_csv(DECAY_DIR / "distance_decay_slopes.tsv", sep="\t")
    expected_pairs = 60 * 59 // 2
    assert len(pairs) == 3 * expected_pairs
    assert set(pairs["compartment"]) == {"Surface", "Deep", "Rhizosphere"}
    assert (pairs.groupby(["site_a", "site_b"]).size() == 3).all()
    assert (pairs["geographic_distance_km"] > 0).all()

    expected = slopes[slopes["family"] == "aitchison"].set_index("response")
    for compartment, part in pairs.groupby("compartment"):
        fitted = np.polyfit(
            part["geographic_distance_km"],
            part["aitchison_dissimilarity"],
            deg=1,
        )[0]
        assert fitted == pytest.approx(
            expected.loc[compartment, "slope_per_km"], abs=1e-10
        )


@pytest.mark.skipif(
    not (DECAY_DIR / "turnover_nestedness_components.tsv").exists(),
    reason="turnover components not staged",
)
def test_turnover_components_are_coverage_standardised():
    frame = pd.read_csv(
        DECAY_DIR / "turnover_nestedness_components.tsv", sep="\t"
    )
    assert len(frame) == 3
    assert frame["standardised_depth"].nunique() == 1
    assert np.allclose(
        frame["mean_sorensen"],
        frame["mean_simpson_turnover"] + frame["mean_nestedness_resultant"],
    )
    assert (frame["turnover_share_of_sorensen"] > 0.5).all()
    # The unstandardised abundance-gradient term tracks library size and is
    # therefore not reportable on its own.
    assert (
        frame["abundance_gradient_vs_log_library_ratio_pearson_r"] > 0.5
    ).all()
