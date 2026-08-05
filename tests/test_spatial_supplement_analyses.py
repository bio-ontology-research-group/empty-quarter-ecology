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
    NEIGHBOUR_COUNTS,
    align_to_reference_groups,
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


def test_align_to_reference_groups_intersects_on_the_reference_keys():
    reference = pd.DataFrame(
        {
            "campaign": [1, 1, 2],
            "site": [1, 2, 1],
            "compartment": ["Surface", "Surface", "Deep"],
        }
    )
    metadata = pd.DataFrame(
        {
            "campaign": [1, 1, 5],
            "site": [1, 2, 9],
            "compartment": ["Surface", "Surface", "Deep"],
        }
    )
    counts = pd.DataFrame(np.arange(9).reshape(3, 3), index=["a", "b", "c"])
    aligned_counts, aligned_metadata, info = align_to_reference_groups(
        counts, metadata, reference
    )
    assert info["groups_dropped_from_asv_cache"] == 1
    assert info["reference_groups_absent_from_asv_cache"] == 1
    assert info["aligned_groups"] == 2
    assert aligned_counts.shape == (3, 2)
    assert len(aligned_metadata) == 2


@pytest.mark.skipif(
    not (ASV_DIR / "asv_resolution_sensitivity.tsv").exists(),
    reason="ASV sensitivity not staged",
)
def test_asv_resolution_supports_the_genus_primary_result():
    frame = pd.read_csv(ASV_DIR / "asv_resolution_sensitivity.tsv", sep="\t")
    genus = frame[frame["resolution"] == "genus"].iloc[0]
    asv = frame[frame["resolution"] == "asv"]
    assert len(asv) >= 2
    assert genus["n_groups"] == 630
    # The ASV arm runs on the intersected canonical cohort, not on the
    # independently grouped 631-profile cache.
    assert (asv["n_groups"] <= 630).all()
    assert (asv["partial_r2"] > 0.5 * genus["partial_r2"]).all()
    assert (asv["permutation_p"] < 0.05).all()


@pytest.mark.skipif(
    not (ASV_DIR / "claim_verdict.json").exists(),
    reason="ASV sensitivity not staged",
)
def test_asv_group_alignment_is_stated_exactly():
    """The two caches do not cover identical groups; say so precisely."""
    frame = pd.read_csv(ASV_DIR / "asv_resolution_sensitivity.tsv", sep="\t")
    verdict = json.loads((ASV_DIR / "claim_verdict.json").read_text())
    alignment = verdict["asv_alignment"]
    assert alignment["reference_groups"] == 630
    assert alignment["aligned_groups"] == 629
    assert alignment["reference_groups_absent_from_asv_cache"] == 1
    assert alignment["groups_dropped_from_asv_cache"] == 2
    asv_groups = set(frame[frame["resolution"] == "asv"]["n_groups"])
    assert asv_groups == {alignment["aligned_groups"]}

    wording = verdict["permitted_wording"]
    assert "629-group intersection" in wording
    assert "629 of the 630 genus-reference groups" in wording
    # The two arms are not on identical groups, so this claim is forbidden.
    assert "on the same groups" not in wording


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
