import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis" / "v3"))

from geographic_prediction import (  # noqa: E402
    CAMPAIGN_SAMPLESHEETS,
    collection_order_alias,
    contiguous_blocks,
    cyclic_reflection_maps,
    load_collection_order,
    training_taxa,
)

RESULTS = ROOT / "analysis/v3/geographic_prediction"
ALIAS = RESULTS / "collection_order_alias.tsv"
NULLS = RESULTS / "prediction_nulls.tsv"
# prediction_folds.tsv holds the primary joint campaign-by-block arm;
# site_level_block_folds.tsv holds the site-block-only sensitivity.
FOLDS = RESULTS / "prediction_folds.tsv"
SITE_FOLDS = RESULTS / "site_level_block_folds.tsv"
VERDICT = RESULTS / "claim_verdict.json"
README = RESULTS / "README.md"


def test_contiguous_blocks_are_equal_sized_and_ordered():
    blocks = contiguous_blocks(np.arange(60), 6)
    assert len(blocks) == 60
    assert list(np.unique(blocks)) == list(range(6))
    assert np.all(np.diff(blocks) >= 0)
    assert list(np.bincount(blocks)) == [10] * 6


def test_cyclic_reflection_maps_are_bijections_and_exclude_identity():
    sites = list(range(1, 11))
    maps = cyclic_reflection_maps(sites)
    assert len(maps) == 2 * len(sites) - 1
    for mapping in maps:
        assert sorted(mapping) == sites
        assert sorted(mapping.values()) == sites
        assert any(key != value for key, value in mapping.items())


def test_training_taxa_never_look_at_held_out_columns():
    # A taxon that is abundant only in the held-out columns must not be
    # selected, otherwise feature selection leaks.
    counts = pd.DataFrame(
        {
            "g0": [100, 100, 100, 0],
            "g1": [100, 100, 100, 0],
            "g2": [1, 1, 1, 1000],
        },
        index=["train_a", "train_b", "train_c", "holdout"],
    ).T
    training = np.array([True, True, True, False])
    selected = training_taxa(counts, None, training, 0.5, 2)
    assert "g2" not in selected
    assert set(selected) == {"g0", "g1"}


def test_load_collection_order_reads_every_campaign():
    order = load_collection_order(ROOT)
    assert set(order["campaign"]) == set(CAMPAIGN_SAMPLESHEETS)
    # trip2-2023.tsv has a nine-name header over ten-field rows; a name-based
    # read silently loses it, so guard the positional read explicitly.
    assert (order["campaign"] == 2).sum() >= 3
    assert order["site"].between(1, 60).all()


@pytest.mark.skipif(not ALIAS.exists(), reason="alias table not staged")
def test_collection_order_is_aliased_with_transect_position():
    alias = pd.read_csv(ALIAS, sep="\t")
    assert len(alias) == len(CAMPAIGN_SAMPLESHEETS)
    assert (alias["abs_spearman_rho"] >= 0.99).all()
    assert alias["n_sites_with_timestamp"].max() == 60


def test_collection_order_alias_is_recomputable():
    from spatial_turnover_rescue import load_coordinates

    rows, summary = collection_order_alias(
        load_collection_order(ROOT), load_coordinates(ROOT)
    )
    assert summary["alias_status"] == (
        "collection_order_aliased_with_transect_position"
    )
    assert summary["min_abs_spearman_rho"] >= 0.99
    assert summary["campaigns_with_abs_rho_at_least_0_99"] == len(rows)


@pytest.mark.skipif(not VERDICT.exists(), reason="verdict not staged")
def test_primary_arm_is_the_joint_campaign_and_block_holdout():
    """The requested design excludes a campaign AND a block together."""
    verdict = json.loads(VERDICT.read_text())
    assert verdict["primary_arm"] == "group_level_campaign_by_block"
    assert verdict["sensitivity_arm"] == "site_level_block"
    definition = verdict["arm_definitions"]["group_level_campaign_by_block"]
    assert "whole campaign" in definition and "contiguous transect" in definition
    sensitivity = verdict["arm_definitions"]["site_level_block"]
    assert "cannot test transport to an unseen campaign" in sensitivity


@pytest.mark.skipif(not VERDICT.exists(), reason="verdict not staged")
def test_overall_status_is_decided_by_the_joint_arm():
    verdict = json.loads(VERDICT.read_text())
    folds = pd.read_csv(FOLDS, sep="\t")
    nulls = pd.read_csv(NULLS, sep="\t")
    joint_nulls = nulls[nulls["arm"] == "group_level_campaign_by_block"]
    assert set(joint_nulls["null"]) == {
        "whole_site_relabelling",
        "cyclic_shift_reflection",
    }
    strict = joint_nulls["p_value"].max()
    assert verdict["strictest_group_level_null_p_value"] == pytest.approx(strict)

    positive = int((folds["fold_skill_r2"] > 0).sum())
    joint_supported = bool(
        verdict["group_level_equal_weight_skill"] > 0
        and positive == len(folds)
        and strict < 0.05
    )
    assert verdict["primary_arm_supported"] is joint_supported
    assert verdict["status"].startswith(
        "joint_campaign_block_supported"
        if joint_supported
        else "joint_campaign_block_not_supported"
    )
    # The site-block sensitivity must never appear before the joint verdict
    # in the status string.
    assert verdict["status"].index("joint_campaign_block") == 0


@pytest.mark.skipif(not VERDICT.exists(), reason="verdict not staged")
def test_joint_arm_currently_fails_and_is_reported_as_failing():
    verdict = json.loads(VERDICT.read_text())
    folds = pd.read_csv(FOLDS, sep="\t")
    assert len(folds) == 18
    assert verdict["group_level_equal_weight_skill"] < 0.05
    assert verdict["group_level_pooled_skill"] < 0.05
    assert int((folds["fold_skill_r2"] > 0).sum()) < len(folds)
    assert verdict["strictest_group_level_null_p_value"] > 0.05
    assert verdict["primary_arm_supported"] is False
    assert verdict["status"] == (
        "joint_campaign_block_not_supported_site_block_sensitivity_supported"
    )
    assert "did not support cross-campaign, cross-block transport" in (
        verdict["permitted_wording"]
    )


@pytest.mark.skipif(not VERDICT.exists(), reason="verdict not staged")
def test_site_block_arm_is_labelled_a_sensitivity_not_the_result():
    verdict = json.loads(VERDICT.read_text())
    folds = pd.read_csv(SITE_FOLDS, sep="\t")
    nulls = pd.read_csv(NULLS, sep="\t")
    assert len(folds) == 6
    assert folds["n_heldout_sites"].sum() == 60
    site_nulls = nulls[nulls["arm"] == "site_level_block"]
    assert verdict["strictest_site_level_null_p_value"] == pytest.approx(
        site_nulls["p_value"].max()
    )
    assert verdict["sensitivity_arm_supported"] is True
    wording = verdict["permitted_wording"]
    # Order matters: the joint arm is stated first.
    assert wording.index("primary test") < wording.index("As a sensitivity")
    assert "does not demonstrate transport to an unseen campaign" in wording


@pytest.mark.skipif(not VERDICT.exists(), reason="verdict not staged")
def test_prohibited_wording_forbids_calling_prediction_successful():
    verdict = json.loads(VERDICT.read_text())
    prohibited = verdict["prohibited_wording"]
    assert "do not describe geographic prediction as" in prohibited.lower()
    assert "site-block-only sensitivity" in prohibited
    assert "collection order" in prohibited
    assert "succeeded" not in verdict["permitted_wording"]
