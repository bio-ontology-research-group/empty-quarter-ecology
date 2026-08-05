from pathlib import Path

import numpy as np

from analysis.v3.pma_endpoint_analysis import (
    analyse_pma,
    expected_rarefied_richness,
    paired_columns,
    write_outputs,
)


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "relic-dna" / "ASV_table.tsv"


def test_expected_rarefied_richness_known_case():
    # At depth two, the doubleton is certain to occur and the singleton has
    # probability 2/3 of occurring.
    assert np.isclose(
        expected_rarefied_richness(np.array([2, 1]), 2),
        1 + 2 / 3,
    )


def test_pair_parser_ignores_controls_and_requires_matching_aliquots():
    pairs = paired_columns(
        ["C1R1T", "C1R1UT", "C2S1UT", "NEGATIVER", "C2S1T"]
    )
    assert [pair["pair_id"] for pair in pairs] == ["C1R1", "C2S1"]
    assert pairs[0]["treated_sample"] == "C1R1T"
    assert pairs[0]["untreated_sample"] == "C1R1UT"


def test_canonical_pma_endpoints_are_bounded_and_reproduced():
    rows, summary = analyse_pma(COUNTS)
    assert len(rows) == 9
    assert summary["rarefaction"]["depth"] == 123_897
    assert np.isclose(
        summary["richness_endpoint"]["treated_mean"],
        608.495324583152,
    )
    assert np.isclose(
        summary["richness_endpoint"]["untreated_mean"],
        771.81834123638,
    )
    assert (
        summary["richness_endpoint"]["wilcoxon_two_sided_exact"]["p_value"]
        == 0.0078125
    )
    assert (
        summary["shannon_endpoint"]["wilcoxon_two_sided_exact"]["p_value"]
        == 0.65234375
    )
    assert summary["status"] == "paired_endpoints_only"
    assert any(
        "relic-DNA fraction" in item
        for item in summary["not_supported"]
    )


def test_output_bundle_is_byte_deterministic(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    write_outputs(ROOT, COUNTS, first, rarefaction_depth=None)
    write_outputs(ROOT, COUNTS, second, rarefaction_depth=None)
    first_files = sorted(path.name for path in first.iterdir())
    second_files = sorted(path.name for path in second.iterdir())
    assert first_files == second_files
    for name in first_files:
        assert (first / name).read_bytes() == (second / name).read_bytes()
