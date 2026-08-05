import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


pytest.importorskip("sklearn")


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis/v3/network_rescue/run_network_rescue.py"
SPEC = importlib.util.spec_from_file_location("network_rescue", SCRIPT)
NETWORK = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = NETWORK
SPEC.loader.exec_module(NETWORK)


def test_sample_metadata_distinguishes_campaign_and_compartment():
    assert NETWORK.sample_metadata("17Sr2") == {
        "sample_id": "17Sr2",
        "campaign": 1,
        "site": 17,
        "compartment": "Surface",
        "replicate": 2,
    }
    assert NETWORK.sample_metadata("V17PRr1")["campaign"] == 5
    assert NETWORK.sample_metadata("V17PRr1")["compartment"] == "Rhizosphere"
    assert NETWORK.sample_metadata("not-a-sample") is None


def test_independent_permutation_preserves_each_campaign_taxon_multiset():
    matrix = np.arange(36, dtype=float).reshape(12, 3)
    campaigns = np.repeat([1, 2, 3], 4)
    permuted = NETWORK.independent_taxon_permutation(
        matrix, campaigns, np.random.default_rng(77)
    )
    for campaign in np.unique(campaigns):
        indices = np.where(campaigns == campaign)[0]
        for column in range(matrix.shape[1]):
            assert sorted(permuted[indices, column]) == sorted(
                matrix[indices, column]
            )


def test_taxa_selection_marks_only_primary_not_sensitivity_only_taxa():
    columns = pd.MultiIndex.from_tuples(
        [(1, 1), (1, 2), (2, 1), (2, 2)],
        names=["campaign", "site"],
    )
    base = pd.DataFrame(
        [
            [40, 41, 42, 43],
            [30, 31, 32, 33],
            [20, 21, 22, 23],
            [10, 11, 12, 13],
        ],
        index=["a", "b", "c", "d"],
        columns=columns,
    )
    grouped = {compartment: base.copy() for compartment in NETWORK.COMPARTMENTS}
    ranking, rows = NETWORK.select_taxa(grouped, 0.5, 4, 2)
    selected = {row["genus"] for row in rows if row["selected_primary"]}
    sensitivity_pool = {
        row["genus"] for row in rows if row["selected_sensitivity_pool"]
    }
    assert ranking == ["a", "b", "c", "d"]
    assert selected == {"a", "b"}
    assert sensitivity_pool == {"a", "b", "c", "d"}


def test_zero_false_fraction_passes_stable_network_gate():
    calibration = {
        "status": "calibrated",
        "selected_alpha": 0.2,
        "selected_combined_null_edge_ratio": 0.03,
    }
    metrics = [
        {
            "stable_edges": 12,
            "expected_false_fraction_among_stable": 0.0,
        }
        for _ in NETWORK.COMPARTMENTS
    ]
    verdict = NETWORK.build_verdict(
        calibration,
        metrics,
        [],
        {1: 10, 2: 2, 3: 10, 4: 10, 5: 4},
        minimum_stable_edges=10,
        maximum_false_fraction=0.1,
    )
    assert verdict["stable_network_gate_passed"] is True
    assert (
        verdict["descriptive_association_status"]
        == "stable_associations_without_robust_density_ordering"
    )


def synthetic_count_table(path):
    rows = {"g1": [], "g2": [], "g3": [], "g4": []}
    columns = []
    for campaign, prefix in ((1, ""), (3, "F")):
        for site in range(1, 7):
            for code in ("S", "D", "P"):
                columns.append(f"{prefix}{site}{code}r1")
                offset = campaign * 7 + site * 3 + ord(code)
                rows["g1"].append(50 + offset % 13)
                rows["g2"].append(35 + (offset * 2) % 17)
                rows["g3"].append(25 + (offset * 3) % 19)
                rows["g4"].append(15 + (offset * 5) % 23)
    pd.DataFrame(rows, index=columns).T.to_csv(path, sep="\t")


def directory_hashes(path):
    return {
        item.name: hashlib.sha256(item.read_bytes()).hexdigest()
        for item in sorted(path.iterdir())
        if item.is_file()
    }


def test_small_analysis_writes_complete_deterministic_artifact_set(tmp_path):
    counts = tmp_path / "counts.tsv"
    synthetic_count_table(counts)
    first = tmp_path / "first"
    second = tmp_path / "second"
    kwargs = dict(
        seed=11,
        minimum_group_reads=1,
        prevalence_threshold=0.25,
        taxa_count=3,
        pseudocount=0.5,
        alpha_grid=(0.2,),
        calibration_null_replicates=2,
        null_replicates=2,
        bootstrap_replicates=2,
        sensitivity_taxa=(3,),
        sensitivity_alpha_multipliers=(1.0,),
        sensitivity_null_replicates=2,
        minimum_stable_edges=1,
    )
    NETWORK.run_analysis(tmp_path, counts, first, **kwargs)
    NETWORK.run_analysis(tmp_path, counts, second, **kwargs)

    expected = {
        "README.md",
        "alpha_calibration.tsv",
        "bootstrap_metrics.tsv",
        "campaign_feasibility.tsv",
        "claim_verdict.json",
        "cohort_accounting.tsv",
        "cohort_groups.tsv",
        "compartment_comparisons.tsv",
        "network_edges.tsv",
        "network_metrics.tsv",
        "network_rescue.tsv",
        "parameters.json",
        "sensitivity_metrics.tsv",
        "taxa_selection.tsv",
    }
    assert set(directory_hashes(first)) == expected
    assert directory_hashes(first) == directory_hashes(second)
    parameters = json.loads((first / "parameters.json").read_text())
    assert parameters["seed"] == 11
    assert parameters["input_sha256"]
    combined_text = "\n".join(
        item.read_text() for item in first.iterdir() if item.is_file()
    ).lower()
    for forbidden in ("key" + "stone", "inter" + "action"):
        assert forbidden not in combined_text
