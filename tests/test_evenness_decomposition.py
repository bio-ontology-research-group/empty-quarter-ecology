import numpy as np
import pandas as pd

from analysis.v3.evenness_decomposition_analysis import (
    PRIMARY_METRICS,
    aggregate_profiles,
    bh_fdr,
    paired_contrasts,
)


def test_bh_fdr_known_values():
    assert np.allclose(
        bh_fdr([0.01, 0.04, 0.03]),
        [0.03, 0.04, 0.04],
    )


def test_aggregation_preserves_site_as_final_unit():
    rows = []
    for site in range(1, 4):
        for trip in (1, 2):
            for position, offset in (
                ("Surface", 0.0),
                ("Deep", 0.1),
                ("Rhizosphere", -0.2),
            ):
                for replicate in (1, 2):
                    rows.append(
                        {
                            "Trip": trip,
                            "Site": site,
                            "Type": position,
                            "depth": 100 + replicate,
                            "shannon": 4.0 + offset,
                            "richness_hurlbert_25000": 500 + 10 * offset,
                            "evenness_h_over_log_hurlbert": 0.8 + offset,
                            "evenness_h_over_log_observed": 0.78 + offset,
                        }
                    )
    frame = pd.DataFrame(rows)
    blocks, site_means = aggregate_profiles(frame)
    assert len(blocks) == 18
    assert len(site_means) == 9
    assert set(site_means["Site"]) == {1, 2, 3}


def test_paired_contrasts_use_three_site_pairs():
    rows = []
    for site in range(1, 4):
        for position, offset in (
            ("Surface", 0.0),
            ("Deep", 0.1),
            ("Rhizosphere", -0.2),
        ):
            row = {"Site": site, "Type": position}
            for metric in PRIMARY_METRICS:
                row[metric] = site + offset
            rows.append(row)
    result = paired_contrasts(pd.DataFrame(rows), PRIMARY_METRICS)
    assert set(result["n_sites"]) == {3}
    root_surface = result[
        (result["metric"] == "shannon")
        & (result["contrast"] == "Rhizosphere-Surface")
    ].iloc[0]
    assert np.isclose(root_surface["mean_difference"], -0.2)
