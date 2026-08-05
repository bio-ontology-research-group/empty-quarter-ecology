from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis" / "v3"))

import spatial_turnover_rescue as spatial


def test_sample_metadata_and_trip1_only_exclusion():
    assert spatial.sample_metadata("e0325_10Dr2") == {
        "sample_id": "e0325_10Dr2",
        "campaign": 1,
        "site": 10,
        "compartment": "Deep",
    }
    assert spatial.sample_metadata("F61PRr3")["site"] == 61


def test_site_level_clr_centres_campaign_compartment_groups():
    columns = pd.MultiIndex.from_tuples(
        [
            (1, 1, "Surface"),
            (1, 2, "Surface"),
            (2, 1, "Surface"),
            (2, 2, "Surface"),
        ],
        names=["campaign", "site", "compartment"],
    )
    counts = pd.DataFrame(
        [[10, 20, 30, 40], [40, 30, 20, 10], [5, 5, 5, 5]],
        index=["a", "b", "c"],
        columns=columns,
    )
    metadata = columns.to_frame(index=False)
    sites, response, n_groups = spatial.site_level_clr(
        counts, metadata, ["a", "b", "c"], None, 0.5
    )
    assert sites.tolist() == [1, 2]
    assert response.shape == (2, 3)
    assert n_groups == 4
    assert np.allclose(response.sum(axis=1), 0)
    assert np.allclose(response.mean(axis=0), 0)


def test_multivariate_model_detects_known_trend():
    rng = np.random.default_rng(4)
    transect = np.linspace(-2, 2, 30)
    design = spatial.design_matrix(transect, 2)
    response = np.column_stack(
        [
            1.5 * transect + rng.normal(0, 0.1, len(transect)),
            -0.5 * transect**2 + rng.normal(0, 0.1, len(transect)),
        ]
    )
    r2, pseudo_f, residual = spatial.fit_multivariate(response, design)
    assert r2 > 0.95
    assert pseudo_f > 100
    assert residual.shape == response.shape


def test_decision_prohibits_process_inference():
    rows = [
        spatial.SpatialResult(
            analysis="primary",
            omitted_campaign=None,
            taxon_count=200,
            trend_degree=2,
            n_sites=60,
            n_groups=100,
            partial_r2=0.1,
            pseudo_f=3.0,
            permutation_p=0.001,
            residual_moran_i=0.01,
            residual_moran_p=0.4,
        )
    ]
    for campaign in range(1, 6):
        rows.append(
            spatial.SpatialResult(
                analysis="leave_one_campaign_out",
                omitted_campaign=campaign,
                taxon_count=200,
                trend_degree=2,
                n_sites=60,
                n_groups=80,
                partial_r2=0.1,
                pseudo_f=3.0,
                permutation_p=0.01,
                residual_moran_i=0.01,
                residual_moran_p=0.4,
            )
        )
    for taxon_count, degree in ((80, 1), (80, 2), (200, 1), (500, 1), (500, 2)):
        rows.append(
            spatial.SpatialResult(
                analysis="primary",
                omitted_campaign=None,
                taxon_count=taxon_count,
                trend_degree=degree,
                n_sites=60,
                n_groups=100,
                partial_r2=0.1,
                pseudo_f=3.0,
                permutation_p=0.01,
                residual_moran_i=0.01,
                residual_moran_p=0.4,
            )
        )
    verdict = spatial.decision(rows)
    assert verdict["status"] == "broad_geographic_structure_supported"
    assert "not dispersal limitation" in verdict["permitted_wording"]
