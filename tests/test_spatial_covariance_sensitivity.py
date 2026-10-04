import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis/v3"))
import spatial_covariance_sensitivity as sensitivity
import render_spatial_covariance_tex as renderer


def test_unequal_group_centering_operator_matches_direct_computation():
    metadata = pd.DataFrame({"site": [1, 2, 3, 1, 3, 2, 3],
                             "campaign": [1, 1, 1, 2, 2, 3, 3],
                             "compartment": ["S"] * 7})
    values = np.arange(21.).reshape(7, 3)
    transform, shared, noise = sensitivity.averaging_operators(metadata, np.arange(1, 4))
    centered = values.copy()
    for indices in metadata.groupby(["campaign", "compartment"]).indices.values():
        centered[indices] -= centered[indices].mean(axis=0)
    expected = np.vstack([centered[metadata.site == site].mean(axis=0) for site in (1, 2, 3)])
    np.testing.assert_allclose(transform @ values, expected)
    np.testing.assert_allclose(noise, transform @ transform.T)
    np.testing.assert_allclose(shared @ np.ones(3), 0, atol=1e-14)


def test_covariance_whitening_and_nuisance_projection():
    x = np.linspace(0, 1, 20)
    covariance = .8 * np.exp(-np.abs(x[:, None] - x[None, :]) / .3) + .2 * np.eye(20)
    design = np.column_stack([np.ones(20), x, x * x])
    transform, tested = sensitivity.null_projection(covariance, design)
    np.testing.assert_allclose(transform @ covariance @ transform.T, np.eye(19), atol=1e-12)
    np.testing.assert_allclose(transform @ np.ones(20), 0, atol=1e-12)
    np.testing.assert_allclose(tested.T @ tested, np.eye(2), atol=1e-12)


def test_scalar_rotation_matches_exact_gls_f_distribution():
    rng = np.random.default_rng(193)
    x = np.linspace(-1, 1, 30)
    covariance = .7 * np.exp(-np.abs(x[:, None] - x[None, :]) / .2) + .3 * np.eye(30)
    design = np.column_stack([np.ones(30), x, x * x])
    response = linalg.cholesky(covariance, lower=True) @ rng.normal(size=(30, 1))
    transform, tested = sensitivity.null_projection(covariance, design)
    result = sensitivity.rotation_test(response, transform, tested, 19999, rng)
    r2 = result["whitened_route_r2"]
    exact = stats.f.sf((r2 / 2) / ((1 - r2) / 27), 2, 27)
    assert abs(exact - result["rotation_p"]) < .025


def test_rotation_preserves_arbitrary_cross_response_dependence():
    rng = np.random.default_rng(17)
    values = rng.normal(size=(13, 5)) @ np.diag([4., 3., 1., 0., 0.])
    rotation = np.linalg.qr(rng.normal(size=(13, 13)))[0]
    np.testing.assert_allclose((rotation @ values).T @ (rotation @ values), values.T @ values, atol=1e-12)


def test_known_centering_constraint_is_whitened_without_jitter():
    n = 10
    center = np.eye(n) - np.ones((n, n)) / n
    x = np.linspace(-1, 1, n)
    design = center @ np.column_stack([np.ones(n), x, x * x])
    transform, tested = sensitivity.null_projection(center, design)
    assert transform.shape == (n - 1, n)
    assert tested.shape == (n - 1, 2)
    np.testing.assert_allclose(transform @ center @ transform.T, np.eye(n - 1), atol=1e-12)
    square_root = sensitivity.covariance_square_root(center)
    np.testing.assert_allclose(square_root @ square_root.T, center, atol=1e-12)


def test_release_grid_and_renderer_retain_every_model_and_calibration_case():
    results = ROOT / "analysis/v3/spatial_covariance_sensitivity_20260909"
    table = pd.read_csv(results / "spatial_covariance_tests.tsv", sep="\t")
    assert len(table) == 74
    assert table.model.is_unique
    assert table.iloc[1:].covariance_support_rank.eq(59).all()
    assert table.tested_rank.eq(2).all()
    text = renderer.render(results)
    for label in renderer.CASES.values():
        assert label in text
    for kind in sensitivity.KERNELS:
        rows = table[table.kernel == kind]
        assert len(rows) == 24
        assert rows.groupby("range_km").size().eq(4).all()
    assert "Gaussian reference is omitted for heavy-tailed" in text
    assert text.count("Not evaluated") == 2
