"""Regression checks for nuisance degrees of freedom in rank association."""
import importlib.util
from pathlib import Path

import numpy as np
from scipy import stats
import statsmodels.api as sm

spec = importlib.util.spec_from_file_location(
    "cross_desert_context", Path(__file__).parents[1] / "analysis/v3/cross_desert_context.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_partial_rank_p_matches_regression_with_nuisance_covariates():
    rng = np.random.default_rng(19)
    nuisance = rng.normal(size=(16, 3))
    x = nuisance[:, 0] + rng.normal(size=16)
    y = .8 * x + nuisance[:, 1] + rng.normal(size=16)
    _, p = module.partial_spearman(x, y, nuisance)
    fitted = sm.OLS(stats.rankdata(y), sm.add_constant(
        np.column_stack([stats.rankdata(x), nuisance])
    )).fit()
    assert np.isclose(p, fitted.pvalues[1], rtol=1e-12)


def test_duplicate_nuisance_column_does_not_consume_an_extra_degree():
    rng = np.random.default_rng(23)
    x, y, z = rng.normal(size=(3, 20))
    unique = module.partial_spearman(x, y, z[:, None])
    duplicate = module.partial_spearman(x, y, np.column_stack([z, z]))
    assert np.allclose(unique, duplicate, rtol=1e-12)
