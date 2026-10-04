"""Regression checks for calendar invariance and fully refitted nuisance terms."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import patsy
import pytest

DIRECTORY = Path(__file__).resolve().parents[1] / "analysis/v3"
sys.path.insert(0, str(DIRECTORY))
import rain_calendar_refit as m


def synthetic():
    rng = np.random.default_rng(7)
    rows = [(site, trip, kind) for site in range(1, 7) for trip in range(1, 5) for kind in ["Deep", "Surface"]
            if not (site == 2 and trip == 3)]
    frame = pd.DataFrame(rows, columns=["Site", "Trip", "Type"])
    frame["longitude_z"] = rng.normal(size=len(frame))
    frame["sampling_day_z"] = rng.normal(size=len(frame))
    z = np.asarray(patsy.dmatrix(m.pulse.GROUPED_NUISANCE_FORMULA, frame))
    x = rng.normal(size=(len(frame), 5)) + z @ rng.normal(size=(z.shape[1], 5))
    y = rng.normal(size=(len(frame), 3)) + z @ rng.normal(size=(z.shape[1], 3))
    return frame, z, x, y


def brute_fit(z, x, y, weights):
    root = np.sqrt(weights)[:, None]
    wz, wx, wy = root * z, root * x, root * y
    rx = wx - wz @ np.linalg.lstsq(wz, wx, rcond=None)[0]
    ry = wy - wz @ np.linalg.lstsq(wz, wy, rcond=None)[0]
    xx, yy, xy = (rx * rx).sum(axis=0), (ry * ry).sum(axis=0), rx.T @ ry
    return xy / xx[:, None], xy / np.sqrt(xx[:, None] * yy)


@pytest.mark.parametrize("site_weights", [np.ones(6), np.array([0, 2, 1, 0, 1, 2]), np.array([4, 0, 0, 1, 1, 0])])
def test_reweighted_fwl_matches_raw_nuisance_refit(site_weights):
    frame, z, x, y = synthetic()
    design = m.RefitDesign(frame, y)
    result = design.fit(design.prepare(x), site_weights)
    beta, rho = brute_fit(z, x, y, site_weights[frame.Site.to_numpy() - 1])
    np.testing.assert_allclose(result["beta"], beta, atol=1e-8)
    np.testing.assert_allclose(result["rho"], rho, atol=1e-8)


def test_refit_is_not_resampling_fixed_nuisance_residuals():
    frame, z, x, y = synthetic()
    weights = np.array([0, 2, 1, 0, 1, 2])[frame.Site.to_numpy() - 1]
    beta, _ = brute_fit(z, x, y, weights)
    rx = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
    ry = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    frozen = (rx.T @ (weights[:, None] * ry)) / ((weights[:, None] * rx ** 2).sum(axis=0)[:, None])
    assert np.max(np.abs(beta - frozen)) > .01


def test_calendar_mapping_preserves_shared_year_and_december_boundary():
    dates = pd.to_datetime(["2023-03-17", "2023-07-14", "2023-12-20", "2024-02-17"])
    mapped = m.map_calendar_dates(dates, {2023: 2025, 2024: 2023})
    assert list(mapped.strftime("%Y-%m-%d")) == ["2025-03-17", "2025-07-14", "2022-12-20", "2023-02-17"]
    with pytest.raises(ValueError, match="Leap-day"):
        m.map_calendar_dates(pd.to_datetime(["2024-02-29"]), {2024: 2023})


def test_first_year_requires_previous_december():
    observations = pd.DataFrame({"Date": pd.to_datetime(["2023-03-17", "2024-02-17", "2025-10-15"])})
    class Weather:
        values = pd.DataFrame(index=pd.date_range("2022-01-01", "2026-01-20"))
    assert m.eligible_years(observations, [Weather()]) == [2023, 2024, 2025]


def test_exhaustive_null_includes_identity_and_ties():
    assert m.exhaustive_p([9, 8, 7, 6, 5, 4]) == 1 / 6
    assert m.exhaustive_p([9, 9, 7, 6, 5, 4]) == 2 / 6
    assert m.exhaustive_p([9, 8, 7, 6, 5, 4], 5) == 1
