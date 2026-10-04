#!/usr/bin/env python3
"""Calendar-preserving rainfall inference and fully refitted site bootstrap.

The seasonal null permutes complete meteorological-year fields (December 1
through November 30), jointly across sites, campaigns and both products.
Only years with every required lookback are eligible. The null assumes their
exchangeability conditional on the observed community responses and design.
The bootstrap resamples original whole-site histories and refits nuisance
coefficients before repeating the entire endpoint/peak search.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
import patsy
import scipy
from scipy import linalg, stats

import rain_pulse_response as pulse
import rain_response_window as base

ENDPOINTS = list(pulse.ENDPOINTS)
PEAKS = pulse.PEAK_GRID_DAYS


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def geodata(directory):
    frames = []
    for trip in range(1, 6):
        frame = pd.read_csv(directory / f"trip{trip}_geodata.tsv", sep="\t")
        frame = frame[frame.Site.astype(str).str.fullmatch(r"[0-9]+")].copy()
        frame["Site"] = frame.Site.astype(int)
        frame = frame[frame.Site.between(1, 60)]
        frame["Trip"] = trip
        frame["Date"] = pd.to_datetime(frame.CenterDate)
        frames.append(frame[["Trip", "Site", "Date", "Latitude", "Longitude"]])
    result = pd.concat(frames, ignore_index=True)
    if result.duplicated(["Trip", "Site"]).any():
        raise ValueError("Duplicate site/campaign coordinate key")
    result["sampling_day"] = (result.Date - result.groupby("Trip").Date.transform("min")).dt.days.astype(float)
    result["longitude_z"] = base.zscore_within(result, "Trip", "Longitude")
    result["sampling_day_z"] = base.zscore_within(result, "Trip", "sampling_day")
    return result.sort_values(["Trip", "Site"]).reset_index(drop=True)


def meteorological_year(date):
    return date.year + int(date.month == 12)


def map_calendar_dates(dates, year_map):
    """Preserve month/day; December belongs to the following seasonal year.

    Leap-day collection dates are explicitly rejected rather than shifted.
    Lookbacks are then the 60 actual complete calendar days before collection.
    """
    mapped = []
    for value in dates:
        date = pd.Timestamp(value)
        if date.month == 2 and date.day == 29:
            raise ValueError("Leap-day collection needs an explicit design rule")
        target = year_map[meteorological_year(date)] - int(date.month == 12)
        mapped.append(date.replace(year=int(target)))
    return pd.DatetimeIndex(mapped)


def eligible_years(observations, weather_products):
    first = max(w.values.index.min().year for w in weather_products)
    last = min(w.values.index.max().year for w in weather_products)
    observed = set(map(meteorological_year, observations.Date))
    eligible = []
    for year in range(first, last + 1):
        dates = map_calendar_dates(observations.Date, {x: year for x in observed})
        valid = all(
            dates.min() - pd.Timedelta(days=60) >= w.values.index.min()
            and dates.max() <= w.values.index.max()
            for w in weather_products
        )
        if valid:
            eligible.append(year)
    if not observed.issubset(eligible):
        raise ValueError("Original assignment is outside the complete calendar orbit")
    if len(eligible) < 2:
        raise ValueError("At least two complete seasonal years are required")
    return eligible


def within_group(values, labels):
    result = np.asarray(values, dtype=float).copy()
    for group in np.unique(labels):
        selected = labels == group
        result[selected] -= result[selected].mean(axis=0)
    return result


class RefitDesign:
    """FWL reduction, exact for weights constant within original site histories.

    Site-by-compartment means are invariant to these weights. Campaign and
    campaign-specific route/day nuisance coefficients are refitted per draw.
    Selection uses signed partial correlation, monotone in the original t
    statistic with the shared degrees of freedom across peaks and endpoints.
    """

    def __init__(self, observations, outcomes):
        self.sites, self.site_index = np.unique(observations.Site, return_inverse=True)
        groups = pd.factorize(pd.MultiIndex.from_frame(observations[["Site", "Type"]]))[0]
        self.groups = groups
        z = np.asarray(patsy.dmatrix(
            "0 + C(Trip) + C(Trip):longitude_z + C(Trip):sampling_day_z", observations
        ))
        self.z = within_group(z, groups)
        self.y = within_group(outcomes, groups)
        self.rows = [np.flatnonzero(self.site_index == k) for k in range(len(self.sites))]
        self.zz = np.array([self.z[i].T @ self.z[i] for i in self.rows])
        self.zy = np.array([self.z[i].T @ self.y[i] for i in self.rows])
        self.yy = np.array([np.square(self.y[i]).sum(axis=0) for i in self.rows])

    def prepare(self, exposure):
        x = within_group(exposure, self.groups)
        return {
            "zx": np.array([self.z[i].T @ x[i] for i in self.rows]),
            "xx": np.array([np.square(x[i]).sum(axis=0) for i in self.rows]),
            "xy": np.array([x[i].T @ self.y[i] for i in self.rows]),
        }

    def fit(self, prepared, weights=None):
        if weights is None:
            weights = np.ones(len(self.sites))
        weights = np.asarray(weights, dtype=float)
        if weights.shape != (len(self.sites),) or np.any(weights < 0) or not weights.any():
            raise ValueError("Site weights must be nonnegative with positive total")
        zz = np.einsum("s,sjk->jk", weights, self.zz)
        zy = np.einsum("s,sjk->jk", weights, self.zy)
        zx = np.einsum("s,sjk->jk", weights, prepared["zx"])
        inv = linalg.pinvh(zz, rtol=1e-10)
        xx = weights @ prepared["xx"] - np.sum(zx * (inv @ zx), axis=0)
        yy = weights @ self.yy - np.sum(zy * (inv @ zy), axis=0)
        xy = np.einsum("s,sjk->jk", weights, prepared["xy"]) - zx.T @ inv @ zy
        if np.any(yy <= 0):
            raise ValueError("Endpoint has zero residual variation")
        valid = xx > 1e-10
        beta = np.divide(xy, xx[:, None], out=np.full_like(xy, np.nan), where=valid[:, None])
        rho = np.divide(xy, np.sqrt(np.maximum(xx[:, None] * yy, 0)),
                        out=np.full_like(xy, np.nan), where=valid[:, None])
        return {"beta": beta, "rho": rho, "partial_r2": rho ** 2}


def winner(fit):
    peak, endpoint = np.unravel_index(np.nanargmax(fit["rho"]), fit["rho"].shape)
    return {"endpoint": ENDPOINTS[endpoint], "peak_days": float(PEAKS[peak]),
            "beta": float(fit["beta"][peak, endpoint]),
            "partial_r2": float(fit["partial_r2"][peak, endpoint]),
            "rho": float(fit["rho"][peak, endpoint]),
            "peak_index": int(peak), "endpoint_index": int(endpoint)}


def exhaustive_p(values, observed_index=0, absolute=False):
    a = np.asarray(values)
    if absolute:
        a = np.abs(a)
    return float(np.mean(a >= a[observed_index] - 1e-12))


def calendar_orbit(observations, products, design, years):
    kernels = pulse.kernel_matrix(np.arange(1, 61))
    rows, exposures = [], {}
    for index, permutation in enumerate(itertools.permutations(years)):
        mapping = dict(zip(years, permutation))
        shifted = observations.copy()
        shifted["Date"] = map_calendar_dates(observations.Date, mapping)
        record = {"assignment": index, "year_mapping": json.dumps(mapping, sort_keys=True)}
        for name, weather in products.items():
            x = pulse.daily_lag_matrix(shifted, weather) @ kernels
            fit = design.fit(design.prepare(x))
            record[f"{name}_positive_max"] = float(np.nanmax(fit["rho"]))
            record[f"{name}_absolute_max"] = float(np.nanmax(np.abs(fit["rho"])))
            if all(k == v for k, v in mapping.items()):
                exposures[name] = x
        record["joint_positive_max"] = max(record[f"{n}_positive_max"] for n in products)
        record["joint_absolute_max"] = max(record[f"{n}_absolute_max"] for n in products)
        rows.append(record)
    return pd.DataFrame(rows), exposures


def bootstrap_refits(design, exposures, draws, seed):
    prepared = {k: design.prepare(v) for k, v in exposures.items()}
    original = {k: winner(design.fit(v)) for k, v in prepared.items()}
    rng = np.random.default_rng(seed)
    rows = []
    for draw in range(draws):
        weights = rng.multinomial(len(design.sites), np.full(len(design.sites), 1 / len(design.sites)))
        for product, sufficient in prepared.items():
            fit = design.fit(sufficient, weights)
            selected = winner(fit)
            fixed = original[product]
            i, j = fixed["peak_index"], fixed["endpoint_index"]
            endpoint_peak = int(np.nanargmax(fit["rho"][:, j]))
            rows.append({"draw": draw + 1, "product": product,
                         "unique_sites": int((weights > 0).sum()),
                         "selected_endpoint": selected["endpoint"],
                         "selected_peak_days": selected["peak_days"],
                         "selected_beta": selected["beta"],
                         "selected_partial_r2": selected["partial_r2"],
                         "fixed_beta": float(fit["beta"][i, j]),
                         "fixed_partial_r2": float(fit["partial_r2"][i, j]),
                         "fixed_endpoint_selected_peak_days": float(PEAKS[endpoint_peak])})
    frame = pd.DataFrame(rows)
    summaries = []
    for product, part in frame.groupby("product"):
        for metric in ["fixed_beta", "fixed_partial_r2", "fixed_endpoint_selected_peak_days", "selected_peak_days"]:
            quantiles = part[metric].quantile([.025, .5, .975]).to_numpy()
            summaries.append({"product": product, "quantity": metric,
                              "lower_2_5": quantiles[0], "median": quantiles[1], "upper_97_5": quantiles[2]})
    return frame, pd.DataFrame(summaries), original


def wet_spell_deletion(observations, products, design):
    """Delete complete regional runs of days with positive recorded rainfall.

    A wet spell is a descriptive calendar cluster, with no independence claim.
    Its dates are removed jointly from all grid cells/compartments in a trip.
    """
    kernels = pulse.kernel_matrix(np.arange(1, 61))
    rows = []
    for product, weather in products.items():
        rain = pulse.daily_lag_matrix(observations, weather)
        for trip, cohort in observations.groupby("Trip"):
            start, stop = cohort.Date.min() - pd.Timedelta(days=60), cohort.Date.max() - pd.Timedelta(days=1)
            subset = weather.values.loc[start:stop]
            active = subset.gt(0).any(axis=1)
            spells = active.ne(active.shift(fill_value=False)).cumsum()
            for _, dates in active[active].groupby(spells[active]):
                selected_dates = dates.index
                changed = rain.copy()
                removed = 0.0
                for i, row in cohort.iterrows():
                    lag_dates = row.Date - pd.to_timedelta(np.arange(1, 61), unit="D")
                    selected = lag_dates.isin(selected_dates)
                    removed += float(changed[i, selected].sum())
                    changed[i, selected] = 0
                if removed == 0:
                    continue
                fit = winner(design.fit(design.prepare(changed @ kernels)))
                rows.append({"product": product, "trip": int(trip),
                             "first_wet_day": selected_dates.min().date(),
                             "last_wet_day": selected_dates.max().date(),
                             "wet_days": len(selected_dates), **fit})
    return pd.DataFrame(rows)


def rank_calibration(observations, products, years, seed, replicates=300):
    """Check rank implementation under the stated randomized-year null.

    Communities have shared site, campaign and spatial campaign errors. The
    assignment is uniform on the complete orbit independently of communities.
    This checks the conditional null algorithm, not real-year exchangeability.
    """
    rng = np.random.default_rng(seed)
    coords = observations[["Site", "Longitude", "Latitude"]].drop_duplicates("Site").set_index("Site").sort_index()
    xy = coords.to_numpy()
    dist = np.sqrt(((xy[:, None] - xy[None, :]) ** 2).sum(axis=2))
    chol = linalg.cholesky(np.exp(-dist / 1.5) + .1 * np.eye(len(xy)), lower=True)
    site_index = observations.Site.to_numpy(int) - 1
    trip_index = observations.Trip.to_numpy(int) - 1
    # Build fixed residualized exposures once; synthetic responses share design.
    z = patsy.dmatrix(pulse.GROUPED_NUISANCE_FORMULA, observations)
    basis, _ = pulse.load_linear_algebra_basis(z)
    candidate = []
    kernels = pulse.kernel_matrix(np.arange(1, 61))
    for permutation in itertools.permutations(years):
        shifted = observations.copy()
        shifted["Date"] = map_calendar_dates(observations.Date, dict(zip(years, permutation)))
        x = np.column_stack([pulse.daily_lag_matrix(shifted, w) @ kernels for w in products.values()])
        residual = pulse.residualize(x, basis)
        candidate.append(residual / np.sqrt(np.square(residual).sum(axis=0)))
    candidate = np.asarray(candidate)
    records = []
    for replicate in range(replicates):
        site_error = rng.normal(size=(60, 3))
        campaign_error = rng.normal(size=(5, 3))
        spatial = np.stack([chol @ rng.normal(size=(60, 3)) for _ in range(5)])
        y = site_error[site_index] + campaign_error[trip_index] + spatial[trip_index, site_index] + rng.normal(size=(len(observations), 3))
        ry = pulse.residualize(y, basis)
        ry /= np.sqrt(np.square(ry).sum(axis=0))
        statistic = np.einsum("onp,ne->ope", candidate, ry).max(axis=(1, 2))
        assigned = int(rng.integers(len(candidate)))
        records.append({"replicate": replicate + 1, "assigned_year_permutation": assigned,
                        "p": exhaustive_p(statistic, assigned)})
    frame = pd.DataFrame(records)
    summaries = []
    for alpha in (.05, .2, .5):
        count = int((frame.p <= alpha + 1e-12).sum())
        ci = stats.binomtest(count, replicates).proportion_ci()
        summaries.append({"nominal_alpha": alpha, "rejections": count,
                          "simulations": replicates, "fraction": count / replicates,
                          "binomial_ci_low": ci.low, "binomial_ci_high": ci.high})
    return frame, summaries


def run(args):
    root, out = args.root.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    alpha = root / "analysis/v2/review/cache/alpha.tsv"
    observations = pulse.load_group_means(alpha, geodata(args.geodata))
    products = {"nasa_power": base.load_nasa(args.nasa), "open_meteo": base.load_open_meteo(args.open_meteo)}
    years = eligible_years(observations, list(products.values()))
    plan = {"status": "fixed_before_observed_tests", "eligible_years": years,
            "calendar": "December-November fields, exact month/day collection, 60 actual complete days",
            "null": "joint exchangeability of eligible whole-year rainfall fields, conditional on community responses and collection design",
            "coupling": "same year map for both products, all sites and campaigns; campaigns sharing a year remain coupled",
            "family": "both rainfall products, all 60 peaks, all 3 endpoints; one-sided positive and two-sided summaries",
            "bootstrap_draws": args.bootstraps, "seed": args.seed,
            "bootstrap": "whole-site multinomial reweighting; site-compartment and campaign-route-day nuisance coefficients refitted; full endpoint/peak search per product",
            "wet_spells": "maximal consecutive calendar days with positive rainfall at any of 60 sites, separately within each campaign lookback; omit whole regional spell"}
    (out / "analysis_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    design = RefitDesign(observations, observations[ENDPOINTS].to_numpy())
    orbit, exposures = calendar_orbit(observations, products, design, years)
    orbit.to_csv(out / "calendar_year_orbit.tsv", sep="\t", index=False)
    boot, summary, observed = bootstrap_refits(design, exposures, args.bootstraps, args.seed)
    boot.to_csv(out / "refitted_site_bootstrap.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    summary.to_csv(out / "refitted_site_bootstrap_summary.tsv", sep="\t", index=False)
    wet_spell_deletion(observations, products, design).to_csv(out / "wet_spell_deletion.tsv", sep="\t", index=False)
    calibration, calibration_summary = rank_calibration(observations, products, years, args.seed + 1, args.simulations)
    calibration.to_csv(out / "calendar_null_calibration.tsv", sep="\t", index=False)
    observations.to_csv(out / "analysis_cohort.tsv", sep="\t", index=False)
    result = {"observed": observed, "year_orbit_size": len(orbit),
              "exact_p": {c: exhaustive_p(orbit[c]) for c in orbit if c.endswith("_max")},
              "bootstrap_endpoint_frequencies": {p: g.selected_endpoint.value_counts(normalize=True).to_dict() for p, g in boot.groupby("product")},
              "calibration": calibration_summary,
              "calibration_scope": "conditional randomized-year null algorithm; exchangeability of the real annual fields is a scientific assumption",
              "interval_scope": "percentile ranges of raw whole-site refits with frozen rainfall; selected-peak ranges summarize reselection"}
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    inputs = [alpha, args.nasa, args.open_meteo, Path(__file__), Path(pulse.__file__), Path(base.__file__)]
    inputs += [args.geodata / f"trip{t}_geodata.tsv" for t in range(1, 6)]
    manifest = {"inputs": [{"path": str(p.resolve()), "sha256": digest(p)} for p in inputs],
                "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                             "scipy": scipy.__version__, "patsy": patsy.__version__},
                "outputs": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "manifest.json"}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--geodata", type=Path, required=True)
    parser.add_argument("--nasa", type=Path, required=True)
    parser.add_argument("--open-meteo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=9999)
    parser.add_argument("--simulations", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260910)
    run(parser.parse_args())
