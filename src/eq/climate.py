"""Climate / weather joins: per-sample rainfall windows and site-level means.

The Open-Meteo-derived daily weather table in
``data/climate/daily_weather.tsv`` has one row per (site, date) with
temperature and precipitation. Every sample has a known trip date which
we join here to produce antecedent-rainfall summaries at multiple
windows for the distributed-lag rainfall-response analysis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import CACHE_DIR, DATA_DIR

DAILY_WEATHER_PATH = DATA_DIR / "climate" / "daily_weather.tsv"
GEODATA_DIR = DATA_DIR / "geodata"

TRIP_CENTER_DATE = {
    1: "2023-03-17",
    2: "2023-07-14",
    3: "2024-02-17",
    4: "2024-08-20",
    5: "2025-10-12",
}


def load_daily_weather() -> pd.DataFrame:
    w = pd.read_csv(DAILY_WEATHER_PATH, sep="\t", parse_dates=["Date"])
    w.columns = [c.strip() for c in w.columns]
    return w.rename(columns={"Site": "site", "Date": "date"})


def rainfall_windows(
    meta: pd.DataFrame, windows_days: tuple[int, ...] = (1, 3, 7, 10, 14, 30)
) -> pd.DataFrame:
    """For each sample, compute cumulative rainfall over each window ending
    on the trip's sampling date. Expects ``meta`` with ``trip`` and ``site``
    columns; adds ``rain_W{N}d``, ``temp_mean_W{N}d`` columns."""
    w = load_daily_weather()
    w["site"] = pd.to_numeric(w["site"], errors="coerce")
    w = w.dropna(subset=["site", "date"]).astype({"site": int})

    out = meta.copy()
    out["trip_date"] = out["trip"].map(
        lambda t: pd.to_datetime(TRIP_CENTER_DATE.get(int(t))) if pd.notna(t) else pd.NaT
    )

    # Per-site sub-frames keyed by date → fast positional lookup
    per_site = {s: g.sort_values("date").set_index("date")
                for s, g in w.groupby("site")}
    for win in windows_days:
        rain = []
        tmean = []
        for (_, r) in out.iterrows():
            site = r.get("site")
            end = r.get("trip_date")
            if pd.isna(site) or pd.isna(end) or int(site) not in per_site:
                rain.append(np.nan); tmean.append(np.nan); continue
            start = end - pd.Timedelta(days=win)
            sub = per_site[int(site)].loc[start:end]
            if sub.empty:
                rain.append(np.nan); tmean.append(np.nan); continue
            rain.append(float(sub["Rain_mm"].sum()))
            tmean.append(float(sub["Mean_Temp_C"].mean()))
        out[f"rain_W{win}d"] = rain
        out[f"temp_mean_W{win}d"] = tmean
    return out


def add_rainfall_to_metadata() -> pd.DataFrame:
    """Read cached metadata, add antecedent-rainfall windows, cache result."""
    meta = pd.read_parquet(CACHE_DIR / "metadata.parquet").reset_index()
    out = rainfall_windows(meta.dropna(subset=["trip", "site"]))
    out.to_parquet(CACHE_DIR / "metadata_with_rainfall.parquet")
    return out
