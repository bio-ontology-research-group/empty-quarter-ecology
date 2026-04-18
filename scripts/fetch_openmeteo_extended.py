"""Fetch Open-Meteo ERA5 daily variables beyond temp/rain.

Pulls uv_index_max, shortwave_radiation_sum, et0_fao_evapotranspiration,
and hourly soil_moisture_0_to_7cm (aggregated to daily mean) for all
60 sampled sites. Saves to cache/climate_extended.parquet.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from eq import CACHE_DIR

GEODATA = REPO / "data" / "geodata"
ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"

DAILY_VARS = [
    "temperature_2m_mean",
    "temperature_2m_max",
    "temperature_2m_min",
    "precipitation_sum",
    "rain_sum",
    "shortwave_radiation_sum",
    "uv_index_max",
    "uv_index_clear_sky_max",
    "et0_fao_evapotranspiration",
    "vapour_pressure_deficit_max",
    "wind_speed_10m_max",
]
HOURLY_VARS = [
    "soil_moisture_0_to_7cm",
    "soil_temperature_0_to_7cm",
]

START = "2022-01-01"
END = "2025-09-30"
TIMEOUT = 60


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60) & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    pooled = (
        pd.concat(frames)
        .groupby("SiteNum")[["Latitude", "Longitude"]]
        .mean()
        .reset_index()
        .sort_values("SiteNum")
        .reset_index(drop=True)
    )
    return pooled


def fetch_site(lat: float, lon: float, site: int) -> pd.DataFrame:
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": START,
        "end_date": END,
        "daily": ",".join(DAILY_VARS),
        "hourly": ",".join(HOURLY_VARS),
        "timezone": "UTC",
    }
    for attempt in range(8):
        try:
            r = requests.get(ARCHIVE, params=params, timeout=TIMEOUT)
            if r.status_code == 429:
                retry = int(r.headers.get("Retry-After", 60))
                print(f"  site {site} 429; sleeping {retry}s", flush=True)
                time.sleep(retry)
                continue
            r.raise_for_status()
            j = r.json()
            break
        except Exception as e:
            wait = min(60, 5 * (attempt + 1))
            print(f"  site {site} attempt {attempt+1} failed ({e}); sleeping {wait}s",
                  flush=True)
            time.sleep(wait)
    else:
        raise RuntimeError(f"site {site} failed after 8 attempts")

    daily = pd.DataFrame({"Date": j["daily"]["time"]})
    for v in DAILY_VARS:
        daily[v] = j["daily"][v]

    hourly = pd.DataFrame({"Time": j["hourly"]["time"]})
    for v in HOURLY_VARS:
        hourly[v] = j["hourly"][v]
    hourly["Date"] = hourly["Time"].str[:10]
    hourly_daily = hourly.groupby("Date")[HOURLY_VARS].mean().reset_index()

    out = daily.merge(hourly_daily, on="Date", how="left")
    out.insert(0, "Site", site)
    return out


def main():
    sites = pooled_sites()
    tmp_dir = CACHE_DIR / "climate_extended_tmp"
    tmp_dir.mkdir(exist_ok=True)
    print(f"Fetching {len(sites)} sites × {len(DAILY_VARS) + len(HOURLY_VARS)} vars "
          f"from {START} to {END} (resume-capable)", flush=True)
    rows = []
    for _, r in sites.iterrows():
        s = int(r.SiteNum)
        ck = tmp_dir / f"site_{s:02d}.parquet"
        if ck.exists():
            print(f"  site {s:02d} cached", flush=True)
            rows.append(pd.read_parquet(ck))
            continue
        print(f"  site {s:02d} ({r.Latitude:.3f}, {r.Longitude:.3f})", flush=True)
        df = fetch_site(r.Latitude, r.Longitude, s)
        df.to_parquet(ck, index=False)
        rows.append(df)
        time.sleep(20.0)  # 429-safe with daily cap: conservative 3/min
    all_df = pd.concat(rows, ignore_index=True)
    out = CACHE_DIR / "climate_extended.parquet"
    all_df.to_parquet(out, index=False)
    print(f"wrote {len(all_df):,} rows → {out}", flush=True)
    summary = (
        all_df.groupby("Site")
        .agg(
            days=("Date", "count"),
            uv_max=("uv_index_max", "max"),
            uv_mean=("uv_index_max", "mean"),
            et0_mean=("et0_fao_evapotranspiration", "mean"),
            sm_mean=("soil_moisture_0_to_7cm", "mean"),
        )
        .round(3)
    )
    print(summary.head(), flush=True)
    print("  ... 60 sites total", flush=True)


if __name__ == "__main__":
    main()
