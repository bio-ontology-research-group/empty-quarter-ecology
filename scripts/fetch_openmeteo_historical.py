"""
Fetch ERA5 daily T_mean and precipitation for each of the 60 EQ sites
from 1980-01-01 through 2024-12-31 via the Open-Meteo archive API.

Used for the per-site historical climate-trend analysis (Hill-threshold
pulse-rainfall extremity, mean annual T trend).

Outputs:
  cache/climate_historical_1980_2024.parquet — site, date, T, P
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
from eq import CACHE_DIR  # noqa: E402

GEODATA = REPO / "data" / "geodata"
ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"

START = "1995-01-01"
END = "2024-12-31"
DAILY_VARS = ["temperature_2m_mean", "precipitation_sum"]
TIMEOUT = 120


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[
            (df.SiteNum >= 1) & (df.SiteNum <= 60)
            & (df.SiteNum == df.SiteNum.astype(int))
        ]
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


def fetch_site(lat: float, lon: float, site: int) -> pd.DataFrame | None:
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": START,
        "end_date": END,
        "daily": ",".join(DAILY_VARS),
        "timezone": "UTC",
    }
    for attempt in range(8):
        try:
            r = requests.get(ARCHIVE, params=params, timeout=TIMEOUT)
            if r.status_code == 200:
                j = r.json()
                d = j.get("daily")
                if d is None:
                    print(f"  site {site}: no daily payload")
                    return None
                df = pd.DataFrame(
                    {
                        "site": site,
                        "date": pd.to_datetime(d["time"]),
                        "T": d.get("temperature_2m_mean"),
                        "P": d.get("precipitation_sum"),
                    }
                )
                return df
            elif r.status_code == 429:
                wait = 30 * (attempt + 1)
                print(f"  site {site}: 429, sleep {wait}s")
                time.sleep(wait)
            else:
                print(f"  site {site}: HTTP {r.status_code}, attempt {attempt}")
                time.sleep(2 ** attempt)
        except (requests.RequestException, ValueError) as e:
            print(f"  site {site}: {e!r}, attempt {attempt}")
            time.sleep(2 ** attempt)
    return None


def main() -> int:
    sites = pooled_sites()
    print(f"pulling 1980–2024 daily T/P for {len(sites)} sites …")
    out_rows = []
    for _, r in sites.iterrows():
        site = int(r.SiteNum)
        lat, lon = float(r.Latitude), float(r.Longitude)
        df = fetch_site(lat, lon, site)
        if df is None or df.empty:
            print(f"  site {site}: FAILED")
            continue
        out_rows.append(df)
        print(
            f"  site {site}: {len(df):,} rows, "
            f"T mean {np.nanmean(df['T']):.1f}°C, "
            f"P sum {np.nansum(df['P']):.0f} mm",
            flush=True,
        )
        # Save partial output every 5 sites so we don't lose progress
        if len(out_rows) % 5 == 0:
            big = pd.concat(out_rows, ignore_index=True)
            big.to_parquet(CACHE_DIR / "climate_historical_1995_2024.parquet")
        time.sleep(8.0)  # be polite — Open-Meteo rate limits
    if not out_rows:
        print("ERROR: no data")
        return 1
    big = pd.concat(out_rows, ignore_index=True)
    out = CACHE_DIR / "climate_historical_1995_2024.parquet"
    big.to_parquet(out)
    print(f"\nwrote {out}  ({len(big):,} rows)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
