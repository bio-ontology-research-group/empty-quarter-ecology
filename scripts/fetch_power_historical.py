"""
Fetch NASA POWER daily T_mean and precipitation for the 60 EQ sites
from 1995-01-01 through 2024-12-31.

NASA POWER (Langley) is free, no key, no aggressive rate limits.
Data resolution is 0.5° (coarser than ERA5's 0.25° via Open-Meteo,
but more than fine enough for site-level climate trends across our
1000-km transect).

Outputs:
  cache/climate_historical_1995_2024.parquet — site, date, T, P
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
POWER = "https://power.larc.nasa.gov/api/temporal/daily/point"

START = "19950101"
END = "20241231"
TIMEOUT = 120


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60)]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (
        pd.concat(frames)
        .groupby("SiteNum")[["Latitude", "Longitude"]]
        .mean()
        .reset_index()
        .sort_values("SiteNum")
        .reset_index(drop=True)
    )


def fetch_site(lat: float, lon: float, site: int) -> pd.DataFrame | None:
    params = {
        "parameters": "T2M,PRECTOTCORR",
        "community": "AG",
        "longitude": lon,
        "latitude": lat,
        "start": START,
        "end": END,
        "format": "JSON",
    }
    for attempt in range(5):
        try:
            r = requests.get(POWER, params=params, timeout=TIMEOUT)
            if r.status_code == 200:
                j = r.json()
                params_data = j.get("properties", {}).get("parameter", {})
                t2m = params_data.get("T2M", {})
                pre = params_data.get("PRECTOTCORR", {})
                if not t2m:
                    print(f"  site {site}: no T2M data", flush=True)
                    return None
                dates = sorted(t2m.keys())
                df = pd.DataFrame(
                    {
                        "site": site,
                        "date": pd.to_datetime(dates, format="%Y%m%d"),
                        "T": [t2m[d] if t2m[d] != -999.0 else np.nan for d in dates],
                        "P": [pre.get(d, np.nan) if pre.get(d, np.nan) != -999.0
                              else np.nan for d in dates],
                    }
                )
                return df
            else:
                print(f"  site {site}: HTTP {r.status_code} attempt {attempt}",
                      flush=True)
                time.sleep(5 * (attempt + 1))
        except (requests.RequestException, ValueError) as e:
            print(f"  site {site}: {e!r} attempt {attempt}", flush=True)
            time.sleep(5 * (attempt + 1))
    return None


def main() -> int:
    sites = pooled_sites()
    print(f"pulling NASA POWER 1995–2024 daily T/P for {len(sites)} sites …",
          flush=True)
    out_rows: list[pd.DataFrame] = []
    for _, r in sites.iterrows():
        site = int(r.SiteNum)
        lat, lon = float(r.Latitude), float(r.Longitude)
        df = fetch_site(lat, lon, site)
        if df is None or df.empty:
            print(f"  site {site}: FAILED", flush=True)
            continue
        out_rows.append(df)
        print(
            f"  site {site}: {len(df):,} rows, "
            f"T mean {np.nanmean(df['T']):.1f}°C, "
            f"P sum {np.nansum(df['P']):.0f} mm",
            flush=True,
        )
        if len(out_rows) % 5 == 0:
            pd.concat(out_rows, ignore_index=True).to_parquet(
                CACHE_DIR / "climate_historical_1995_2024.parquet"
            )
        time.sleep(1.5)
    if not out_rows:
        print("ERROR: no data", flush=True)
        return 1
    big = pd.concat(out_rows, ignore_index=True)
    out = CACHE_DIR / "climate_historical_1995_2024.parquet"
    big.to_parquet(out)
    print(f"\nwrote {out}  ({len(big):,} rows; {big.site.nunique()} sites)",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
