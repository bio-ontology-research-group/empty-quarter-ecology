"""Fetch NASA POWER daily UV/radiation/soil-moisture variables for 60 sites.

POWER API (NASA/Langley) is free, no key, no aggressive rate limits.
Variables:
- ALLSKY_SFC_UV_INDEX         daily max UV index
- ALLSKY_SFC_SW_DWN           daily shortwave radiation (MJ/m²/d)
- EVPTRNS                     evapotranspiration (mm/d)
- GWETTOP                     0-10 cm top-layer soil moisture (fraction)
- TS                          earth skin temperature (°C)
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
POWER = "https://power.larc.nasa.gov/api/temporal/daily/point"

VARS = [
    "ALLSKY_SFC_UV_INDEX",
    "ALLSKY_SFC_SW_DWN",
    "EVPTRNS",
    "GWETTOP",
    "TS",
    "PS",  # surface pressure
    "WS2M",
]
START = "20220101"
END = "20250930"


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60) & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (
        pd.concat(frames)
        .groupby("SiteNum")[["Latitude", "Longitude"]]
        .mean()
        .reset_index()
        .sort_values("SiteNum")
        .reset_index(drop=True)
    )


def fetch_site(lat: float, lon: float, site: int) -> pd.DataFrame:
    params = {
        "parameters": ",".join(VARS),
        "community": "AG",
        "longitude": lon, "latitude": lat,
        "start": START, "end": END,
        "format": "JSON",
    }
    for attempt in range(5):
        try:
            r = requests.get(POWER, params=params, timeout=60)
            r.raise_for_status()
            j = r.json()
            data = j["properties"]["parameter"]
            df = pd.DataFrame({v: pd.Series(data[v]) for v in VARS})
            df.index.name = "YYYYMMDD"
            df = df.reset_index()
            df["Date"] = pd.to_datetime(df["YYYYMMDD"], format="%Y%m%d")
            df.insert(0, "Site", site)
            return df.drop(columns="YYYYMMDD")
        except Exception as e:
            wait = 5 * (attempt + 1)
            print(f"  site {site} attempt {attempt+1} failed ({e}); sleep {wait}s",
                  flush=True)
            time.sleep(wait)
    raise RuntimeError(f"site {site} failed after 5 attempts")


def main():
    sites = pooled_sites()
    tmp_dir = CACHE_DIR / "nasa_power_tmp"
    tmp_dir.mkdir(exist_ok=True)
    print(f"Fetching {len(sites)} sites × {len(VARS)} NASA-POWER daily vars", flush=True)
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
        time.sleep(1.5)
    all_df = pd.concat(rows, ignore_index=True)
    out = CACHE_DIR / "nasa_power_daily.parquet"
    # replace NASA POWER sentinel -999 with NaN
    for v in VARS:
        all_df.loc[all_df[v] <= -900, v] = np.nan
    all_df.to_parquet(out, index=False)
    print(f"wrote {len(all_df):,} rows → {out}", flush=True)
    print(all_df.groupby("Site")[VARS].mean().head().round(2))


if __name__ == "__main__":
    main()
