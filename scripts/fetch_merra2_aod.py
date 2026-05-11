"""Fetch MERRA-2 monthly aerosol optical depth (M2TMNXAER) for the EQ region
2022-2025 via earthaccess + Earthdata bearer-token auth in ~/.netrc.

Variables extracted (per-site monthly mean):
  TOTEXTTAU  Total aerosol optical thickness at 550 nm
  TOTSCATAU  Total aerosol scattering optical thickness
  DUEXTTAU   Dust extinction optical thickness
  DUSCATAU   Dust scattering optical thickness
  DUFLUXU/V  Dust column flux U/V (zonal/meridional, kg/m/s)
  SSEXTTAU   Sea-salt extinction (sanity baseline; should be small in EQ)

Output:
  cache/merra2_aod/raw/MERRA2_400.tavgM_2d_aer_Nx.YYYYMM.nc4   raw downloads
  cache/merra2_aod/per_site_monthly.csv                          60 sites x 48 months
  cache/merra2_aod/region_monthly.csv                            grid-region averages
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import earthaccess

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

OUT = CACHE_DIR / "merra2_aod"
RAW = OUT / "raw"
RAW.mkdir(parents=True, exist_ok=True)

VARS = ["TOTEXTTAU", "TOTSCATAU", "DUEXTTAU", "DUSCATAU",
        "DUFLUXU", "DUFLUXV", "SSEXTTAU"]
EQ_BBOX = (43.0, 18.0, 56.0, 23.0)  # (lon_min, lat_min, lon_max, lat_max)


def pooled_sites() -> pd.DataFrame:
    GEODATA = REPO / "data" / "geodata"
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60)
                & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (pd.concat(frames).groupby("SiteNum")[["Latitude", "Longitude"]].mean()
            .reset_index().sort_values("SiteNum").reset_index(drop=True))


def main():
    import time
    print("[auth] Earthdata via netrc", flush=True)
    auth = None
    last_err = None
    for attempt in range(6):
        try:
            auth = earthaccess.login(strategy="netrc")
            if auth.authenticated:
                break
        except Exception as e:
            last_err = e
            print(f"  attempt {attempt+1} failed: {e}; sleep 15s", flush=True)
            time.sleep(15)
    if auth is None or not auth.authenticated:
        raise RuntimeError(f"Earthdata auth failed after retries: {last_err}")

    print("[search] M2TMNXAER 2022-2025...", flush=True)
    granules = earthaccess.search_data(
        short_name="M2TMNXAER",
        version="5.12.4",
        bounding_box=EQ_BBOX,
        temporal=("2022-01-01", "2025-12-31"),
    )
    print(f"  {len(granules)} granules", flush=True)

    # Download (skips already-downloaded)
    print("[download] (skips existing)...", flush=True)
    files = earthaccess.download(granules, str(RAW))
    print(f"  {len(files)} files in {RAW}", flush=True)

    # Read each, slice to EQ region, average across region
    sites = pooled_sites()
    print(f"  sites: {len(sites)}", flush=True)
    site_rows = []
    region_rows = []
    for f in sorted(files):
        f = str(f)
        try:
            ds = xr.open_dataset(f, engine="h5netcdf")
        except Exception as e:
            print(f"  ! could not open {f}: {e}", flush=True)
            continue
        # Subset region
        ds_eq = ds.sel(lat=slice(EQ_BBOX[1], EQ_BBOX[3]),
                        lon=slice(EQ_BBOX[0], EQ_BBOX[2]))
        # Time is monthly: one timestep per file
        ts = pd.to_datetime(ds["time"].values[0])
        # Per-site nearest-grid-cell extraction
        for _, r in sites.iterrows():
            sub = ds_eq.sel(lat=float(r["Latitude"]),
                             lon=float(r["Longitude"]),
                             method="nearest")
            row = {"site": int(r["SiteNum"]), "date": ts}
            for v in VARS:
                if v in sub:
                    row[v] = float(sub[v].values[0])
            site_rows.append(row)
        # Regional mean
        rrec = {"date": ts}
        for v in VARS:
            if v in ds_eq:
                rrec[v] = float(ds_eq[v].mean().values)
        region_rows.append(rrec)
        ds.close()

    site_df = pd.DataFrame(site_rows).sort_values(["site", "date"])
    region_df = pd.DataFrame(region_rows).sort_values("date")
    site_df.to_csv(OUT / "per_site_monthly.csv", index=False)
    region_df.to_csv(OUT / "region_monthly.csv", index=False)
    print(f"\nWrote per-site (n={len(site_df)}) and region (n={len(region_df)}) "
          f"AOD time series.", flush=True)
    print("\nRegion summary (mean / max DUEXTTAU per year):")
    region_df["year"] = region_df["date"].dt.year
    print(region_df.groupby("year")[["DUEXTTAU", "TOTEXTTAU"]].agg(["mean", "max"]).to_string())


if __name__ == "__main__":
    main()
