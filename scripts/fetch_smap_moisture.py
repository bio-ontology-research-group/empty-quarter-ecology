"""Fetch SMAP L3 daily 9-km soil moisture for 60 Empty Quarter sites.

Uses earthaccess with ~/.netrc auth. Downloads daily granules for
the full study period (2022-01-01 to 2025-09-30) and extracts
per-site values.

SMAP L3 passive SM (SPL3SMP) short-name: SPL3SMP
See: https://nsidc.org/data/spl3smp
"""

from __future__ import annotations

import sys
from pathlib import Path

import earthaccess
import numpy as np
import pandas as pd
import h5py

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

GEODATA = REPO / "data" / "geodata"


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in (1, 2, 3, 4, 5):
        df = pd.read_csv(GEODATA / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"]).drop_duplicates("SiteNum")
        df = df[(df.SiteNum >= 1) & (df.SiteNum <= 60) & (df.SiteNum == df.SiteNum.astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (
        pd.concat(frames).groupby("SiteNum")[["Latitude", "Longitude"]].mean()
        .reset_index().sort_values("SiteNum").reset_index(drop=True)
    )


def latlon_to_ease2_grid_idx(lat: float, lon: float,
                              cell_deg: float = 0.09) -> tuple[int, int]:
    """Approximate EASE2 9km grid index; exact extraction uses file lat/lon arrays."""
    # placeholder — real extraction uses h5 lat/lon arrays
    return int((90 - lat) / cell_deg), int((lon + 180) / cell_deg)


def extract_from_file(path: Path, sites: pd.DataFrame) -> list[dict]:
    """Return list of {Site, Date, SM_am, SM_pm} dicts from one SPL3SMP h5."""
    rows = []
    with h5py.File(path, "r") as f:
        # SMAP SPL3SMP structure: /Soil_Moisture_Retrieval_Data_AM and _PM
        date_str = path.name.split("_")[4]   # SMAP_L3_SM_P_YYYYMMDD_...
        date = pd.to_datetime(date_str, format="%Y%m%d")
        for half, gkey in [("am", "Soil_Moisture_Retrieval_Data_AM"),
                            ("pm", "Soil_Moisture_Retrieval_Data_PM")]:
            try:
                grp = f[gkey]
                lat_arr = grp["latitude" if half == "am" else "latitude_pm"][:]
                lon_arr = grp["longitude" if half == "am" else "longitude_pm"][:]
                sm_arr = grp["soil_moisture" if half == "am"
                             else "soil_moisture_pm"][:]
            except KeyError:
                continue
            fill = -9999.0
            sm_arr = np.where(sm_arr <= fill + 1, np.nan, sm_arr)

            # For each site, find nearest valid cell
            for _, s in sites.iterrows():
                # Grid is regular EASE2 at 9 km; find nearest by cell
                # Build a small dist mask
                dlat = np.abs(lat_arr - s.Latitude)
                dlon = np.abs(lon_arr - s.Longitude)
                # find within 0.2° box
                mask = (dlat < 0.15) & (dlon < 0.2)
                if not mask.any():
                    continue
                vals = sm_arr[mask]
                good = vals[np.isfinite(vals)]
                if len(good) == 0:
                    continue
                rows.append({"Site": int(s.SiteNum), "Date": date,
                             "half": half, "SM": float(good.mean())})
    return rows


def main():
    sites = pooled_sites()
    print(f"auth earthdata ...", flush=True)
    auth = earthaccess.login(strategy="netrc")
    if not auth.authenticated:
        raise RuntimeError("earthdata auth failed — check ~/.netrc")
    print("  authenticated OK")

    # Empty Quarter bbox + 0.1° pad
    bbox = (sites.Longitude.min() - 0.1, sites.Latitude.min() - 0.1,
            sites.Longitude.max() + 0.1, sites.Latitude.max() + 0.1)
    print(f"  bbox = {bbox}")

    # SMAP L3 passive daily (9 km)
    short_name = "SPL3SMP"

    # Search in monthly windows to avoid huge result lists
    all_rows = []
    out_tmp = CACHE_DIR / "smap_tmp"
    out_tmp.mkdir(exist_ok=True)

    # Download in month chunks
    starts = pd.date_range("2022-01-01", "2025-09-01", freq="MS")
    for s_date in starts:
        e_date = (s_date + pd.offsets.MonthEnd(1))
        month_key = s_date.strftime("%Y%m")
        month_rows_file = CACHE_DIR / "smap_tmp" / f"{month_key}.parquet"
        if month_rows_file.exists():
            print(f"  {month_key} cached")
            df_m = pd.read_parquet(month_rows_file)
            all_rows.append(df_m)
            continue

        print(f"  searching {month_key} ({s_date.date()} - {e_date.date()})",
              flush=True)
        try:
            results = earthaccess.search_data(
                short_name=short_name,
                bounding_box=bbox,
                temporal=(s_date.strftime("%Y-%m-%d"),
                          e_date.strftime("%Y-%m-%d")),
            )
        except Exception as e:
            print(f"    search failed: {e}")
            continue
        print(f"    {len(results)} granules", flush=True)
        if not results:
            continue
        try:
            files = earthaccess.download(results, str(out_tmp))
        except Exception as e:
            print(f"    download failed: {e}"); continue
        month_rows = []
        for fp in files:
            fp = Path(fp)
            if not fp.exists() or fp.stat().st_size < 1000:
                continue
            try:
                month_rows.extend(extract_from_file(fp, sites))
            except Exception as e:
                print(f"      skip {fp.name}: {e}"); continue
            # delete file to save space
            try: fp.unlink()
            except Exception: pass
        df_m = pd.DataFrame(month_rows)
        df_m.to_parquet(month_rows_file, index=False)
        all_rows.append(df_m)
        print(f"    extracted {len(df_m)} site-day rows", flush=True)

    full = pd.concat(all_rows, ignore_index=True)
    # collapse AM/PM to daily mean
    daily = (full.groupby(["Site","Date"])["SM"].mean()
             .reset_index())
    out = CACHE_DIR / "smap_daily.parquet"
    daily.to_parquet(out, index=False)
    print(f"\nwrote {len(daily):,} site-day rows → {out}")
    print(daily.groupby("Site").SM.mean().describe().round(3))


if __name__ == "__main__":
    main()
