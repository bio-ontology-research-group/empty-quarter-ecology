"""Fetch Copernicus GLO-30 DEM for Empty Quarter sites — 30 m resolution,
free, no authentication required, hosted as cloud-optimised GeoTIFFs
on AWS Open Data.

Extracts elevation + slope + aspect + topographic wetness index at
each of the 60 pooled-site coordinates.

Source: https://registry.opendata.aws/copernicus-dem/
Tiles: s3://copernicus-dem-30m/Copernicus_DSM_COG_10_N{NN}_00_E{EEE}_00_DEM/
"""

from __future__ import annotations

import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from eq import CACHE_DIR

GEODATA = REPO / "data" / "geodata"
AWS_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"


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


def tile_name(lat: float, lon: float) -> str:
    """Copernicus DEM tiling: 1° × 1° tiles named by SW corner."""
    lat_floor = int(np.floor(lat))
    lon_floor = int(np.floor(lon))
    ns = "N" if lat_floor >= 0 else "S"
    ew = "E" if lon_floor >= 0 else "W"
    return (f"Copernicus_DSM_COG_10_{ns}{abs(lat_floor):02d}_00_"
            f"{ew}{abs(lon_floor):03d}_00_DEM")


def download_tile(name: str, out_dir: Path) -> Path | None:
    out = out_dir / f"{name}.tif"
    if out.exists():
        return out
    url = f"{AWS_BASE}/{name}/{name}.tif"
    r = requests.get(url, timeout=300)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    out.write_bytes(r.content)
    return out


def sample_tile(tif: Path, lat: float, lon: float) -> dict:
    import rasterio
    with rasterio.open(tif) as src:
        elev = src.read(1).astype(float)
        elev[elev < -1000] = np.nan
        transform = src.transform
    # 1° tile: ~3600 × 3600 px for 30 m at equator; at ~20°N, 30 m = 30/(111320*cos(20))
    xres_m = abs(transform.a) * 111_320 * np.cos(np.deg2rad(lat))
    yres_m = abs(transform.e) * 111_320
    gy, gx = np.gradient(elev, yres_m, xres_m)
    slope = np.degrees(np.arctan(np.sqrt(gx**2 + gy**2)))
    aspect = (np.degrees(np.arctan2(-gx, gy)) + 360) % 360
    tan_slope = np.tan(np.deg2rad(slope)) + 1e-3
    twi = np.log(xres_m / tan_slope)
    import rasterio
    with rasterio.open(tif) as src:
        row, col = src.index(lon, lat)
    r = max(0, min(elev.shape[0] - 1, row))
    c = max(0, min(elev.shape[1] - 1, col))
    return {
        "elev_m": float(elev[r, c]),
        "slope_deg": float(slope[r, c]),
        "aspect_deg": float(aspect[r, c]),
        "twi": float(twi[r, c]),
    }


def main():
    sites = pooled_sites()
    cache_dir = CACHE_DIR / "cop_dem_tiles"
    cache_dir.mkdir(exist_ok=True)

    # Which tiles do we need?
    tiles = {tile_name(s.Latitude, s.Longitude) for _, s in sites.iterrows()}
    print(f"need {len(tiles)} tiles: {sorted(tiles)}")

    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(download_tile, t, cache_dir): t for t in tiles}
        tile_paths = {}
        for f in as_completed(futs):
            t = futs[f]
            try:
                p = f.result()
                if p is None:
                    print(f"  {t}: NOT FOUND (ocean/outside coverage)")
                else:
                    tile_paths[t] = p
                    print(f"  {t}: cached")
            except Exception as e:
                print(f"  {t}: FAIL {e}")

    print(f"\nextracting at {len(sites)} sites...")
    rows = []
    for _, s in sites.iterrows():
        tn = tile_name(s.Latitude, s.Longitude)
        if tn not in tile_paths:
            rows.append({"site": int(s.SiteNum), "elev_m": np.nan,
                         "slope_deg": np.nan, "aspect_deg": np.nan,
                         "twi": np.nan})
            continue
        vals = sample_tile(tile_paths[tn], s.Latitude, s.Longitude)
        vals["site"] = int(s.SiteNum)
        rows.append(vals)

    topo = pd.DataFrame(rows)
    out = CACHE_DIR / "topography.tsv"
    topo.to_csv(out, sep="\t", index=False)
    print(f"wrote → {out}")
    print(topo.describe().round(2))


if __name__ == "__main__":
    main()
