#!/usr/bin/env python3
"""Build pairwise inter-site geometry: haversine distance and bearing
(azimuth from i to j, degrees from north, meteorological convention).

Inputs:
  data/geodata/trip{1..5}_geodata.tsv   (60 sites, lat/lon)

Outputs:
  cache/pairwise_geometry.tsv  columns: site_i, site_j, dist_km, bearing_ij_deg
                              (also bearing_ji_deg = (bearing_ij_deg+180)%360)
"""
from __future__ import annotations

from pathlib import Path
import math
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
GEO = REPO / "data" / "geodata"
OUT = REPO / "cache" / "pairwise_geometry.tsv"


def pooled_sites() -> pd.DataFrame:
    frames = []
    for t in range(1, 6):
        df = pd.read_csv(GEO / f"trip{t}_geodata.tsv", sep="\t")
        df["SiteNum"] = pd.to_numeric(df["Site"], errors="coerce")
        df = df.dropna(subset=["SiteNum"])
        df = df[(df["SiteNum"] >= 1) & (df["SiteNum"] <= 60)
                & (df["SiteNum"] == df["SiteNum"].astype(int))]
        frames.append(df[["SiteNum", "Latitude", "Longitude"]])
    return (pd.concat(frames)
            .groupby("SiteNum")[["Latitude", "Longitude"]].mean()
            .reset_index().sort_values("SiteNum").reset_index(drop=True))


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    p1 = math.radians(lat1); p2 = math.radians(lat2)
    dp = math.radians(lat2 - lat1); dl = math.radians(lon2 - lon1)
    a = math.sin(dp/2)**2 + math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2 * R * math.asin(math.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    """Initial bearing from (lat1,lon1) to (lat2,lon2). Met convention:
    0=N, 90=E, 180=S, 270=W."""
    p1 = math.radians(lat1); p2 = math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    x = math.sin(dl) * math.cos(p2)
    y = math.cos(p1)*math.sin(p2) - math.sin(p1)*math.cos(p2)*math.cos(dl)
    return (math.degrees(math.atan2(x, y)) + 360) % 360


def main():
    sites = pooled_sites()
    print(f"sites: {len(sites)}", flush=True)
    rows = []
    for i in range(len(sites)):
        si, lai, loi = int(sites.SiteNum[i]), sites.Latitude[i], sites.Longitude[i]
        for j in range(i + 1, len(sites)):
            sj, laj, loj = int(sites.SiteNum[j]), sites.Latitude[j], sites.Longitude[j]
            d = haversine_km(lai, loi, laj, loj)
            bij = bearing_deg(lai, loi, laj, loj)
            bji = (bij + 180) % 360
            rows.append({"site_i": si, "site_j": sj, "dist_km": round(d, 3),
                         "bearing_ij_deg": round(bij, 2),
                         "bearing_ji_deg": round(bji, 2)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT.name}: {len(df)} pairs", flush=True)
    print(f"  distance range: {df.dist_km.min():.1f} – {df.dist_km.max():.1f} km")
    print(f"  median: {df.dist_km.median():.1f} km")
    print(df.head(3).to_string(index=False))


if __name__ == "__main__":
    main()
