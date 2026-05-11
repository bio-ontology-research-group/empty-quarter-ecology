#!/usr/bin/env python3
"""Geographic / climate / chemistry correlates of per-site relic load.
Identifies environmental drivers of relic accumulation.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "relic_population"


def main():
    rec = pd.read_csv(OUT / "per_sample_relic_fraction.tsv", sep="\t")

    # Per-site median relic_frac
    per_site = (rec.groupby("site")
                .agg(median_relic_frac=("relic_frac", "median"),
                      mean_relic_frac=("relic_frac", "mean"),
                      n_samples=("sample", "count"))
                .reset_index())

    # Geodata (trip 1 baseline)
    geo = pd.read_csv(REPO / "data" / "geodata" / "trip1_geodata.tsv",
                       sep="\t")
    geo = geo.rename(columns={"Site": "site"})
    geo["site"] = pd.to_numeric(geo["site"], errors="coerce")
    geo = geo.dropna(subset=["site"])
    geo["site"] = geo["site"].astype(int)
    sites = per_site.merge(geo[["site", "Latitude", "Longitude",
                                   "AnnualMeanTemp", "AnnualTotalPrecip"]],
                              on="site", how="left")

    # XRF
    xrf = pd.read_csv(REPO / "data" / "geochemistry" /
                       "xrf_lab_table_all_trips.tsv", sep="\t")
    print(f"XRF columns: {xrf.columns.tolist()[:10]}", flush=True)
    if "site" in xrf.columns:
        xrf["site"] = pd.to_numeric(xrf["site"], errors="coerce")
        xrf = xrf.dropna(subset=["site"])
        xrf["site"] = xrf["site"].astype(int)
    xrf_num = xrf.select_dtypes(include="number").drop(columns=["site"],
                                                            errors="ignore")
    xrf_num["site"] = xrf["site"].values
    xrf_per_site = (xrf_num.groupby("site").median().reset_index()
                     if "site" in xrf_num.columns else None)
    if xrf_per_site is not None:
        sites = sites.merge(xrf_per_site, on="site", how="left",
                              suffixes=("", "_xrf"))

    sites.to_csv(OUT / "per_site_relic_with_env.tsv", sep="\t", index=False)
    print(f"\nSites with env data: {sites.dropna(subset=['Latitude']).shape[0]}",
          flush=True)

    # Correlations with relic_frac
    target = "median_relic_frac"
    print(f"\n=== Spearman correlations with {target} ===")
    rows = []
    for col in sites.select_dtypes(include="number").columns:
        if col in (target, "mean_relic_frac", "n_samples", "site"):
            continue
        s = sites[[target, col]].dropna()
        if len(s) < 10: continue
        r, p = spearmanr(s[target], s[col])
        rows.append({"variable": col, "n": len(s), "rho": float(r),
                      "p": float(p)})
    corr = pd.DataFrame(rows).sort_values("rho", key=lambda x: x.abs(),
                                              ascending=False)
    print(corr.head(20).round(4).to_string(index=False))
    corr.to_csv(OUT / "per_site_relic_env_corr.tsv", sep="\t", index=False)

    # Latitude/longitude gradient
    if "Latitude" in sites.columns:
        sites_geo = sites.dropna(subset=["Latitude", target])
        print(f"\n=== Spatial / climate ===")
        for col in ("Latitude", "Longitude", "AnnualMeanTemp",
                     "AnnualTotalPrecip"):
            if col in sites_geo.columns:
                r, p = spearmanr(sites_geo[target], sites_geo[col])
                print(f"  {col}: rho={r:+.3f}  p={p:.3g}  n={len(sites_geo)}")


if __name__ == "__main__":
    main()
