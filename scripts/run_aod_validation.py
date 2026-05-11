#!/usr/bin/env python3
"""Validate wind-Mantel results against MERRA-2 AOD (dust optical depth).

Two analyses:
  A. Threshold validation: per (site, month), does WS10M_MAX > 7 m/s day count
     correlate with monthly DUEXTTAU? If yes -> our threshold-based dust-uplift
     proxy is meaningful.

  B. AOD interaction: stratify pair-trip Mantel results by trip-window dust
     intensity. In high-DUEXTTAU periods, does wind connectivity explain more
     of community dissimilarity than in low-DUEXTTAU periods?

Inputs:
  cache/merra2_aod/per_site_monthly.csv     (60 sites x 48 months)
  data/climate/daily_weather_full.csv        (NASA POWER daily wind)
  cache/wind_dispersal/sweep_mantel_full.tsv (existing sweep)

Outputs:
  cache/merra2_aod/threshold_validation.tsv
  cache/merra2_aod/site_dust_threshold_correlation.tsv
  cache/merra2_aod/aod_x_wind_interaction.tsv
  cache/merra2_aod/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "merra2_aod"


def main():
    aod = pd.read_csv(OUT / "per_site_monthly.csv")
    aod["date"] = pd.to_datetime(aod["date"])
    aod["month"] = aod["date"].dt.to_period("M")
    print(f"AOD rows: {len(aod)}, sites: {aod['site'].nunique()}, "
          f"months: {aod['month'].nunique()}", flush=True)

    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather["site"] = weather["site"].astype(int)
    weather["month"] = weather["date"].dt.to_period("M")

    # A. Threshold validation: WS10M_MAX>thr days vs DUEXTTAU per (site, month)
    site_month = []
    for thr in [5.0, 7.0, 10.0]:
        wm = (weather.assign(dust_day=(weather["WS10M_MAX"] > thr).astype(int))
                     .groupby(["site", "month"])
                     .agg(dust_days=("dust_day", "sum"),
                          mean_WS10M_MAX=("WS10M_MAX", "mean"))
                     .reset_index())
        wm["thr"] = thr
        site_month.append(wm)
    sm = pd.concat(site_month)

    merged = sm.merge(aod[["site", "month", "DUEXTTAU", "TOTEXTTAU"]],
                      on=["site", "month"], how="inner")
    print(f"\n[validation] Wind dust-day count vs MERRA-2 DUEXTTAU "
          f"(per site-month, n={len(merged)//3} pairs):")
    rows = []
    for thr, grp in merged.groupby("thr"):
        rho_dust, p_dust = spearmanr(grp["dust_days"], grp["DUEXTTAU"])
        rho_tot,  p_tot  = spearmanr(grp["dust_days"], grp["TOTEXTTAU"])
        rho_mean, p_mean = spearmanr(grp["mean_WS10M_MAX"], grp["DUEXTTAU"])
        rows.append({"threshold_ms": thr,
                     "rho_dust_days_vs_DUEXTTAU": rho_dust, "p_dust_vs_DU": p_dust,
                     "rho_dust_days_vs_TOTEXTTAU": rho_tot, "p_dust_vs_TOT": p_tot,
                     "rho_meanWS_vs_DUEXTTAU": rho_mean, "p_mean_vs_DU": p_mean,
                     "n": len(grp)})
        print(f"  thr={thr:>4.1f}: rho(dust_days, DUEXTTAU)={rho_dust:+.3f} "
              f"(p={p_dust:.2g}); rho(dust_days, TOTEXTTAU)={rho_tot:+.3f}",
              flush=True)
    pd.DataFrame(rows).to_csv(OUT / "threshold_validation.tsv", sep="\t", index=False)

    # Per-site correlation
    sit_rows = []
    for site, grp in merged[merged["thr"] == 7.0].groupby("site"):
        if len(grp) < 6:
            continue
        rho, p = spearmanr(grp["dust_days"], grp["DUEXTTAU"])
        sit_rows.append({"site": int(site), "n_months": len(grp),
                         "rho": rho, "p": p,
                         "mean_DUEXTTAU": float(grp["DUEXTTAU"].mean()),
                         "mean_dust_days": float(grp["dust_days"].mean())})
    sit_df = pd.DataFrame(sit_rows)
    sit_df.to_csv(OUT / "site_dust_threshold_correlation.tsv",
                  sep="\t", index=False)

    print(f"\n[per-site] median rho={sit_df['rho'].median():+.3f}, "
          f"frac_positive={(sit_df['rho']>0).mean():.0%}, "
          f"frac_p<0.05={(sit_df['p']<0.05).mean():.0%}", flush=True)

    # B. AOD-stratified Mantel interaction:
    # For each trip's CenterDate, compute trip-window mean DUEXTTAU per site,
    # then stratify Mantel results by high-vs-low dust trip-windows.
    print("\n[interaction] not yet wired -- requires trip-resolved AOD windows",
          flush=True)

    # Summary
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("MERRA-2 AOD validation of wind-based dust proxy\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Per-(site, month) records: {len(merged)//3}\n\n")
        fh.write("Spearman correlation: dust-day count vs DUEXTTAU:\n")
        fh.write(pd.DataFrame(rows).to_string(index=False))
        fh.write("\n\nPer-site correlations (thr=7 m/s):\n")
        fh.write(f"  median rho: {sit_df['rho'].median():+.3f}\n")
        fh.write(f"  fraction with rho>0: {(sit_df['rho']>0).mean():.0%}\n")
        fh.write(f"  fraction with p<0.05: {(sit_df['p']<0.05).mean():.0%}\n")
        fh.write(f"  range mean_DUEXTTAU: "
                 f"{sit_df['mean_DUEXTTAU'].min():.3f} – "
                 f"{sit_df['mean_DUEXTTAU'].max():.3f}\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
