#!/usr/bin/env python3
"""Tier-1 #4: Pulse-reserve hypothesis at microbial scale (trip-resolved).

For each sample, look up its (site, trip) -> sample_date via geodata,
then summarise precip windows (7 / 30 / 90 days) prior to that date
from NASA POWER. Classify samples as `pulse` (recent precip) vs
`reserve` (none). For each abundant ASV, Mann-Whitney U test pulse vs
reserve abundance, and aggregate to genus.

Inputs:
  cache/feature_table.parquet
  cache/taxonomy.parquet
  data/climate/daily_weather_full.csv     (full NASA POWER, 19 vars)
  data/geodata/trip{1..5}_geodata.tsv     (CenterDate per site/trip)

Outputs:
  cache/pulse_reserve/sample_pulse_state.tsv
  cache/pulse_reserve/asv_pulse_response.tsv
  cache/pulse_reserve/genus_pulse_response.tsv
  cache/pulse_reserve/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "pulse_reserve"
OUT.mkdir(parents=True, exist_ok=True)

WINDOWS = {"d7": 7, "d30": 30, "d90": 90}
PULSE_MM_THRESHOLD = {"d7": 1.0, "d30": 5.0, "d90": 10.0}
PREV_MIN = 0.05


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet")
    print(f"feature_table: {ft.shape}", flush=True)

    # Build sample metadata via trip parser
    smeta = parse_samples_to_df(ft.columns)
    print(f"parsed: {len(smeta)} (per-trip: {smeta['trip'].value_counts().sort_index().to_dict()})",
          flush=True)

    # Build (site, trip) -> sample_date map from geodata
    geo_frames = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t")
            g["trip"] = trip
            g["CenterDate"] = pd.to_datetime(g["CenterDate"])
            g["Site_int"] = pd.to_numeric(g["Site"], errors="coerce")
            g = g.dropna(subset=["Site_int", "CenterDate"])
            g["site"] = g["Site_int"].astype(int)
            geo_frames.append(g[["site", "trip", "CenterDate"]])
    geo = pd.concat(geo_frames).drop_duplicates(["site", "trip"])
    smeta = smeta.merge(geo, on=["site", "trip"], how="left")
    smeta = smeta.dropna(subset=["CenterDate"])
    smeta["sample_date"] = smeta["CenterDate"]
    print(f"after date join: {len(smeta)} samples", flush=True)

    # Load full weather (with wind)
    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather["site"] = weather["site"].astype(int)
    print(f"weather: {weather.shape}, sites covered: {weather['site'].nunique()}",
          flush=True)

    # Pivot weather to (site, date) -> Precip for fast slicing
    pmap = weather[["site", "date", "PRECTOTCORR"]].rename(
        columns={"PRECTOTCORR": "precip_mm"})

    # For each sample compute trailing-window precip
    rec = []
    pmap_indexed = pmap.set_index(["site", "date"]).sort_index()
    for _, r in smeta.iterrows():
        site = int(r["site"])
        d = pd.Timestamp(r["sample_date"])
        row = {"sample": r["sample"], "site": site, "trip": int(r["trip"]),
               "compartment": r["compartment"], "sample_date": d}
        try:
            site_w = pmap_indexed.loc[site]
        except KeyError:
            for k in WINDOWS:
                row[f"precip_{k}"] = np.nan
            rec.append(row)
            continue
        for k, win in WINDOWS.items():
            mask = (site_w.index >= d - pd.Timedelta(days=win - 1)) & \
                   (site_w.index <= d)
            row[f"precip_{k}"] = float(site_w.loc[mask, "precip_mm"].sum())
        rec.append(row)
    pulse_df = pd.DataFrame(rec)
    pulse_df["pulse_d30"] = pulse_df["precip_d30"] >= PULSE_MM_THRESHOLD["d30"]
    pulse_df.to_csv(OUT / "sample_pulse_state.tsv", sep="\t", index=False)

    print("\nPulse fraction by compartment x trip:")
    print(pulse_df.groupby(["compartment", "trip"])["pulse_d30"].agg(
        ["count", "sum", "mean"]).to_string())
    print("\nOverall pulse fraction by compartment:")
    print(pulse_df.groupby("compartment")["pulse_d30"].mean().to_string())

    # Per-ASV Mann-Whitney
    rel = ft.div(ft.sum(axis=0).replace(0, 1), axis=1)
    prev = (ft > 0).sum(axis=1) / ft.shape[1]
    keep = prev[prev >= PREV_MIN].index
    print(f"\nASVs with prevalence >= {PREV_MIN}: {len(keep)}", flush=True)
    rel_k = rel.loc[keep]

    pulse_samples = pulse_df.loc[pulse_df["pulse_d30"], "sample"].tolist()
    reserve_samples = pulse_df.loc[~pulse_df["pulse_d30"], "sample"].tolist()
    print(f"pulse n={len(pulse_samples)}, reserve n={len(reserve_samples)}",
          flush=True)
    if min(len(pulse_samples), len(reserve_samples)) < 10:
        print("WARNING: small group size; thresholds may need tuning", flush=True)

    rows = []
    for asv in keep:
        x = rel_k.loc[asv]
        p = x.reindex(pulse_samples).dropna().values
        r = x.reindex(reserve_samples).dropna().values
        if len(p) < 5 or len(r) < 5:
            continue
        try:
            U, pv = mannwhitneyu(p, r, alternative="two-sided")
        except ValueError:
            continue
        mp = float(np.mean(p) + 1e-9)
        mr = float(np.mean(r) + 1e-9)
        rows.append({"ASV": asv, "n_pulse": len(p), "n_reserve": len(r),
                     "mean_pulse": mp, "mean_reserve": mr,
                     "log2_FC": float(np.log2(mp / mr)),
                     "U": float(U), "p_value": float(pv)})
    asv_df = pd.DataFrame(rows)
    if len(asv_df):
        asv_df["p_BH"] = multipletests(asv_df["p_value"].fillna(1.0),
                                        method="fdr_bh")[1]
        asv_df = asv_df.sort_values("p_value")
    asv_df.to_csv(OUT / "asv_pulse_response.tsv", sep="\t", index=False)

    # Aggregate to genus
    asv_to_genus = tax["genus"].fillna("Unclassified")
    asv_df["genus"] = asv_df["ASV"].map(asv_to_genus)
    g = (asv_df.groupby("genus")
              .agg(n_asv=("ASV", "count"),
                   median_log2FC=("log2_FC", "median"),
                   n_signif=("p_BH", lambda s: int((s < 0.05).sum())))
              .reset_index()
              .sort_values("median_log2FC", ascending=False))
    g.to_csv(OUT / "genus_pulse_response.tsv", sep="\t", index=False)

    print(f"\nTop 10 pulse-responding genera:")
    print(g.head(10).to_string(index=False))
    print("\nTop 10 reserve-favoured genera:")
    print(g.tail(10).to_string(index=False))

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Pulse-reserve analysis (Tier-1 #4)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Trip prefix mapping in use: 1=(none) 2=T 3=F 4=S 5=V\n\n")
        fh.write(f"Samples (with date join): {len(pulse_df)}\n")
        fh.write(f"  pulse  (precip30d >= {PULSE_MM_THRESHOLD['d30']} mm): "
                 f"{len(pulse_samples)}\n")
        fh.write(f"  reserve: {len(reserve_samples)}\n\n")
        fh.write("Per-compartment pulse fraction:\n")
        fh.write(pulse_df.groupby("compartment")["pulse_d30"].mean().to_string())
        fh.write(f"\n\nASVs tested (prev>={PREV_MIN}): {len(asv_df)}\n")
        if len(asv_df):
            fh.write(f"  significant (p_BH<0.05): {int((asv_df['p_BH']<0.05).sum())}\n")
        fh.write("\nTop 15 pulse-responsive genera:\n")
        fh.write(g.head(15).to_string(index=False))
        fh.write("\n\nTop 15 reserve-favoured genera:\n")
        fh.write(g.tail(15).to_string(index=False))
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
