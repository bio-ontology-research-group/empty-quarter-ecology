#!/usr/bin/env python3
"""Temporal dynamics of the two ecological strategies (DOM-cycling vs
halophile-spore) across the 5 trips.

Strategy A (DOM): Bact_DOM guild + Massilia + Lysobacter + Acidibacter
Strategy B (Halo): Bacilli halo + Halomonas

For each (site, compartment, trip):
  - Strategy A relabund (sum of A members)
  - Strategy B relabund (sum of B members)
  - Strategy ratio log2(A/B)

Then:
  1. Per-site temporal trajectories
  2. Strategy stability vs switching sites
  3. Correlation with precipitation lag (NASA POWER)
  4. Compartment-specific dynamics
  5. Trip-mean comparison

Outputs in cache/two_strategy_temporal/.
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
OUT = CACHE / "two_strategy_temporal"
OUT.mkdir(parents=True, exist_ok=True)


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


# Strategy membership (data-driven from earlier analyses)
STRATEGY_A = ["Nibribacter", "Daejeonella", "Cytophaga", "Cnuella",
                  "Niastella", "Flavitalea", "Ohtaekwangia", "Tellurirhabdus",
                  "Massilia", "Lysobacter", "Acidibacter"]
STRATEGY_B = ["Aquibacillus", "Sediminibacillus", "Litchfieldia",
                  "Tumebacillus", "Gracilibacillus", "Oceanobacillus",
                  "Polygonibacillus", "Salirhabdus", "Alkalihalobacillus",
                  "Halomonas"]


def main():
    print("Loading inputs ...", flush=True)
    p = pd.read_csv(CACHE / "relic_priors" /
                     "relic_score_with_mag_prior.tsv", sep="\t")
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    alive_mag = set(p.loc[p["relic_score_with_mag"] <= 0.3, "asv_id"])
    ft_a = ft.loc[ft.index.isin(alive_mag)]

    # Per-genus aggregation
    ft2 = ft_a.copy()
    ft2.index = ft2.index.rename("asv_id")
    m = ft2.reset_index().merge(tax[["asv_id", "phylum", "genus"]],
                                  on="asv_id", how="left")
    m = m.dropna(subset=["genus"])
    m = m[~m["genus"].astype(str).isin(["NA", ""])]
    sample_cols = [c for c in m.columns if c not in
                       ("asv_id", "phylum", "genus")]
    gen = m.groupby("genus")[sample_cols].sum()
    rel = relabund(gen)

    # Per-sample strategy abundances
    A_present = [g for g in STRATEGY_A if g in rel.index]
    B_present = [g for g in STRATEGY_B if g in rel.index]
    print(f"  Strategy A members present: {len(A_present)}/{len(STRATEGY_A)}",
          flush=True)
    print(f"  Strategy B members present: {len(B_present)}/{len(STRATEGY_B)}",
          flush=True)

    sA = rel.loc[A_present].sum(axis=0)
    sB = rel.loc[B_present].sum(axis=0)

    df = pd.DataFrame({"sample": sA.index,
                          "strategy_A": sA.values,
                          "strategy_B": sB.values}).merge(smeta, on="sample")
    df["log2_A_over_B"] = np.log2((df["strategy_A"] + 1e-6) /
                                       (df["strategy_B"] + 1e-6))
    df["dominant"] = np.where(df["strategy_A"] > df["strategy_B"], "A", "B")
    df.to_csv(OUT / "per_sample_strategy.tsv", sep="\t", index=False)

    # ----- 1. Per-(site, compartment, trip) summary -----
    print("\n[1] Per (site, comp, trip) strategy summary ...", flush=True)
    gby = (df.groupby(["site", "compartment", "trip"])
           .agg(median_A=("strategy_A", "median"),
                median_B=("strategy_B", "median"),
                median_log2_AoB=("log2_A_over_B", "median"),
                n=("sample", "count"))
           .reset_index())
    gby.to_csv(OUT / "per_site_comp_trip_strategy.tsv", sep="\t", index=False)

    # ----- 2. Per-(site, comp) temporal trajectory + stability -----
    print("\n[2] Stability vs switching sites ...", flush=True)
    sc = (gby.groupby(["site", "compartment"])
           .agg(n_trips=("trip", "nunique"),
                trips_dominant_A=("median_log2_AoB",
                                       lambda x: int((x > 0).sum())),
                trips_dominant_B=("median_log2_AoB",
                                       lambda x: int((x <= 0).sum())),
                AoB_min=("median_log2_AoB", "min"),
                AoB_max=("median_log2_AoB", "max"),
                AoB_range=("median_log2_AoB",
                                lambda x: float(x.max() - x.min())),
                AoB_std=("median_log2_AoB", "std"))
           .reset_index())
    sc["always_A"] = (sc["trips_dominant_A"] == sc["n_trips"]).astype(int)
    sc["always_B"] = (sc["trips_dominant_B"] == sc["n_trips"]).astype(int)
    sc["switching"] = ((sc["trips_dominant_A"] > 0) &
                          (sc["trips_dominant_B"] > 0)).astype(int)
    sc.to_csv(OUT / "per_site_comp_stability.tsv", sep="\t", index=False)

    print(f"  Total (site, comp) cells with >=2 trips: "
          f"{(sc['n_trips']>=2).sum()}", flush=True)
    print(f"  Always-A cells:    {sc['always_A'].sum()}", flush=True)
    print(f"  Always-B cells:    {sc['always_B'].sum()}", flush=True)
    print(f"  Switching cells:   {sc['switching'].sum()}", flush=True)
    print(f"  Median |AoB range| in switching cells: "
          f"{sc.loc[sc['switching']==1, 'AoB_range'].median():.2f}",
          flush=True)

    # By compartment
    print("\n  Stability by compartment:", flush=True)
    print(sc.groupby("compartment")[["always_A", "always_B",
                                          "switching"]].sum().to_string())

    # ----- 3. Per-trip global summary -----
    print("\n[3] Per-trip strategy means ...", flush=True)
    tt = (df.groupby(["trip", "compartment"])
           .agg(median_A=("strategy_A", "median"),
                median_B=("strategy_B", "median"),
                mean_log2_AoB=("log2_A_over_B", "mean"),
                pct_dominant_A=("dominant",
                                     lambda x: float((x == "A").mean() * 100)),
                n=("sample", "count"))
           .reset_index())
    print(tt.round(3).to_string(index=False))
    tt.to_csv(OUT / "per_trip_strategy.tsv", sep="\t", index=False)

    # ----- 4. Precipitation correlation -----
    print("\n[4] Precipitation correlation ...", flush=True)
    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather = weather.rename(columns={"PRECTOTCORR": "precip_mm"})

    geo_frames = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t")
            g["trip"] = trip
            geo_frames.append(g)
    geo = pd.concat(geo_frames, ignore_index=True)
    geo = geo.rename(columns={"Site": "site"})
    geo["site"] = pd.to_numeric(geo["site"], errors="coerce")
    geo = geo.dropna(subset=["site"])
    geo["site"] = geo["site"].astype(int)
    geo["sample_date"] = pd.to_datetime(geo["CenterDate"], errors="coerce")

    df = df.merge(geo[["site", "trip", "sample_date"]],
                    on=["site", "trip"], how="left")

    # Compute precip in different windows before sample_date
    rec = []
    for _, r in df.iterrows():
        if pd.isna(r["sample_date"]): continue
        ww = weather[(weather["site"] == r["site"])]
        for win in (7, 30, 90, 180, 365):
            cutoff = r["sample_date"] - pd.Timedelta(days=win)
            sub = ww[(ww["date"] <= r["sample_date"]) &
                       (ww["date"] >= cutoff)]
            rec.append({"sample": r["sample"],
                          "window": f"d{win}",
                          "precip_mm": float(sub["precip_mm"].sum())})
    pp = pd.DataFrame(rec)
    pp_pivot = pp.pivot_table(index="sample", columns="window",
                                  values="precip_mm").reset_index()
    df_w = df.merge(pp_pivot, on="sample", how="left")

    print("  Strategy A vs precip (Spearman):")
    for win in ("d7", "d30", "d90", "d180", "d365"):
        if win not in df_w.columns: continue
        sub = df_w[["strategy_A", win]].dropna()
        if len(sub) < 30: continue
        rA, _ = spearmanr(sub["strategy_A"], sub[win])
        sub2 = df_w[["strategy_B", win]].dropna()
        rB, _ = spearmanr(sub2["strategy_B"], sub2[win])
        sub3 = df_w[["log2_A_over_B", win]].dropna()
        rAoB, _ = spearmanr(sub3["log2_A_over_B"], sub3[win])
        print(f"    {win:<6}  A: rho={rA:+.3f}  B: rho={rB:+.3f}  "
              f"log2(A/B): rho={rAoB:+.3f}", flush=True)

    df_w.to_csv(OUT / "per_sample_strategy_with_precip.tsv",
                sep="\t", index=False)

    # ----- 5. Site classification map -----
    print("\n[5] Mean strategy per site (across all trips, all compartments) "
          "...", flush=True)
    per_site = (df.groupby("site")
                 .agg(median_A=("strategy_A", "median"),
                      median_B=("strategy_B", "median"),
                      median_log2_AoB=("log2_A_over_B", "median"),
                      n_samples=("sample", "count"))
                 .reset_index())
    per_site["dominant_strategy"] = np.where(
        per_site["median_log2_AoB"] > 0, "A_DOM", "B_HALO")
    per_site["strength"] = per_site["median_log2_AoB"].abs()
    per_site = per_site.sort_values("median_log2_AoB", ascending=False)
    per_site.to_csv(OUT / "per_site_strategy.tsv", sep="\t", index=False)
    print(f"  Sites dominantly A: {(per_site['dominant_strategy']=='A_DOM').sum()}")
    print(f"  Sites dominantly B: {(per_site['dominant_strategy']=='B_HALO').sum()}")
    print("  Top 10 most A-dominant sites:")
    print(per_site.head(10)[["site", "median_A", "median_B",
                                  "median_log2_AoB"]].round(3)
          .to_string(index=False))
    print("  Top 10 most B-dominant sites:")
    print(per_site.tail(10)[["site", "median_A", "median_B",
                                  "median_log2_AoB"]].round(3)
          .to_string(index=False))

    # ----- 6. Switching events -----
    print("\n[6] Switching events (sites where AoB changes sign across trips) "
          "...", flush=True)
    switching = sc[sc["switching"] == 1]
    print(f"  N switching (site, comp) cells: {len(switching)}", flush=True)
    if len(switching):
        print("  Top 10 by AoB range:")
        print(switching.sort_values("AoB_range", ascending=False)
              .head(10)[["site", "compartment", "n_trips",
                            "trips_dominant_A", "trips_dominant_B",
                            "AoB_min", "AoB_max", "AoB_range"]]
              .round(2).to_string(index=False))


if __name__ == "__main__":
    main()
