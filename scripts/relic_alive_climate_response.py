#!/usr/bin/env python3
"""Climate-trend / time-for-space / pulse-reserve genus-level reanalysis on
alive vs relic vs all subsets. Adds:

  A. Per-site per-trip alive Shannon vs aridity / temp / precip
  B. Time-for-space: T1-4 mean alive composition predicts T5? (alive-only)
  C. Pulse-reserve genus-level (top responders in alive vs all)
  D. Climate-change-projection-relevant: per-site relic_frac vs annual T

Outputs:
  cache/relic_alive_subset/climate_response_<pool>.tsv
  cache/relic_alive_subset/time_for_space_<pool>.tsv
  cache/relic_alive_subset/genus_pulse_responders_<pool>.tsv
  cache/relic_alive_subset/climate_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "relic_alive_subset"
OUT.mkdir(parents=True, exist_ok=True)


def relabund(M):
    return M.div(M.sum(axis=0).replace(0, 1), axis=1)


def shannon(arr):
    a = arr[arr > 0]
    if len(a) == 0: return 0.0
    p = a / a.sum()
    return float(-(p * np.log(p)).sum())


def main():
    print("Loading inputs ...", flush=True)
    ft_all = pd.read_parquet(CACHE / "feature_table.parquet")
    ft_alive = pd.read_parquet(CACHE / "feature_table_alive.parquet")
    ft_relic = pd.read_parquet(CACHE / "feature_table_relic.parquet")
    smeta = parse_samples_to_df(ft_all.columns)
    smeta["site"] = smeta["site"].astype(int)

    # Per-site trip-1 climate
    geo_frames = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t"); g["trip"] = trip
            geo_frames.append(g)
    geo = pd.concat(geo_frames, ignore_index=True)
    geo = geo.rename(columns={"Site": "site"})
    geo["site"] = pd.to_numeric(geo["site"], errors="coerce")
    geo = geo.dropna(subset=["site"]); geo["site"] = geo["site"].astype(int)
    geo_t1 = geo[geo["trip"] == 1].drop_duplicates("site")[
        ["site", "Latitude", "Longitude", "AnnualMeanTemp",
          "AnnualTotalPrecip"]]

    summary_rows = []

    # ------------------------------------------------------------------
    # A. Per-(sample, pool) Shannon vs climate
    # ------------------------------------------------------------------
    print("\n=== A. Sample Shannon vs climate ===", flush=True)
    for label, ft in (("all", ft_all), ("alive", ft_alive),
                          ("relic", ft_relic)):
        sh = pd.DataFrame({"sample": ft.columns,
                              "shannon": [shannon(ft[c].values)
                                            for c in ft.columns]})
        sh = sh.merge(smeta, on="sample", how="left")
        sh = sh.merge(geo_t1, on="site", how="left")
        sh.to_csv(OUT / f"climate_response_{label}.tsv", sep="\t", index=False)

        for v in ("Latitude", "Longitude", "AnnualMeanTemp",
                    "AnnualTotalPrecip"):
            ss = sh[["shannon", v]].dropna()
            if len(ss) < 30: continue
            r, p = spearmanr(ss["shannon"], ss[v])
            summary_rows.append({"analysis": "shannon_vs_climate",
                                    "pool": label, "variable": v,
                                    "rho": float(r), "p": float(p),
                                    "n": len(ss)})
            print(f"  {label:<5}  shannon ~ {v:<20}  rho={r:+.3f}  p={p:.3g}",
                  flush=True)

    # ------------------------------------------------------------------
    # B. Time-for-space: T1-4 mean predicts T5 (Bray-Curtis distance)
    # ------------------------------------------------------------------
    from scipy.spatial.distance import braycurtis
    print("\n=== B. Time-for-space (T1-4 -> T5) ===", flush=True)
    for label, ft in (("all", ft_all), ("alive", ft_alive),
                          ("relic", ft_relic)):
        rows = []
        Mr = relabund(ft)
        for site in sorted(smeta["site"].unique()):
            for comp in ("rhizosphere", "surface", "deep"):
                t14 = smeta[(smeta["site"] == site) &
                              (smeta["compartment"] == comp) &
                              (smeta["trip"].between(1, 4))]
                t5  = smeta[(smeta["site"] == site) &
                              (smeta["compartment"] == comp) &
                              (smeta["trip"] == 5)]
                t14_s = list(set(t14["sample"]) & set(Mr.columns))
                t5_s  = list(set(t5["sample"]) & set(Mr.columns))
                if len(t14_s) < 2 or len(t5_s) < 1: continue
                t14_mean = Mr[t14_s].mean(axis=1).values
                for s5 in t5_s:
                    bc = braycurtis(t14_mean, Mr[s5].values)
                    rows.append({"site": site, "compartment": comp,
                                  "t5_sample": s5, "n_t14_samples": len(t14_s),
                                  "bc_t5_vs_t14_mean": float(bc)})
        df = pd.DataFrame(rows)
        if len(df) == 0: continue
        df.to_csv(OUT / f"time_for_space_{label}.tsv", sep="\t", index=False)
        med = float(df["bc_t5_vs_t14_mean"].median())
        print(f"  {label:<5}  median BC(T5 vs mean T1-4) = {med:.3f}  "
              f"(n={len(df)})", flush=True)
        summary_rows.append({"analysis": "time_for_space", "pool": label,
                                "metric": "median_BC_t5_vs_t14",
                                "value": med, "n": len(df)})

    # ------------------------------------------------------------------
    # C. Pulse-reserve genus-level (alive vs all)
    # ------------------------------------------------------------------
    print("\n=== C. Pulse-reserve genus-level ===", flush=True)
    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather = weather.rename(columns={"PRECTOTCORR": "precip_mm"})

    sm = smeta.merge(geo[["site", "trip", "CenterDate"]],
                       on=["site", "trip"], how="left")
    sm["sample_date"] = pd.to_datetime(sm["CenterDate"], errors="coerce")
    sm = sm.dropna(subset=["sample_date"])
    state_rows = []
    for _, r in sm.iterrows():
        ww = weather[(weather["site"] == r["site"]) &
                       (weather["date"] <= r["sample_date"]) &
                       (weather["date"] >= r["sample_date"] -
                          pd.Timedelta(days=30))]
        precip30 = float(ww["precip_mm"].sum())
        state_rows.append({"sample": r["sample"], "site": r["site"],
                            "compartment": r["compartment"],
                            "trip": r["trip"], "precip30": precip30,
                            "state": "pulse" if precip30 > 1 else "reserve"})
    state = pd.DataFrame(state_rows)

    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})
    for label, ft in (("all", ft_all), ("alive", ft_alive)):
        # Aggregate to genus relabund
        ft2 = ft.copy()
        ft2.index = ft2.index.rename("asv_id")
        m = ft2.reset_index()
        m = m.merge(tax[["asv_id", "genus"]], on="asv_id", how="left")
        m = m.dropna(subset=["genus"])
        sample_cols = [c for c in m.columns if c not in ("asv_id", "genus")]
        gen = m.groupby("genus")[sample_cols].sum()
        gen_rel = relabund(gen)

        rows = []
        for g in gen_rel.index:
            v = gen_rel.loc[g]
            pulse_v = v[state.loc[state["state"] == "pulse",
                                    "sample"]].dropna()
            reserve_v = v[state.loc[state["state"] == "reserve",
                                      "sample"]].dropna()
            if len(pulse_v) < 30 or len(reserve_v) < 30: continue
            try:
                U, p = mannwhitneyu(pulse_v, reserve_v,
                                       alternative="greater")
            except Exception:
                continue
            mp = float(pulse_v.median())
            mr = float(reserve_v.median())
            rows.append({"genus": g, "median_pulse": mp,
                          "median_reserve": mr,
                          "log2_pulse_reserve": float(np.log2(
                              (mp + 1e-9) / (mr + 1e-9))),
                          "p": float(p)})
        df = pd.DataFrame(rows).sort_values("p")
        if len(df) == 0: continue
        n = len(df); p = df["p"].values; idx = np.argsort(p)
        ranked = np.empty(n); ranked[idx] = np.arange(1, n + 1)
        df["q"] = np.minimum(p * n / ranked, 1.0)
        df.to_csv(OUT / f"genus_pulse_responders_{label}.tsv",
                   sep="\t", index=False)
        sig = df[df["q"] < 0.05]
        print(f"  {label:<5} genus tested={len(df)}  q<0.05={len(sig)}",
              flush=True)
        summary_rows.append({"analysis": "genus_pulse_responders",
                                "pool": label, "n_tested": len(df),
                                "n_sig": len(sig)})
        # Top 10
        print(f"  Top 10 pulse-responder genera ({label}):")
        print(df.head(10)[["genus", "median_pulse", "median_reserve",
                              "log2_pulse_reserve", "q"]].round(5)
              .to_string(index=False), flush=True)

    # ------------------------------------------------------------------
    # D. Per-site relic_frac vs annual climate
    # ------------------------------------------------------------------
    print("\n=== D. Per-site relic_frac vs climate (already done in "
          "geographic_correlates) ===", flush=True)

    # Save summary
    smdf = pd.DataFrame(summary_rows)
    smdf.to_csv(OUT / "climate_response_summary.tsv", sep="\t", index=False)
    print(f"\nWrote summary to {OUT}/climate_response_summary.tsv")


if __name__ == "__main__":
    main()
