#!/usr/bin/env python3
"""Tier-1 #2: Time-for-space substitution test (trip-resolved).

Per-ASV environment-response model fit on T1-T4 samples (4 seasonal
"snapshots"), then evaluated on T5 to assess transferability. Per-genus
predicted shifts under CMIP6 SSP1-2.6 vs SSP5-8.5 deltas (already in
cache/cmip6_interventions.tsv) compared to observed seasonal range.

Inputs:
  cache/feature_table.parquet
  cache/taxonomy.parquet
  data/climate/daily_weather_full.csv     full NASA POWER w/ wind, RH
  data/geodata/trip{1..5}_geodata.tsv     CenterDate per (site, trip)
  cache/cmip6_interventions.tsv           projected delta_T_C, delta_P_pct

Outputs:
  cache/tfs/asv_env_coefs.tsv
  cache/tfs/cv_t14_to_t5_per_asv.tsv
  cache/tfs/cmip6_predicted_shifts_per_genus.tsv
  cache/tfs/seasonal_range_observed_vs_cmip6.tsv
  cache/tfs/summary.txt
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
OUT = CACHE / "tfs"
OUT.mkdir(parents=True, exist_ok=True)

PREV_MIN = 0.05


def load_geo() -> pd.DataFrame:
    geo = []
    for trip in range(1, 6):
        gp = DATA / "geodata" / f"trip{trip}_geodata.tsv"
        if gp.exists():
            g = pd.read_csv(gp, sep="\t")
            g["trip"] = trip
            g["CenterDate"] = pd.to_datetime(g["CenterDate"])
            g["Site_int"] = pd.to_numeric(g["Site"], errors="coerce")
            g = g.dropna(subset=["Site_int", "CenterDate"])
            g["site"] = g["Site_int"].astype(int)
            geo.append(g[["site", "trip", "CenterDate"]])
    return pd.concat(geo).drop_duplicates(["site", "trip"])


def attach_env(smeta: pd.DataFrame, weather: pd.DataFrame) -> pd.DataFrame:
    """Attach per-sample 30-day mean T2M, RH2M, sum PRECTOTCORR before sample_date."""
    w = weather.set_index(["site", "date"]).sort_index()
    rec = []
    for _, r in smeta.iterrows():
        site = int(r["site"])
        d = pd.Timestamp(r["sample_date"])
        try:
            site_w = w.loc[site]
        except KeyError:
            rec.append({"sample": r["sample"], "T_30d": np.nan, "RH_30d": np.nan,
                        "P_30d": np.nan})
            continue
        mask = (site_w.index >= d - pd.Timedelta(days=29)) & (site_w.index <= d)
        win = site_w.loc[mask]
        rec.append({"sample": r["sample"],
                    "T_30d":  float(win["T2M"].mean()) if len(win) else np.nan,
                    "RH_30d": float(win["RH2M"].mean()) if len(win) else np.nan,
                    "P_30d":  float(win["PRECTOTCORR"].sum()) if len(win) else np.nan})
    return smeta.merge(pd.DataFrame(rec), on="sample", how="left")


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet")
    print(f"feature_table: {ft.shape}", flush=True)

    smeta = parse_samples_to_df(ft.columns)
    geo = load_geo()
    smeta = smeta.merge(geo, on=["site", "trip"], how="left")
    smeta = smeta.dropna(subset=["CenterDate"])
    smeta["sample_date"] = smeta["CenterDate"]
    print(f"samples after date join: {len(smeta)}", flush=True)

    weather = pd.read_csv(DATA / "climate" / "daily_weather_full.csv")
    weather["date"] = pd.to_datetime(weather["date"], format="%Y%m%d")
    weather["site"] = weather["site"].astype(int)

    smeta = attach_env(smeta, weather[["site", "date", "T2M", "RH2M", "PRECTOTCORR"]])
    smeta = smeta.dropna(subset=["T_30d", "RH_30d", "P_30d"])
    print(f"samples after env attach: {len(smeta)}", flush=True)
    print(smeta[["T_30d", "RH_30d", "P_30d"]].describe().to_string())

    rel = ft.div(ft.sum(axis=0).replace(0, 1), axis=1)
    prev = (ft > 0).sum(axis=1) / ft.shape[1]
    keep = prev[prev >= PREV_MIN].index
    print(f"ASVs prev>={PREV_MIN}: {len(keep)}", flush=True)
    rel_k = rel.loc[keep]

    # ---- Per-ASV linear model trained on ALL samples (with compartment FE)
    eps = 1e-6
    coef_rows = []
    env_idx = smeta.set_index("sample")
    for asv in keep:
        x = rel_k.loc[asv].reindex(env_idx.index)
        df = pd.DataFrame({"y": np.log10(x.values + eps),
                           "T": env_idx["T_30d"].values,
                           "RH": env_idx["RH_30d"].values,
                           "P": env_idx["P_30d"].values,
                           "comp": env_idx["compartment"].values}).dropna()
        if len(df) < 30:
            continue
        for comp in ["rhizosphere", "surface", "deep"]:
            sub = df.loc[df["comp"] == comp]
            if len(sub) < 20:
                continue
            X = np.column_stack([np.ones(len(sub)),
                                  sub["T"].values, sub["RH"].values, sub["P"].values])
            try:
                beta, *_ = np.linalg.lstsq(X, sub["y"].values, rcond=None)
                pred = X @ beta
                ss_res = np.sum((sub["y"].values - pred) ** 2)
                ss_tot = np.sum((sub["y"].values - sub["y"].mean()) ** 2)
                r2 = 1 - ss_res / max(ss_tot, 1e-12)
                coef_rows.append({"ASV": asv, "compartment": comp, "n_obs": len(sub),
                                  "beta_T": float(beta[1]),
                                  "beta_RH": float(beta[2]),
                                  "beta_P": float(beta[3]),
                                  "intercept": float(beta[0]),
                                  "R2": float(r2)})
            except Exception:
                continue
    coefs = pd.DataFrame(coef_rows)
    coefs.to_csv(OUT / "asv_env_coefs.tsv", sep="\t", index=False)
    print(f"per-ASV coefs: {len(coefs)} rows", flush=True)

    # ---- Cross-validation: train on T1-T4, predict T5
    train_samples = smeta.loc[smeta["trip"] != 5, "sample"].tolist()
    test_samples  = smeta.loc[smeta["trip"] == 5, "sample"].tolist()
    print(f"\nCV split: train T1-4 n={len(train_samples)}, test T5 n={len(test_samples)}",
          flush=True)
    env_train = smeta.set_index("sample").reindex(train_samples)
    env_test  = smeta.set_index("sample").reindex(test_samples)

    cv_rows = []
    for asv in keep:
        y_tr = np.log10(rel_k.loc[asv].reindex(train_samples).values + eps)
        y_te = np.log10(rel_k.loc[asv].reindex(test_samples).values + eps)
        X_tr = np.column_stack([np.ones(len(env_train)),
                                 env_train["T_30d"].values,
                                 env_train["RH_30d"].values,
                                 env_train["P_30d"].values])
        X_te = np.column_stack([np.ones(len(env_test)),
                                 env_test["T_30d"].values,
                                 env_test["RH_30d"].values,
                                 env_test["P_30d"].values])
        try:
            beta, *_ = np.linalg.lstsq(X_tr, y_tr, rcond=None)
            pred_te = X_te @ beta
            r, p = spearmanr(pred_te, y_te)
            cv_rows.append({"ASV": asv, "rho_pred_obs": float(r) if r == r else np.nan,
                            "p_value": float(p) if p == p else np.nan,
                            "n_train": len(train_samples), "n_test": len(test_samples)})
        except Exception:
            continue
    cv_df = pd.DataFrame(cv_rows)
    cv_df.to_csv(OUT / "cv_t14_to_t5_per_asv.tsv", sep="\t", index=False)
    if len(cv_df):
        med_rho = cv_df["rho_pred_obs"].median()
        frac_p = (cv_df["p_value"] < 0.05).mean()
        print(f"CV median rho(pred,obs)={med_rho:+.3f}; "
              f"frac p<0.05 = {frac_p:.1%}", flush=True)

    # ---- CMIP6 predicted shifts per genus
    cmip = pd.read_csv(CACHE / "cmip6_interventions.tsv", sep="\t")
    coefs_g = coefs.copy()
    coefs_g["genus"] = coefs_g["ASV"].map(tax["genus"].fillna("Unclassified"))
    g_summary = (coefs_g.groupby(["compartment", "genus"])
                 .agg(n_asv=("ASV", "count"),
                      median_betaT=("beta_T", "median"),
                      median_betaRH=("beta_RH", "median"),
                      median_betaP=("beta_P", "median"),
                      median_R2=("R2", "median"))
                 .reset_index())

    pred = []
    for _, c in cmip.iterrows():
        comp, scen, h = c["compartment"], c["scenario"], c["horizon"]
        dT, dPpct = float(c["delta_T_C"]), float(c["delta_P_pct"])
        for _, r in g_summary[g_summary["compartment"] == comp].iterrows():
            dlog10 = r["median_betaT"] * dT  # P channel ignored unless we scale
            pred.append({"compartment": comp, "scenario": scen, "horizon": h,
                         "genus": r["genus"], "n_asv": int(r["n_asv"]),
                         "predicted_log10_shift": float(dlog10),
                         "median_betaT": float(r["median_betaT"]),
                         "median_R2": float(r["median_R2"])})
    pred_df = pd.DataFrame(pred)
    pred_df.to_csv(OUT / "cmip6_predicted_shifts_per_genus.tsv", sep="\t", index=False)

    # ---- Observed seasonal range vs CMIP6 magnitude
    season_range = []
    smeta_t = smeta.set_index("sample")
    trip_temps = smeta.groupby("trip")["T_30d"].mean()
    obs_T_range = float(trip_temps.max() - trip_temps.min())
    for _, c in cmip.iterrows():
        season_range.append({"scenario": c["scenario"], "horizon": c["horizon"],
                             "compartment": c["compartment"],
                             "delta_T_C": float(c["delta_T_C"]),
                             "obs_seasonal_T_range_C": obs_T_range,
                             "ratio_cmip_to_seasonal": float(c["delta_T_C"]) / obs_T_range
                             if obs_T_range > 0 else np.nan})
    sr = pd.DataFrame(season_range)
    sr.to_csv(OUT / "seasonal_range_observed_vs_cmip6.tsv", sep="\t", index=False)

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Time-for-space substitution test (Tier-1 #2)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Trip prefix mapping in use: 1=(none) 2=T 3=F 4=S 5=V\n\n")
        fh.write(f"samples after date join + env attach: {len(smeta)}\n")
        fh.write(f"per-trip: {smeta['trip'].value_counts().sort_index().to_dict()}\n")
        fh.write(f"per-compartment: "
                 f"{smeta['compartment'].value_counts().to_dict()}\n\n")
        fh.write(f"observed seasonal T range (mean per trip): "
                 f"{obs_T_range:.2f} °C\n")
        if len(cv_df):
            fh.write(f"\nT1-4 -> T5 cross-validation:\n")
            fh.write(f"  ASVs tested: {len(cv_df)}\n")
            fh.write(f"  median rho(pred, obs): "
                     f"{cv_df['rho_pred_obs'].median():+.3f}\n")
            fh.write(f"  fraction p<0.05: "
                     f"{(cv_df['p_value']<0.05).mean():.1%}\n")
            fh.write(f"  fraction with rho>0: "
                     f"{(cv_df['rho_pred_obs']>0).mean():.1%}\n")
        fh.write("\nTop 10 CMIP6-projected genus winners (SSP3-7.0, "
                 "rhizosphere, 2050):\n")
        sub = pred_df[(pred_df["scenario"] == "SSP3-7.0") &
                      (pred_df["horizon"].astype(str) == "2050") &
                      (pred_df["compartment"] == "rhizosphere")] \
                .sort_values("predicted_log10_shift", ascending=False)
        fh.write(sub.head(10).to_string(index=False))
        fh.write("\n\nTop 10 losers:\n")
        fh.write(sub.tail(10).to_string(index=False))
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
