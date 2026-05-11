#!/usr/bin/env python3
"""Climate projection: P(B-dominant) under CMIP6 SSP scenarios.

Steps:
  1. Build logit model: dominant ~ climate + chemistry, fit on per-cell
     observations
  2. Apply CMIP6 SSP1-2.6 / SSP2-4.5 / SSP3-7.0 deltas at 2050 + 2100
  3. Per-site projected P(B), aggregated and per-site flip estimates
  4. Identify which sites flip from A to B under each scenario

Outputs in cache/two_strategy_projection/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "two_strategy_projection"
OUT.mkdir(parents=True, exist_ok=True)

# CMIP6 deltas extracted from cache/cmip6_interventions.tsv
SCENARIOS = {
    "Historical": (0.0, 0.0),
    "SSP1-2.6_2050": (1.2, -2.0),
    "SSP1-2.6_2100": (1.6, -3.0),
    "SSP2-4.5_2050": (1.7, -3.0),
    "SSP2-4.5_2100": (2.7, -5.0),
    "SSP3-7.0_2050": (2.5, -5.0),
    "SSP3-7.0_2100": (4.0, -10.0),
}


def main():
    print("Loading inputs ...", flush=True)
    pst = pd.read_csv(CACHE / "two_strategy_temporal" /
                       "per_sample_strategy_with_precip.tsv", sep="\t")
    pst["dominant"] = np.where(pst["log2_A_over_B"] > 0, "A", "B")

    geo_t1 = pd.read_csv(DATA / "geodata" / "trip1_geodata.tsv", sep="\t")
    geo_t1 = geo_t1.rename(columns={"Site": "site"})
    geo_t1["site"] = pd.to_numeric(geo_t1["site"], errors="coerce")
    geo_t1 = geo_t1.dropna(subset=["site"])
    geo_t1["site"] = geo_t1["site"].astype(int)

    xrf = pd.read_csv(DATA / "geochemistry" /
                       "xrf_lab_table_all_trips.tsv", sep="\t")
    xrf["compartment"] = xrf["compartment"].str.lower()
    xrf["site"] = pd.to_numeric(xrf["site"], errors="coerce")
    xrf = xrf.dropna(subset=["site"])
    xrf["site"] = xrf["site"].astype(int)
    # Per-(site, comp) median XRF (sabkha proxies)
    sabkha_cols = ["SO3", "Na", "Cl", "Ca", "Mg"]
    sabkha_cols = [c for c in sabkha_cols if c in xrf.columns]
    xrf_p = (xrf[sabkha_cols + ["site", "compartment"]]
              .groupby(["site", "compartment"]).median().reset_index())

    # Build the modeling table: per (sample) with climate + chemistry
    df = (pst.merge(geo_t1[["site", "Latitude", "Longitude",
                                  "AnnualMeanTemp", "AnnualTotalPrecip"]],
                       on="site", how="left")
            .merge(xrf_p, on=["site", "compartment"], how="left"))
    df["sabkha_score"] = (df[sabkha_cols].fillna(0).sum(axis=1)
                              if sabkha_cols else 0)
    df = df.dropna(subset=["AnnualMeanTemp", "AnnualTotalPrecip",
                              "sabkha_score"])
    print(f"  Modeling samples: {len(df)}", flush=True)
    print(f"  A: {(df['dominant']=='A').sum()}  "
          f"B: {(df['dominant']=='B').sum()}", flush=True)

    # ----- Build logit model -----
    feat_cols = ["AnnualMeanTemp", "AnnualTotalPrecip", "sabkha_score",
                   "Latitude", "Longitude"]
    X = df[feat_cols].values
    y = (df["dominant"] == "B").astype(int).values
    print(f"\n  P(B) baseline: {y.mean():.3f}", flush=True)

    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    lr = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                              random_state=42)
    p_cv = cross_val_predict(lr, Xs, y, cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(y, p_cv)
    print(f"  P(B|climate+chemistry) GLM CV AUC: {auc:.3f}", flush=True)
    lr.fit(Xs, y)
    print(f"  Coefficients (z-scaled):")
    for col, c in zip(feat_cols, lr.coef_[0]):
        print(f"    {col:<22}: {c:+.3f}", flush=True)

    # ----- Apply CMIP6 deltas -----
    print("\n=== CMIP6 projections ===", flush=True)
    rec = []
    rec_per_site = []
    for scen, (dT, dP_pct) in SCENARIOS.items():
        df_proj = df.copy()
        df_proj["AnnualMeanTemp"] = df_proj["AnnualMeanTemp"] + dT
        df_proj["AnnualTotalPrecip"] = (
            df_proj["AnnualTotalPrecip"] * (1 + dP_pct / 100))
        Xp = scaler.transform(df_proj[feat_cols].values)
        p_proj = lr.predict_proba(Xp)[:, 1]
        df_proj["p_B"] = p_proj
        df_proj["scenario"] = scen

        # Per-sample projected dominant
        df_proj["proj_dominant"] = np.where(p_proj > 0.5, "B", "A")
        # Compare to current dominant
        flip_AtoB = ((df["dominant"] == "A") & (df_proj["proj_dominant"] == "B")).sum()
        flip_BtoA = ((df["dominant"] == "B") & (df_proj["proj_dominant"] == "A")).sum()
        n_total = len(df)
        n_proj_B = (p_proj > 0.5).sum()
        n_curr_B = (y == 1).sum()
        rec.append({"scenario": scen, "delta_T": dT, "delta_P_pct": dP_pct,
                      "n_samples": n_total,
                      "n_currently_B": int(n_curr_B),
                      "n_projected_B": int(n_proj_B),
                      "n_flip_AtoB": int(flip_AtoB),
                      "n_flip_BtoA": int(flip_BtoA),
                      "frac_projected_B": float(n_proj_B / n_total),
                      "median_p_B": float(np.median(p_proj))})

        # Per-site projection
        per_site = (df_proj.groupby("site")
                      .agg(p_B_median=("p_B", "median"),
                            curr_dominant_mode=("dominant",
                                                  lambda x: x.mode().iloc[0]),
                            proj_dominant_mode=("proj_dominant",
                                                  lambda x: x.mode().iloc[0]))
                      .reset_index())
        per_site["scenario"] = scen
        per_site["flip"] = ((per_site["curr_dominant_mode"] == "A") &
                              (per_site["proj_dominant_mode"] == "B")).astype(int)
        rec_per_site.append(per_site)

    summary = pd.DataFrame(rec)
    summary.to_csv(OUT / "scenario_summary.tsv", sep="\t", index=False)
    print(summary.round(3).to_string(index=False))

    per_site_all = pd.concat(rec_per_site, ignore_index=True)
    per_site_all.to_csv(OUT / "per_site_projection.tsv", sep="\t", index=False)

    # Per-site: which sites flip across scenarios?
    print("\n=== Per-site flip risk ===")
    flip_pivot = per_site_all.pivot_table(index="site", columns="scenario",
                                                values="proj_dominant_mode",
                                                aggfunc="first")
    flip_pivot["scenarios_flipped"] = flip_pivot.apply(
        lambda row: (row[1:] == "B").sum() if row.iloc[0] == "A" else 0,
        axis=1)
    print(f"\n  Sites currently A that flip in >=1 scenario:")
    high_risk = flip_pivot[flip_pivot["scenarios_flipped"] >= 1]
    print(high_risk.to_string())
    print(f"\n  N sites currently A: "
          f"{(flip_pivot['Historical']=='A').sum()}", flush=True)
    print(f"  N sites flipping under SSP1-2.6_2100: "
          f"{((flip_pivot['Historical']=='A') & (flip_pivot['SSP1-2.6_2100']=='B')).sum()}",
          flush=True)
    print(f"  N sites flipping under SSP2-4.5_2100: "
          f"{((flip_pivot['Historical']=='A') & (flip_pivot['SSP2-4.5_2100']=='B')).sum()}",
          flush=True)
    print(f"  N sites flipping under SSP3-7.0_2100: "
          f"{((flip_pivot['Historical']=='A') & (flip_pivot['SSP3-7.0_2100']=='B')).sum()}",
          flush=True)
    flip_pivot.to_csv(OUT / "per_site_scenarios.tsv", sep="\t")


if __name__ == "__main__":
    main()
