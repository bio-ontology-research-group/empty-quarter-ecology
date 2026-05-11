#!/usr/bin/env python3
"""Climate projection v2 — use TRANSITION-LEVEL longitudinal evidence
instead of cross-sectional model.

The cross-sectional model (v1) confounded geographic features with climate.
Here we use the 453 trip-to-trip transitions as the training set:
  - Outcome: did the cell switch dominance (or by how much did log2(A/B) shift)?
  - Predictors: delta_precip_d7, delta_precip_d365, delta_temp
  - Application: project under CMIP6 deltas (uniform across sites)

Two formulations:
  (i)  Probabilistic A->B switch model trained on transitions
  (ii) ΔP(B) per cell trained on per-cell precip changes

Outputs in cache/two_strategy_projection_v2/.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "two_strategy_projection_v2"
OUT.mkdir(parents=True, exist_ok=True)

# CMIP6 deltas (T in °C, precip in %)
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
    tr = pd.read_csv(CACHE / "transition_asymmetry" / "all_transitions.tsv",
                      sep="\t")
    print(f"Transitions: {len(tr)}", flush=True)
    print(f"  Distribution: {tr['transition'].value_counts().to_dict()}",
          flush=True)

    # ============================================================
    # (i) Probabilistic A->B switch model
    # ============================================================
    print("\n=== (i) Logit: P(A->B | delta_climate) on A-start transitions ===")
    A_start = tr[tr["from_dom"] == "A"].copy()
    A_start["AtoB"] = (A_start["transition"] == "A->B").astype(int)
    print(f"  N A-start transitions: {len(A_start)}", flush=True)
    print(f"  Of which A->B: {A_start['AtoB'].sum()} ({A_start['AtoB'].mean()*100:.1f}%)",
          flush=True)

    feat_cols_AtoB = ["delta_d7", "delta_d30", "delta_d90", "delta_d180",
                         "delta_d365"]
    X = A_start[feat_cols_AtoB].fillna(0).values
    y = A_start["AtoB"].values
    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    lr_AtoB = LogisticRegression(C=1.0, max_iter=2000,
                                       class_weight="balanced",
                                       random_state=42).fit(Xs, y)
    print(f"  Coefficients (z-scaled, positive = pushes toward B):")
    for c, w in zip(feat_cols_AtoB, lr_AtoB.coef_[0]):
        print(f"    {c:<15}: {w:+.3f}", flush=True)

    print("\n=== (i') Logit: P(B->A | delta_climate) on B-start transitions ===")
    B_start = tr[tr["from_dom"] == "B"].copy()
    B_start["BtoA"] = (B_start["transition"] == "B->A").astype(int)
    print(f"  N B-start transitions: {len(B_start)}", flush=True)
    print(f"  Of which B->A: {B_start['BtoA'].sum()} ({B_start['BtoA'].mean()*100:.1f}%)",
          flush=True)
    Xb = B_start[feat_cols_AtoB].fillna(0).values
    yb = B_start["BtoA"].values
    scaler_b = StandardScaler().fit(Xb)
    Xbs = scaler_b.transform(Xb)
    lr_BtoA = LogisticRegression(C=1.0, max_iter=2000,
                                       class_weight="balanced",
                                       random_state=42).fit(Xbs, yb)
    print(f"  Coefficients (z-scaled, positive = pushes toward A):")
    for c, w in zip(feat_cols_AtoB, lr_BtoA.coef_[0]):
        print(f"    {c:<15}: {w:+.3f}", flush=True)

    # ============================================================
    # (ii) Linear: delta_log2_AoB ~ delta_climate
    # ============================================================
    print("\n=== (ii) Linear: delta_log2(A/B) ~ delta_climate ===")
    Xall = tr[feat_cols_AtoB].fillna(0).values
    yall = tr["delta_log2_AoB"].fillna(0).values
    sc = StandardScaler().fit(Xall)
    XalI = sc.transform(Xall)
    lin = LinearRegression().fit(XalI, yall)
    print(f"  Coefficients (z-scaled, positive = increases A/B ratio):")
    for c, w in zip(feat_cols_AtoB, lin.coef_):
        print(f"    {c:<15}: {w:+.3f}", flush=True)
    print(f"  R^2 (in-sample): {lin.score(XalI, yall):.3f}", flush=True)

    # ============================================================
    # CMIP6 application
    # ============================================================
    print("\n=== CMIP6 projections: change in P(switch) per A-start cell ===")
    # Baseline: per-cell observed delta_d365 distribution
    # Project: delta_d365_new = delta_d365_observed + scenario_delta_P
    # We model the perturbation: dP(switch) = dP given dD365 = delta_P_pct * mean(d365_baseline)

    # Mean baseline d365 across observed cells
    mean_baseline_d365 = tr["to_d365"].dropna().mean()
    print(f"  Mean baseline d365 precip: {mean_baseline_d365:.1f} mm", flush=True)
    print(f"  Mean d7: "
          f"{tr['to_d7'].dropna().mean():.2f} mm", flush=True)

    rec = []
    for scen, (dT, dP_pct) in SCENARIOS.items():
        # Convert relative dP_pct to absolute dD365 (mm) per cell
        # dD365 = baseline_d365 * dP_pct / 100
        # We add this Δ to every cell's observed delta_d365 trajectory
        dD7 = tr["to_d7"].dropna().mean() * dP_pct / 100  # mm change
        dD365 = mean_baseline_d365 * dP_pct / 100

        # P(A->B) under projected delta climate: shift the delta vector by dD7, dD365
        X_shifted = A_start[feat_cols_AtoB].fillna(0).copy()
        X_shifted["delta_d7"] = X_shifted["delta_d7"] + dD7
        X_shifted["delta_d365"] = X_shifted["delta_d365"] + dD365
        X_shifted_arr = scaler.transform(X_shifted.values)
        p_AtoB_proj = lr_AtoB.predict_proba(X_shifted_arr)[:, 1]
        # Baseline P(A->B) without shift
        p_AtoB_base = lr_AtoB.predict_proba(scaler.transform(
            A_start[feat_cols_AtoB].fillna(0).values))[:, 1]

        # P(B->A) under projected
        Xb_shifted = B_start[feat_cols_AtoB].fillna(0).copy()
        Xb_shifted["delta_d7"] = Xb_shifted["delta_d7"] + dD7
        Xb_shifted["delta_d365"] = Xb_shifted["delta_d365"] + dD365
        Xb_shifted_arr = scaler_b.transform(Xb_shifted.values)
        p_BtoA_proj = lr_BtoA.predict_proba(Xb_shifted_arr)[:, 1]
        p_BtoA_base = lr_BtoA.predict_proba(scaler_b.transform(
            B_start[feat_cols_AtoB].fillna(0).values))[:, 1]

        # Net effect on A-start: increase in A->B rate
        rec.append({"scenario": scen, "delta_T": dT, "delta_P_pct": dP_pct,
                      "delta_d7_mm": dD7, "delta_d365_mm": dD365,
                      "P_AtoB_baseline_mean": float(p_AtoB_base.mean()),
                      "P_AtoB_projected_mean": float(p_AtoB_proj.mean()),
                      "delta_P_AtoB": float(p_AtoB_proj.mean() -
                                                p_AtoB_base.mean()),
                      "P_BtoA_baseline_mean": float(p_BtoA_base.mean()),
                      "P_BtoA_projected_mean": float(p_BtoA_proj.mean()),
                      "delta_P_BtoA": float(p_BtoA_proj.mean() -
                                                p_BtoA_base.mean())})
    summary = pd.DataFrame(rec)
    summary.to_csv(OUT / "scenario_summary_v2.tsv", sep="\t", index=False)
    print(summary.round(4).to_string(index=False))

    # ============================================================
    # Equilibrium fraction-B from the 2-state Markov chain
    # ============================================================
    print("\n=== Equilibrium B-fraction (2-state Markov from baseline rates)===")
    # Equilibrium: pi_B = P(A->B) / (P(A->B) + P(B->A))
    for _, row in summary.iterrows():
        p_AtoB = row["P_AtoB_projected_mean"]
        p_BtoA = row["P_BtoA_projected_mean"]
        pi_B = p_AtoB / (p_AtoB + p_BtoA) if p_AtoB + p_BtoA > 0 else np.nan
        print(f"  {row['scenario']:<18}  P(A->B)={p_AtoB:.3f}  "
              f"P(B->A)={p_BtoA:.3f}  ->  pi_B = "
              f"{pi_B:.3f}", flush=True)

    # Per-site projection using the baseline cell-level switching rates
    print("\n=== Per-site projected B-fraction shift ===")
    # For each cell, model uses the cell's own historical delta vector
    # We can compute per-site equilibrium pi_B at baseline vs projected
    # Group A-start by site, compute mean P(A->B) per site at baseline + each scenario
    rec_per_site = []
    for scen, (dT, dP_pct) in SCENARIOS.items():
        dD7 = tr["to_d7"].dropna().mean() * dP_pct / 100
        dD365 = mean_baseline_d365 * dP_pct / 100
        X_shifted = A_start[feat_cols_AtoB].fillna(0).copy()
        X_shifted["delta_d7"] = X_shifted["delta_d7"] + dD7
        X_shifted["delta_d365"] = X_shifted["delta_d365"] + dD365
        A_start_proj = A_start.copy()
        A_start_proj["p_AtoB"] = lr_AtoB.predict_proba(
            scaler.transform(X_shifted.values))[:, 1]
        per_site = (A_start_proj.groupby(["site", "compartment"])
                     .agg(mean_P_AtoB=("p_AtoB", "mean"),
                          n=("transition", "count"))
                     .reset_index())
        per_site["scenario"] = scen
        rec_per_site.append(per_site)
    per_site_all = pd.concat(rec_per_site, ignore_index=True)
    per_site_all.to_csv(OUT / "per_site_AtoB_risk.tsv", sep="\t", index=False)

    # Sites with highest projected risk under SSP3-7.0_2100
    high_risk = (per_site_all[per_site_all["scenario"] == "SSP3-7.0_2100"]
                  .sort_values("mean_P_AtoB", ascending=False)
                  .head(15))
    print(f"\n  Top 15 (site, comp) cells by projected P(A->B) under SSP3-7.0_2100:")
    print(high_risk.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
