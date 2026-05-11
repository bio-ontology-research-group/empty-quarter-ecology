#!/usr/bin/env python3
"""Composite relic-likelihood indicator (Path B).

Trains a model on PMA-anchored ASVs (T/UT viability as ground truth) and
extrapolates a relic-likelihood score to all ~75k ASVs.

Features (per EQ ASV):
  - persistence_max:   max n_trips_present across (site, comp) cells
  - persistence_mean:  mean n_trips_present
  - log_mean_abund:    log10 of mean per-sample relative abundance (where >0)
  - log_max_abund:     log10 of max abundance
  - n_sites_detected:  in how many of 60 sites was the ASV ever detected
  - frac_deep:         fraction of detections that occurred in deep compartment
  - frac_surface:      "" surface
  - frac_rhizo:        "" rhizosphere
  - emp_cosmo_90:      ASV has a min-25 EMP match at 90bp (0/1)
  - emp_cosmo_100:     "" 100bp (0/1)
  - emp_cosmo_150:     "" 150bp (0/1)

Ground truth: weighted_median_ratio from per_eq_asv_viability.tsv.
  Relic-positive label  if T/UT < 0.1
  Relic-negative label  if T/UT > 0.5
  ASVs with 0.1 <= T/UT <= 0.5 are "ambiguous" -- excluded from training.

Model: logistic regression (L2) + 5-fold CV.
Also fits a gradient-boosted classifier for a nonlinear comparison.

Outputs:
  cache/test6_disconfirmation/relic_features.parquet
  cache/test6_disconfirmation/relic_model_metrics.txt
  cache/test6_disconfirmation/relic_indicator_per_asv.tsv
  cache/test6_disconfirmation/relic_calibration.tsv
  cache/test6_disconfirmation/relic_feature_importance.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, average_precision_score,
                                brier_score_loss, log_loss,
                                roc_curve)
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df
CACHE = REPO / "cache"
OUT = CACHE / "test6_disconfirmation"

RELIC_T_LO = 0.1   # T/UT < 0.1  -> relic
ALIVE_T_HI = 0.5   # T/UT > 0.5  -> alive


def build_features() -> pd.DataFrame:
    """Per-ASV feature table covering all detected ASVs."""
    print("Loading feature table ...", flush=True)
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    print(f"  shape: {ft.shape}", flush=True)
    smeta = parse_samples_to_df(ft.columns)
    # Strip header cols if any
    smeta = smeta.set_index("sample")

    # Per-ASV: relative abundance per sample
    print("Normalizing to relative abundance per sample ...", flush=True)
    sample_sums = ft.sum(axis=0).replace(0, 1)
    ft_rel = ft.div(sample_sums, axis=1)

    # Mean / max relative abundance (where >0)
    print("Computing per-ASV abundance summaries ...", flush=True)
    pos_count = (ft_rel > 0).sum(axis=1)
    sum_rel = ft_rel.sum(axis=1)
    max_rel = ft_rel.max(axis=1)
    mean_when_present = sum_rel / pos_count.replace(0, 1)

    # Compartment fractions: of the times this ASV is detected, what fraction
    # of detections were in each compartment?
    print("Computing compartment fractions ...", flush=True)
    comp_map = smeta["compartment"]
    presence = (ft > 0).astype(np.int8)
    site_map = smeta["site"].astype(int)

    # Sites detected (uniqueness)
    site_per_sample = site_map.values
    # For each ASV: set of sites where presence>0
    print("  building site/comp fraction arrays ...", flush=True)
    feats = pd.DataFrame(index=ft.index)
    feats["mean_when_present"] = mean_when_present
    feats["max_abund"] = max_rel
    feats["n_detections"] = pos_count
    feats["log_mean_abund"] = np.log10(mean_when_present.replace(0, np.nan))
    feats["log_max_abund"] = np.log10(max_rel.replace(0, np.nan))

    # Fraction of detections in each compartment
    for comp in ("rhizosphere", "surface", "deep"):
        cols = comp_map[comp_map == comp].index
        cols = [c for c in cols if c in presence.columns]
        if not cols:
            feats[f"frac_{comp}"] = 0.0
            continue
        det = presence[cols].sum(axis=1)
        feats[f"frac_{comp}"] = (det / pos_count.replace(0, 1)).fillna(0)

    # n distinct sites
    print("  counting distinct sites per ASV ...", flush=True)
    site_per_col = pd.Series(site_per_sample, index=presence.columns)
    # build sparse-friendly: for each ASV, n_unique sites in the cols where >0
    # Do it column-major: aggregate site presence indicators
    # Simpler: for each site, compute presence-any across its samples
    by_site = {}
    for s, cols in site_per_col.groupby(site_per_col):
        cols_list = list(cols.index)
        by_site[s] = (presence[cols_list].sum(axis=1) > 0).astype(np.int8)
    site_present = pd.DataFrame(by_site)
    feats["n_sites_detected"] = site_present.sum(axis=1)
    print(f"  feats shape: {feats.shape}", flush=True)
    return feats


def add_persistence(feats: pd.DataFrame) -> pd.DataFrame:
    print("Adding persistence ...", flush=True)
    persist = pd.read_parquet(CACHE / "test6_persistence" /
                                "per_OTU_site_persistence.parquet")
    grp = (persist.groupby("ASV")
           .agg(persistence_max=("n_trips_present", "max"),
                persistence_mean=("n_trips_present", "mean"),
                n_site_records=("n_trips_present", "count"))
           .reset_index().rename(columns={"ASV": "asv_id"}))
    feats = feats.copy()
    feats.index = feats.index.rename("asv_id")
    feats = feats.reset_index()
    feats = feats.merge(grp, on="asv_id", how="left")
    feats[["persistence_max", "persistence_mean"]] = (
        feats[["persistence_max", "persistence_mean"]].fillna(1))
    feats["n_site_records"] = feats["n_site_records"].fillna(0)
    return feats


def add_emp(feats: pd.DataFrame) -> pd.DataFrame:
    print("Adding EMP min25 cosmopolitanism flags ...", flush=True)
    for length in (90, 100, 150):
        path = CACHE / "emp_cosmopolitanism" / f"eq_vs_emp_{length}bp.tsv"
        if not path.exists():
            feats[f"emp_cosmo_{length}"] = 0
            continue
        d = pd.read_csv(path, sep="\t", header=None,
                         usecols=[0], names=["asv_id"])
        d = d["asv_id"].astype(str).unique()
        feats[f"emp_cosmo_{length}"] = (
            feats["asv_id"].isin(d).astype(int))
        print(f"  emp_cosmo_{length}bp hits: "
              f"{feats[f'emp_cosmo_{length}'].sum()}", flush=True)
    return feats


def add_pma_viability(feats: pd.DataFrame) -> pd.DataFrame:
    print("Adding PMA viability ground truth ...", flush=True)
    viab = pd.read_csv(OUT / "per_eq_asv_viability.tsv", sep="\t")
    viab = viab.rename(columns={"eq_asv": "asv_id"})
    feats = feats.merge(viab[["asv_id", "weighted_median_ratio",
                                "median_ratio", "n_pma_proxies", "max_pid"]],
                          on="asv_id", how="left")
    return feats


def main():
    feats = build_features()
    feats = add_persistence(feats)
    feats = add_emp(feats)
    feats = add_pma_viability(feats)
    feats.to_parquet(OUT / "relic_features.parquet")
    print(f"\nFeature table written: {feats.shape}", flush=True)

    # Define training set: ASVs with PMA viability and falling outside the
    # ambiguous middle band
    train_mask = (feats["weighted_median_ratio"].notna() &
                   ((feats["weighted_median_ratio"] < RELIC_T_LO) |
                    (feats["weighted_median_ratio"] > ALIVE_T_HI)))
    df = feats[train_mask].copy()
    df["y_relic"] = (df["weighted_median_ratio"] < RELIC_T_LO).astype(int)
    print(f"\nTraining set:", flush=True)
    print(f"  total PMA-anchored ASVs:        "
          f"{feats['weighted_median_ratio'].notna().sum()}", flush=True)
    print(f"  outside ambiguous band (kept):  {len(df)}", flush=True)
    print(f"  relic class (y=1) count:        {int(df['y_relic'].sum())}",
          flush=True)
    print(f"  alive class (y=0) count:        "
          f"{int((1-df['y_relic']).sum())}", flush=True)

    # Feature columns
    fcols = [
        "persistence_max", "persistence_mean", "n_site_records",
        "log_mean_abund", "log_max_abund", "n_detections", "n_sites_detected",
        "frac_deep", "frac_surface", "frac_rhizosphere",
        "emp_cosmo_90", "emp_cosmo_100", "emp_cosmo_150",
    ]
    df = df.dropna(subset=fcols).copy()
    print(f"  after dropna on features:       {len(df)}", flush=True)

    X = df[fcols].values
    y = df["y_relic"].values

    # Sanitize any inf
    X = np.where(np.isfinite(X), X, 0.0)

    # Scaler + LR
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    print("\n=== 5-fold CV: Logistic Regression (L2) ===", flush=True)
    lr = LogisticRegression(penalty="l2", C=1.0, max_iter=2000,
                              class_weight="balanced", random_state=42)
    p_lr = cross_val_predict(lr, Xs, y, cv=cv, method="predict_proba")[:, 1]
    auc_lr = roc_auc_score(y, p_lr)
    ap_lr = average_precision_score(y, p_lr)
    brier_lr = brier_score_loss(y, p_lr)
    ll_lr = log_loss(y, np.clip(p_lr, 1e-6, 1-1e-6))
    print(f"  AUC = {auc_lr:.3f}", flush=True)
    print(f"  AP  = {ap_lr:.3f}", flush=True)
    print(f"  Brier = {brier_lr:.3f}", flush=True)
    print(f"  LogLoss = {ll_lr:.3f}", flush=True)

    print("\n=== 5-fold CV: Gradient Boosting (nonlinear) ===", flush=True)
    gb = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                       learning_rate=0.05,
                                       random_state=42)
    p_gb = cross_val_predict(gb, X, y, cv=cv, method="predict_proba")[:, 1]
    auc_gb = roc_auc_score(y, p_gb)
    ap_gb = average_precision_score(y, p_gb)
    brier_gb = brier_score_loss(y, p_gb)
    ll_gb = log_loss(y, np.clip(p_gb, 1e-6, 1-1e-6))
    print(f"  AUC = {auc_gb:.3f}", flush=True)
    print(f"  AP  = {ap_gb:.3f}", flush=True)
    print(f"  Brier = {brier_gb:.3f}", flush=True)
    print(f"  LogLoss = {ll_gb:.3f}", flush=True)

    # Refit on all training data, then score every ASV
    print("\nFitting final models on all training data ...", flush=True)
    lr_final = LogisticRegression(penalty="l2", C=1.0, max_iter=2000,
                                       class_weight="balanced",
                                       random_state=42)
    lr_final.fit(Xs, y)
    gb_final = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                              learning_rate=0.05,
                                              random_state=42)
    gb_final.fit(X, y)

    # Score the entire ASV table
    print("Scoring all ASVs ...", flush=True)
    full = feats.dropna(subset=fcols).copy()
    Xall = full[fcols].values
    Xall = np.where(np.isfinite(Xall), Xall, 0.0)
    Xall_s = sc.transform(Xall)
    full["relic_score_lr"] = lr_final.predict_proba(Xall_s)[:, 1]
    full["relic_score_gb"] = gb_final.predict_proba(Xall)[:, 1]

    out_cols = ["asv_id"] + fcols + [
        "weighted_median_ratio", "n_pma_proxies",
        "relic_score_lr", "relic_score_gb",
    ]
    full[out_cols].to_csv(OUT / "relic_indicator_per_asv.tsv",
                          sep="\t", index=False)
    print(f"  wrote {len(full)} ASV rows", flush=True)

    # Distribution of relic_score for PMA-anchored set vs all
    print("\nrelic_score_gb distribution (all ASVs):")
    for q in (10, 25, 50, 75, 90):
        print(f"  p{q}: {np.percentile(full['relic_score_gb'], q):.3f}",
              flush=True)
    print("\nrelic_score_gb distribution (PMA-anchored, training subset):")
    sub = full[full["asv_id"].isin(df["asv_id"])]
    for q in (10, 25, 50, 75, 90):
        print(f"  p{q}: {np.percentile(sub['relic_score_gb'], q):.3f}",
              flush=True)

    # Calibration: bin scores into deciles, report mean predicted vs mean
    # observed (only on training set where ground truth exists)
    df_score = df.merge(full[["asv_id", "relic_score_lr", "relic_score_gb"]],
                          on="asv_id", how="inner")
    cal = []
    for k in range(10):
        lo, hi = k / 10, (k + 1) / 10
        m = (df_score["relic_score_gb"] >= lo) & \
            (df_score["relic_score_gb"] < hi if k < 9 else
             df_score["relic_score_gb"] <= hi)
        if m.sum() == 0: continue
        cal.append({"bin_lo": lo, "bin_hi": hi, "n": int(m.sum()),
                     "mean_pred": float(df_score.loc[m, "relic_score_gb"].mean()),
                     "mean_obs":  float(df_score.loc[m, "y_relic"].mean())})
    cal_df = pd.DataFrame(cal)
    cal_df.to_csv(OUT / "relic_calibration.tsv", sep="\t", index=False)
    print("\nCalibration (gradient boost):")
    print(cal_df.to_string(index=False), flush=True)

    # Feature importance: GB feature importances + LR coefficients
    fi = pd.DataFrame({
        "feature": fcols,
        "lr_coef": lr_final.coef_[0],
        "lr_coef_zscaled": lr_final.coef_[0],
        "gb_importance": gb_final.feature_importances_,
    }).sort_values("gb_importance", ascending=False)
    fi.to_csv(OUT / "relic_feature_importance.tsv", sep="\t", index=False)
    print("\nFeature importance:")
    print(fi.to_string(index=False), flush=True)

    # Summary
    with open(OUT / "relic_model_metrics.txt", "w") as fh:
        fh.write("Composite relic-likelihood indicator (Path B)\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Training set: {len(df)} ASVs (y=1 relic: "
                  f"{int(df['y_relic'].sum())}, y=0 alive: "
                  f"{int((1-df['y_relic']).sum())})\n")
        fh.write(f"Application:  {len(full)} ASVs scored\n\n")

        fh.write("--- Cross-validated performance ---\n\n")
        fh.write("Logistic Regression (L2, balanced):\n")
        fh.write(f"  AUC     = {auc_lr:.3f}\n")
        fh.write(f"  AP      = {ap_lr:.3f}\n")
        fh.write(f"  Brier   = {brier_lr:.3f}\n")
        fh.write(f"  LogLoss = {ll_lr:.3f}\n\n")

        fh.write("Gradient Boosting (200 trees, depth 3):\n")
        fh.write(f"  AUC     = {auc_gb:.3f}\n")
        fh.write(f"  AP      = {ap_gb:.3f}\n")
        fh.write(f"  Brier   = {brier_gb:.3f}\n")
        fh.write(f"  LogLoss = {ll_gb:.3f}\n\n")

        fh.write("--- Feature importance ---\n")
        fh.write(fi.to_string(index=False))
        fh.write("\n\n")

        fh.write("--- Calibration table (GB) ---\n")
        fh.write(cal_df.to_string(index=False))
        fh.write("\n\n")

        fh.write("--- Distribution of relic_score_gb ---\n")
        fh.write("All ASVs:\n")
        for q in (10, 25, 50, 75, 90):
            fh.write(f"  p{q}: "
                      f"{np.percentile(full['relic_score_gb'], q):.3f}\n")
        fh.write("PMA-anchored training subset:\n")
        for q in (10, 25, 50, 75, 90):
            fh.write(f"  p{q}: "
                      f"{np.percentile(sub['relic_score_gb'], q):.3f}\n")

        fh.write("\n--- HONEST CAVEATS ---\n")
        fh.write(
            "  - Training data (PMA T/UT) covers only sites C1+C2 at Trip 5,\n"
            "    rhizosphere + C2 surface compartments. The score generalizes\n"
            "    by leveraging persistence/abundance/cosmopolitanism, but its\n"
            "    *absolute* calibration (i.e., 'p=0.7 means 70% chance of\n"
            "    being relic') reflects C1+C2 conditions, not biome-wide.\n"
            "  - Use as a relic-likelihood ranking, not a binary classifier.\n"
            "  - Compartment generalization: from Test 6C QC, C1R vs C2R\n"
            "    cross-site rho was 0.34 -- moderate. The score should be\n"
            "    reported with this caveat.\n")
    print(f"\nWrote {OUT}/relic_model_metrics.txt", flush=True)


if __name__ == "__main__":
    main()
