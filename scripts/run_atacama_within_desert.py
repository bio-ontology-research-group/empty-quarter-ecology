#!/usr/bin/env python3
"""Atacama re-analysis: identify positive variables that explain CSP1-2 presence.

Within-desert logistic regression on Atacama 53 samples:
  P(CSP1-2 detected) ~ elevation + soil_RH + soil_T + EC + pH

Outputs:
  cache/atacama_within_desert_logit.tsv
  cache/atacama_within_desert_summary.txt
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

df = pd.read_csv(CACHE / "crossdesert" / "per_sample.tsv", sep="\t")
print(f"Total cross-desert samples: {len(df)}")

ata = df[df["desert"] == "Atacama"].copy()
print(f"Atacama samples: {len(ata)}")

ata["has_csp"] = (ata["csp_rel_85"] > 0).astype(int)
print(f"CSP+ samples (85% id): {ata['has_csp'].sum()} / {len(ata)} "
      f"({ata['has_csp'].mean():.1%})")

# Numeric coercion
for col in ["elevation_m","avg_soil_rh","avg_soil_temp",
            "electrical_conductivity","ph","soil_organic_carbon"]:
    ata[col] = pd.to_numeric(ata[col], errors="coerce")

# Quick descriptive: CSP+ vs CSP-
print("\n=== Descriptive: CSP+ vs CSP- ===")
desc = ata.groupby("has_csp").agg(
    n=("run","count"),
    elevation_med=("elevation_m","median"),
    elevation_iqr=("elevation_m", lambda x: x.quantile(0.75) - x.quantile(0.25)),
    RH_med=("avg_soil_rh","median"),
    T_med=("avg_soil_temp","median"),
    EC_med=("electrical_conductivity","median"),
    pH_med=("ph","median"),
    SOC_med=("soil_organic_carbon","median"),
)
print(desc.to_string())

# Mann-Whitney U for each variable
from scipy.stats import mannwhitneyu

results = []
for col in ["elevation_m","avg_soil_rh","avg_soil_temp",
            "electrical_conductivity","ph","soil_organic_carbon"]:
    a = ata[ata["has_csp"]==1][col].dropna().values
    b = ata[ata["has_csp"]==0][col].dropna().values
    if len(a) >= 3 and len(b) >= 3:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        results.append({
            "variable": col,
            "n_csp_pos": len(a),
            "n_csp_neg": len(b),
            "median_pos": np.median(a),
            "median_neg": np.median(b),
            "mannwhitney_p": p,
        })
        print(f"  {col}: median+ = {np.median(a):.2f}, median- = {np.median(b):.2f}, p = {p:.3g}")

res_df = pd.DataFrame(results)
res_df.to_csv(CACHE / "atacama_within_desert_logit.tsv", sep="\t", index=False)

# Logistic regression on most informative variables (drop nas)
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Standardise predictors for interpretable coefficients
ata2 = ata.dropna(subset=["elevation_m","avg_soil_rh","avg_soil_temp"]).copy()
print(f"\nLogistic regression sample size: {len(ata2)}")

for col in ["elevation_m","avg_soil_rh","avg_soil_temp"]:
    if col in ata2.columns:
        m, s = ata2[col].mean(), ata2[col].std()
        ata2[col + "_z"] = (ata2[col] - m) / s

try:
    mod = smf.logit("has_csp ~ elevation_m_z + avg_soil_rh_z + avg_soil_temp_z",
                    data=ata2).fit(disp=False, maxiter=200)
    print("\n=== Logistic regression: P(CSP+) ~ z(elev) + z(RH) + z(T) ===")
    print(mod.summary())
    summary_path = CACHE / "atacama_within_desert_summary.txt"
    with open(summary_path, "w") as fh:
        fh.write(str(mod.summary()))
        fh.write("\n\nDescriptive (CSP+ vs CSP-):\n")
        fh.write(desc.to_string())
        fh.write("\n\nMann-Whitney p-values:\n")
        fh.write(res_df.to_string(index=False))
        fh.write("\n")
    print(f"\nWrote {summary_path}")
except Exception as e:
    print(f"Logit failed: {e}")

# Cross-tabulate by elevation tertile and RH tertile
print("\n=== Atacama prevalence by elevation tertile ===")
ata["elev_t"] = pd.qcut(ata["elevation_m"], q=3,
                         labels=["low","mid","high"], duplicates="drop")
print(ata.groupby("elev_t", observed=True).agg(
    n=("run","count"),
    csp_prev=("has_csp","mean"),
).to_string())

print("\n=== Atacama prevalence by RH tertile (where available) ===")
ata_rh = ata.dropna(subset=["avg_soil_rh"]).copy()
ata_rh["RH_t"] = pd.qcut(ata_rh["avg_soil_rh"], q=3,
                          labels=["low","mid","high"], duplicates="drop")
print(ata_rh.groupby("RH_t", observed=True).agg(
    n=("run","count"),
    csp_prev=("has_csp","mean"),
).to_string())

# 2x2 contingency: high elev (>3000m) vs CSP detection
ata["high_elev"] = (ata["elevation_m"] >= 3000).astype(int)
ct = pd.crosstab(ata["high_elev"], ata["has_csp"])
print(f"\n=== 2x2 contingency: elev >= 3000 m vs CSP+ ===")
print(ct.to_string())

from scipy.stats import fisher_exact
if ct.shape == (2,2):
    odds, p_fisher = fisher_exact(ct.values)
    print(f"Fisher's exact: OR = {odds:.2f}, p = {p_fisher:.3g}")

print("\nDone.")
