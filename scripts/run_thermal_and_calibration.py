#!/usr/bin/env python3
"""Two analyses for the trans-biome generalisation:
  (1) Cross-desert salinity proxy calibration — maps Atacama EC and
      EQ XRF S onto a common per-desert log-z scale and reports the
      threshold S_crit where Shannon falls below the per-desert
      median, with bootstrap CIs.
  (2) Thermal performance curve — CSP1-2 prevalence vs mean annual
      temperature with binomial 95% CI, fit a logistic to formalise
      the ``thermally bounded'' claim.

Outputs:
  cache/cross_desert_salinity_calibration.tsv
  cache/thermal_performance_curve.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import beta as beta_dist

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
CD = CACHE / "crossdesert"

# Approximate mean annual temperature (°C) per desert
DESERT_MAT = {
    "Namib":         18.0,    # Hewlett 2011 Namib MAT
    "Gurbantunggut":  7.5,    # Cold-temperate Gurbantunggut
    "EmptyQuarter":  28.5,    # Rub' al-Khali (NOAA AVHRR site mean)
    "Atacama":       12.0,    # Yungay-area MAT
    "McMurdo":      -16.0,    # Dry Valleys MAT
}

# ----------------------------------------------------------------------------
# (1) Cross-desert salinity calibration — derive comparable threshold
# ----------------------------------------------------------------------------
xs = pd.read_csv(CD / "per_sample.tsv", sep="\t")
xs = xs[xs["shannon_otu97"].notna()].copy()
print(f"Cross-desert samples: {len(xs)}")

# Pull the EQ frame too
frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")

def per_desert_zscore(s):
    s = s.astype(float)
    s = s[(s > 0) & np.isfinite(s)]
    z = (np.log10(s) - np.log10(s).mean()) / np.log10(s).std()
    return z

calib = []

# EQ — XRF S
eq_S = frame.loc[frame["S"].notna() & (frame["S"] > 0), "S"]
eq_sh = frame.loc[frame["S"].notna() & (frame["S"] > 0), "shannon"]
eq_z = per_desert_zscore(eq_S)
eq_med = float(eq_sh.median())
# threshold = z at which Shannon falls to 1 SD below desert median
sh_thresh = eq_med - eq_sh.std()
above = eq_z[eq_sh.values < sh_thresh]
threshold_z_eq = float(above.median()) if len(above) else float("nan")
calib.append({
    "desert": "EmptyQuarter", "proxy": "XRF S (% dry mass)",
    "n": int(len(eq_z)),
    "median_S_proxy": float(eq_S.median()),
    "shannon_median": eq_med,
    "shannon_sd": float(eq_sh.std()),
    "z_at_low_diversity": threshold_z_eq,
})

# Atacama — EC
atac = xs[(xs["desert"] == "Atacama") & xs["electrical_conductivity"].notna()].copy()
atac["EC"] = pd.to_numeric(atac["electrical_conductivity"], errors="coerce")
atac = atac.dropna(subset=["EC"])
atac = atac[atac["EC"] > 0]
if not atac.empty:
    atac_z = per_desert_zscore(atac["EC"])
    atac_sh = atac["shannon_otu97"]
    atac_med = float(atac_sh.median())
    sh_t = atac_med - atac_sh.std()
    above = atac_z[atac_sh.values < sh_t]
    th = float(above.median()) if len(above) else float("nan")
    calib.append({
        "desert": "Atacama", "proxy": "EC (mS/cm)",
        "n": int(len(atac_z)),
        "median_S_proxy": float(atac["EC"].median()),
        "shannon_median": atac_med,
        "shannon_sd": float(atac_sh.std()),
        "z_at_low_diversity": th,
    })

calib_df = pd.DataFrame(calib)
print("\nSalinity proxy threshold (z-score where Shannon drops 1 SD):")
print(calib_df.to_string(index=False))

# Bootstrap CI on the threshold for each desert
def bootstrap_threshold(s, sh, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(s)
    out = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        sb = s.iloc[idx]; shb = sh.iloc[idx]
        z = (np.log10(sb.clip(1e-6)) - np.log10(sb.clip(1e-6)).mean()) / \
            np.log10(sb.clip(1e-6)).std()
        thr = shb.median() - shb.std()
        above = z[shb.values < thr]
        if len(above):
            out.append(float(above.median()))
    return np.array(out) if out else np.array([np.nan])

if not atac.empty:
    eq_boot = bootstrap_threshold(eq_S, eq_sh)
    atac_boot = bootstrap_threshold(atac["EC"], atac_sh)
    print(f"\nEQ      threshold z = {eq_boot.mean():+.2f} "
          f"[{np.quantile(eq_boot, 0.025):+.2f}, "
          f"{np.quantile(eq_boot, 0.975):+.2f}]")
    print(f"Atacama threshold z = {atac_boot.mean():+.2f} "
          f"[{np.quantile(atac_boot, 0.025):+.2f}, "
          f"{np.quantile(atac_boot, 0.975):+.2f}]")
    calib_df.loc[calib_df["desert"] == "EmptyQuarter",
                 "z_threshold_lo"] = float(np.quantile(eq_boot, 0.025))
    calib_df.loc[calib_df["desert"] == "EmptyQuarter",
                 "z_threshold_hi"] = float(np.quantile(eq_boot, 0.975))
    calib_df.loc[calib_df["desert"] == "Atacama",
                 "z_threshold_lo"] = float(np.quantile(atac_boot, 0.025))
    calib_df.loc[calib_df["desert"] == "Atacama",
                 "z_threshold_hi"] = float(np.quantile(atac_boot, 0.975))

calib_df.to_csv(CACHE / "cross_desert_salinity_calibration.tsv",
                sep="\t", index=False)
print(f"\nwrote {CACHE / 'cross_desert_salinity_calibration.tsv'}")

# ----------------------------------------------------------------------------
# (2) Thermal performance curve — CSP1-2 prevalence vs MAT
# ----------------------------------------------------------------------------
print("\n" + "="*60)
print("Thermal performance curve")
print("="*60)
summ = pd.read_csv(CD / "comparison_summary.tsv", sep="\t")
gurb_stats = (CD / "gurbantunggut_final_stats.txt").read_text()
# Parse Gurb prevalence
import re
m = re.search(r"CSP1-2 prevalence at 85% V4: (\d+)/(\d+) \(([\d.]+)%\)",
              gurb_stats)
gurb_present, gurb_n, gurb_pct = int(m.group(1)), int(m.group(2)), float(m.group(3))

rows = []
for d in ("Namib", "EmptyQuarter", "Atacama", "McMurdo"):
    rec = summ[summ["desert"].str.replace(" (baseline)", "", regex=False) == d]
    if rec.empty:
        rec = summ[summ["desert"].str.contains(d)]
    if rec.empty:
        continue
    n_total = int(rec.iloc[0]["n_samples"])
    f_present = float(rec.iloc[0]["frac_samples_with_CSP85"])
    n_present = int(round(f_present * n_total))
    rows.append({"desert": d, "n_total": n_total, "n_present": n_present,
                 "fraction": f_present, "MAT_C": DESERT_MAT[d]})
rows.append({"desert": "Gurbantunggut",
             "n_total": gurb_n, "n_present": gurb_present,
             "fraction": gurb_pct / 100.0,
             "MAT_C": DESERT_MAT["Gurbantunggut"]})

# Wilson CI for prevalence
from scipy.stats import beta as betadist
def wilson_ci(k, n, alpha=0.05):
    if n == 0:
        return (np.nan, np.nan)
    lo = betadist.ppf(alpha/2, k, n - k + 1) if k > 0 else 0
    hi = betadist.ppf(1 - alpha/2, k + 1, n - k) if k < n else 1
    return float(lo), float(hi)

for r in rows:
    lo, hi = wilson_ci(r["n_present"], r["n_total"])
    r["ci_lo"] = lo
    r["ci_hi"] = hi
    print(f"  {r['desert']:13s}  MAT = {r['MAT_C']:+5.1f}°C  "
          f"prev = {r['fraction']*100:5.1f}%  "
          f"CI = [{lo*100:5.1f}, {hi*100:5.1f}]  n = {r['n_total']}")

# Logistic fit: prevalence vs MAT
def logi(x, T_50, k):
    return 1.0 / (1.0 + np.exp(-k * (x - T_50)))

x = np.array([r["MAT_C"] for r in rows], dtype=float)
y = np.array([r["fraction"] for r in rows], dtype=float)
try:
    popt, pcov = curve_fit(logi, x, y, p0=[5.0, 0.3], maxfev=10000)
    se = np.sqrt(np.diag(pcov))
    print(f"\nLogistic fit: T_50 = {popt[0]:+.1f}°C ± {se[0]:.1f},  "
          f"k = {popt[1]:.2f} per °C")
    for r in rows:
        r["fit_logi_T50"] = float(popt[0])
        r["fit_logi_k"] = float(popt[1])
except Exception as e:
    print(f"logistic fit failed: {e}")

therm_df = pd.DataFrame(rows)
therm_df.to_csv(CACHE / "thermal_performance_curve.tsv",
                sep="\t", index=False)
print(f"\nwrote {CACHE / 'thermal_performance_curve.tsv'}")
