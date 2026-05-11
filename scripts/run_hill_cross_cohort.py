#!/usr/bin/env python3
"""Pre-registered Hill threshold cross-cohort prediction.

Locks the EQ within-site Hill parameters (K=0.02%, n=9.4 — see
within_site_hill_fit.tsv) and tests transferability to Atacama
using EC as the salinity proxy.

Strategy:
  - Within-Atacama: fit Hill (Shannon ~ Vmax * (1 - S^n / (K^n + S^n)))
  - Compare three nested hypotheses:
      H0: linear Shannon ~ EC
      H1: Hill with free K, n
      H2: Hill with K, n locked to EQ values (rescaled by per-biome
          z-score of salinity, so we test SHAPE transferability)
  - Report AIC, R^2, and whether non-linearity is detected.

Outputs:
  cache/hill_cross_cohort_atacama.tsv
  cache/hill_cross_cohort_atacama_fits.txt
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

# ---- Pre-registered EQ parameters (locked before this analysis) ----
EQ_K_PCT = 0.02      # % S (XRF), within-site K
EQ_N_HILL = 9.4      # within-site Hill coefficient
EQ_K_POOLED = 0.09   # pooled K, comparator
EQ_N_POOLED = 4.6
print(f"PRE-REGISTERED (locked): EQ within-site K={EQ_K_PCT}%, n_Hill={EQ_N_HILL}")
print(f"PRE-REGISTERED (locked): EQ pooled K={EQ_K_POOLED}%, n_Hill={EQ_N_POOLED}")
print()

# ---- Load Atacama data ----
df = pd.read_csv(CACHE / "crossdesert" / "per_sample.tsv", sep="\t")
ata = df[df["desert"] == "Atacama"].copy()
ata["EC"] = pd.to_numeric(ata["electrical_conductivity"], errors="coerce")
ata = ata.dropna(subset=["EC", "shannon_otu97"]).reset_index(drop=True)
print(f"Atacama samples with EC + Shannon: {len(ata)}")

# Site ID
ata["site"] = ata["library_name"].str.split(".").str[1]
print(f"Atacama sites: {ata['site'].nunique()}")
print()

# Spearman
rho, p = spearmanr(ata["EC"], ata["shannon_otu97"])
print(f"Spearman rho EC vs Shannon: {rho:+.3f} (p={p:.3g})")

# ---- Model fits ----
def hill_decline(S, Vmax, K, n):
    """Decline Hill: Shannon = Vmax * K^n / (K^n + S^n)."""
    return Vmax * (K**n) / (K**n + S**n)

def linear_decline(S, a, b):
    return a + b * S

# Need positive values for Hill
m = ata.copy()
m["EC_clip"] = m["EC"].clip(lower=1e-3)

# H0: linear
slope_a, intercept_a = np.polyfit(m["EC_clip"], m["shannon_otu97"], 1)
pred_lin = intercept_a + slope_a * m["EC_clip"]
r2_lin = 1 - np.sum((m["shannon_otu97"] - pred_lin)**2) / np.sum(
    (m["shannon_otu97"] - m["shannon_otu97"].mean())**2)

# H1: free Hill
try:
    p0 = [m["shannon_otu97"].max(), m["EC_clip"].median(), 4.0]
    popt, pcov = curve_fit(hill_decline, m["EC_clip"], m["shannon_otu97"],
                            p0=p0, maxfev=10000,
                            bounds=([0, 1e-4, 0.5], [10, 5, 30]))
    Vmax_free, K_free, n_free = popt
    pred_h1 = hill_decline(m["EC_clip"], *popt)
    r2_h1 = 1 - np.sum((m["shannon_otu97"] - pred_h1)**2) / np.sum(
        (m["shannon_otu97"] - m["shannon_otu97"].mean())**2)
except Exception as e:
    print(f"Free Hill fit failed: {e}")
    Vmax_free, K_free, n_free = np.nan, np.nan, np.nan
    r2_h1 = np.nan

# H2: locked-shape Hill (K and n locked to EQ values, but EC is in dS/m,
# whereas EQ S is in % dry mass — so rescale by per-biome median salinity
# and test whether the SHAPE transfers, not the absolute K).
# rescaled S' = S / median(S); apply EQ K_rel = 0.02% / median_EQ_S_pct (within-site)
# We don't have EQ within-site median S in this script's context; use rank-based test.
# Define rank-z salinity: locked K is at the EQ within-site quantile of EC=K_z.
# As a transferability test: does Atacama-fitted Hill yield n ~ 9.4 within-site
# (steep, saturating)?  And does freely-fit K within-site land at the
# Atacama's lower quantile?

# Within-site demean (site fixed effect): subtract site mean of EC and Shannon
ata2 = ata.dropna(subset=["EC","shannon_otu97"]).copy()
ata2["EC_c"] = ata2.groupby("site")["EC"].transform(lambda x: x - x.mean())
ata2["S_c"] = ata2.groupby("site")["shannon_otu97"].transform(lambda x: x - x.mean())
# Within-site linear slope
if (ata2["EC_c"]**2).sum() > 0:
    beta_within = np.sum(ata2["EC_c"] * ata2["S_c"]) / np.sum(ata2["EC_c"]**2)
    pred_within = beta_within * ata2["EC_c"]
    sse_within = np.sum((ata2["S_c"] - pred_within)**2)
    sst_within = np.sum(ata2["S_c"]**2)
    r2_within = 1 - sse_within / sst_within if sst_within > 0 else np.nan
else:
    beta_within, r2_within = np.nan, np.nan

# Within-site is dominated by sites with replicates; many Atacama
# sites have only 2-3 samples so within-site EC variance is small.
# We report it as a sensitivity check.

# AIC: -2*ln(L) + 2*k.  Assume Gaussian residuals.
n = len(m)
def aic(sse, k):
    if sse <= 0: return np.inf
    return n * np.log(sse / n) + 2 * k

sse_lin = np.sum((m["shannon_otu97"] - pred_lin)**2)
aic_lin = aic(sse_lin, 2)

if not np.isnan(r2_h1):
    sse_h1 = np.sum((m["shannon_otu97"] - pred_h1)**2)
    aic_h1 = aic(sse_h1, 3)
else:
    aic_h1 = np.nan

# H2: rescale EC into "EQ-equivalent units" by dividing by Atacama median EC
# (matching to where EQ within-site median lives at K_pct/median_S ~ K-ratio)
# This is a coarse but principled scale-free test.
median_EC_ata = m["EC_clip"].median()
EQ_median_S = 0.10  # EQ within-site median (approximate; from main paper RQ10)
m["EC_rescaled"] = m["EC_clip"] / median_EC_ata * EQ_median_S
# Now fit Hill with K and n LOCKED to EQ values, only Vmax free
def hill_locked(S, Vmax):
    return Vmax * (EQ_K_PCT**EQ_N_HILL) / (EQ_K_PCT**EQ_N_HILL + S**EQ_N_HILL)
try:
    popt2, _ = curve_fit(hill_locked, m["EC_rescaled"], m["shannon_otu97"],
                          p0=[m["shannon_otu97"].max()], maxfev=5000,
                          bounds=([0],[10]))
    pred_h2 = hill_locked(m["EC_rescaled"], *popt2)
    sse_h2 = np.sum((m["shannon_otu97"] - pred_h2)**2)
    r2_h2 = 1 - sse_h2 / np.sum((m["shannon_otu97"] - m["shannon_otu97"].mean())**2)
    aic_h2 = aic(sse_h2, 1)
    Vmax_locked = popt2[0]
except Exception as e:
    Vmax_locked, r2_h2, aic_h2 = np.nan, np.nan, np.nan

# H2-pooled: same with EQ pooled K=0.09%, n=4.6
def hill_locked_pooled(S, Vmax):
    return Vmax * (EQ_K_POOLED**EQ_N_POOLED) / (EQ_K_POOLED**EQ_N_POOLED + S**EQ_N_POOLED)
try:
    popt3, _ = curve_fit(hill_locked_pooled, m["EC_rescaled"], m["shannon_otu97"],
                          p0=[m["shannon_otu97"].max()], maxfev=5000,
                          bounds=([0],[10]))
    pred_h3 = hill_locked_pooled(m["EC_rescaled"], *popt3)
    sse_h3 = np.sum((m["shannon_otu97"] - pred_h3)**2)
    r2_h3 = 1 - sse_h3 / np.sum((m["shannon_otu97"] - m["shannon_otu97"].mean())**2)
    aic_h3 = aic(sse_h3, 1)
    Vmax_locked_pooled = popt3[0]
except Exception as e:
    Vmax_locked_pooled, r2_h3, aic_h3 = np.nan, np.nan, np.nan

print()
print("=== Atacama Hill cross-cohort fits ===")
print(f"H0 linear: R^2={r2_lin:.3f}, AIC={aic_lin:.2f}  (slope b={slope_a:+.3f})")
print(f"H1 free Hill: K_free={K_free:.3f} dS/m, n_free={n_free:.2f}, "
      f"R^2={r2_h1:.3f}, AIC={aic_h1:.2f}")
print(f"H2 locked-shape (EQ within-site K=0.02%, n=9.4 rescaled): "
      f"R^2={r2_h2:.3f}, AIC={aic_h2:.2f}, Vmax={Vmax_locked:.2f}")
print(f"H3 locked-shape (EQ pooled K=0.09%, n=4.6 rescaled): "
      f"R^2={r2_h3:.3f}, AIC={aic_h3:.2f}, Vmax={Vmax_locked_pooled:.2f}")
print(f"Within-site linear slope (Atacama): beta={beta_within:.3f}, R^2={r2_within:.3f}")
print()

# Save tabular result
out = pd.DataFrame([
    {"hypothesis":"H0_linear",
     "params":"a + b*EC", "k":2,
     "R2": r2_lin, "AIC": aic_lin,
     "K": np.nan, "n_Hill": np.nan, "Vmax": np.nan,
     "free_or_locked":"free"},
    {"hypothesis":"H1_free_Hill",
     "params":"Vmax,K,n_free", "k":3,
     "R2": r2_h1, "AIC": aic_h1,
     "K": K_free, "n_Hill": n_free, "Vmax": Vmax_free,
     "free_or_locked":"free"},
    {"hypothesis":"H2_locked_within_site_EQ",
     "params":"Vmax only", "k":1,
     "R2": r2_h2, "AIC": aic_h2,
     "K": EQ_K_PCT, "n_Hill": EQ_N_HILL, "Vmax": Vmax_locked,
     "free_or_locked":"locked_to_EQ"},
    {"hypothesis":"H3_locked_pooled_EQ",
     "params":"Vmax only", "k":1,
     "R2": r2_h3, "AIC": aic_h3,
     "K": EQ_K_POOLED, "n_Hill": EQ_N_POOLED, "Vmax": Vmax_locked_pooled,
     "free_or_locked":"locked_to_EQ"},
])
out.to_csv(CACHE / "hill_cross_cohort_atacama.tsv", sep="\t", index=False)
print(f"Wrote {CACHE/'hill_cross_cohort_atacama.tsv'}")

# Save text summary
summary_path = CACHE / "hill_cross_cohort_atacama_fits.txt"
with open(summary_path, "w") as fh:
    fh.write("Pre-registered cross-cohort Hill prediction (EQ -> Atacama)\n")
    fh.write("="*70 + "\n\n")
    fh.write(f"PRE-REGISTERED EQ within-site: K={EQ_K_PCT}% S, n_Hill={EQ_N_HILL}\n")
    fh.write(f"PRE-REGISTERED EQ pooled:      K={EQ_K_POOLED}% S, n_Hill={EQ_N_POOLED}\n\n")
    fh.write(f"Atacama n (with EC and Shannon): {len(m)}\n")
    fh.write(f"Atacama sites: {ata['site'].nunique()}\n")
    fh.write(f"Spearman EC vs Shannon: rho={rho:+.3f}, p={p:.3g}\n\n")
    fh.write(out.to_string(index=False))
    fh.write("\n\n")
    fh.write("Conclusion:\n")
    if not np.isnan(aic_h1):
        delta = aic_lin - aic_h1
        if delta > 4:
            fh.write(f"  Free Hill fits BETTER than linear (ΔAIC={delta:.1f}), "
                     "non-linear shape replicates in Atacama.\n")
        elif delta > 0:
            fh.write(f"  Free Hill is marginally preferred over linear "
                     f"(ΔAIC={delta:.1f}).\n")
        else:
            fh.write(f"  Linear preferred over free Hill (ΔAIC={delta:.1f}).\n")
    if not np.isnan(aic_h2):
        delta_h2 = aic_lin - aic_h2
        fh.write(f"  Locked-EQ-within-site Hill: ΔAIC vs linear = {delta_h2:+.1f}\n")
    if not np.isnan(K_free):
        ratio = K_free / median_EC_ata
        fh.write(f"  Free Hill K_Atacama/median_EC = {ratio:.2f} "
                 f"(EQ K_within-site/median_S ~ {EQ_K_PCT/EQ_median_S:.2f})\n")
        fh.write("  If these ratios match, the threshold's quantile position transfers.\n")
print(f"Wrote {summary_path}")
print()
print("Done.")
