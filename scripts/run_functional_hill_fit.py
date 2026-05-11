#!/usr/bin/env python3
"""B: Functional Hill fit using K00108 density (replaces 16S CSP1-2 abundance).

Strategy:
  - Aggregate per-sample K00108 hits and total CDS.
  - Compute K00108 density (hits per million CDS) as functional capacity.
  - Match samples to XRF S (per site×compartment).
  - Match samples to community Shannon (per sample).
  - Re-fit:
      H_baseline:  Shannon ~ a + b*S
      H_functional: Shannon ~ Hill_decline(S; V, K, n_Hill)  -- pooled
      H_functional_within_site: Shannon ~ Hill_decline within-site
      H_K00108_density: Shannon ~ Hill(K00108_density; ...) — alt mediator
  - Compare to within-EQ 16S CSP1-2 Hill (cache/within_site_hill_fit.tsv).

Outputs:
  cache/functional_hill_fit.tsv
  cache/functional_hill_fit.txt
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

# ---- Load functional census ----
samp = pd.read_csv(CACHE / "betA_per_sample_summary.tsv", sep="\t")
print(f"Per-sample census: {len(samp)} samples")
print(f"  median K00108/sample: {samp['n_K00108_strict'].median()}")
print(f"  median K00108_density: {samp['K00108_density_per_M_CDS'].median():.2f}")

# Parse trip/site from sample_id (e.g., 10PRr2 -> site 10, rhiz, rep 2)
import re
def parse_sid(sid):
    s = sid.split("_")[-1] if "_" in sid else sid
    m = re.match(r"^([0-9]+)([A-Z]+)r([0-9]+)$", s)
    if not m: return None, None, None
    return int(m.group(1)), {"PR":"rhizosphere","S":"surface","D":"deep"}[m.group(2)], int(m.group(3))

samp[["site","compartment","replicate"]] = samp["sample_id"].apply(
    lambda x: pd.Series(parse_sid(x)))

# ---- Join with XRF S (per site × compartment) and per-sample Shannon ----
xrf = pd.read_csv(CACHE/"xrf_site_compartment_panel.tsv", sep="\t")
print(f"\nXRF panel: {len(xrf)} (site, compartment) records")
shannon = pd.read_csv(CACHE/"per_sample_shannon.tsv", sep="\t")
shannon["mg_token"] = shannon["sample"].str.split("_").str[1]
print(f"Per-sample Shannon: {len(shannon)} samples")

# Match metagenomic sample_id to 16S sample via mg_token
samp_x = samp.merge(shannon[["mg_token","shannon"]], left_on="sample_id",
                     right_on="mg_token", how="left")

# Also bring in XRF S per (site, compartment) — site as int
samp_x = samp_x.merge(xrf[["site","compartment","S"]],
                      on=["site","compartment"], how="left")
merged = samp_x.dropna(subset=["S","shannon","K00108_density_per_M_CDS"])
print(f"Merged samples with S, Shannon, K00108_density: {len(merged)}")
print(f"Sites covered: {merged['site'].nunique()}")
print(f"Compartments: {merged['compartment'].value_counts().to_dict()}")

# Spearman: K00108 density vs S
rho_S, p_S = spearmanr(merged["K00108_density_per_M_CDS"], merged["S"])
print(f"\nSpearman K00108_density vs S: rho={rho_S:+.3f} (p={p_S:.3g})")
rho_Sh, p_Sh = spearmanr(merged["K00108_density_per_M_CDS"], merged["shannon"])
print(f"Spearman K00108_density vs Shannon: rho={rho_Sh:+.3f} (p={p_Sh:.3g})")
rho_lin_S, p_lin_S = spearmanr(merged["S"], merged["shannon"])
print(f"Spearman S vs Shannon: rho={rho_lin_S:+.3f} (p={p_lin_S:.3g})")

# ---- Hill fit on functional capacity vs salinity ----
def hill_decline(S, Vmax, K, n):
    return Vmax * (K**n) / (K**n + S**n)

m = merged.copy()
m["S_clip"] = m["S"].clip(lower=1e-3)

# Pre-registered EQ within-site (from within_site_hill_fit.tsv)
EQ_K_PCT = 0.02
EQ_N_HILL = 9.4

# H1: Free Hill on Shannon ~ S
try:
    popt, _ = curve_fit(hill_decline, m["S_clip"], m["shannon"],
                        p0=[m["shannon"].max(), m["S_clip"].median(), 4.0],
                        maxfev=10000,
                        bounds=([0, 1e-4, 0.5], [10, 5, 30]))
    Vmax, K, n_h = popt
    pred = hill_decline(m["S_clip"], *popt)
    sse = np.sum((m["shannon"] - pred)**2)
    sst = np.sum((m["shannon"] - m["shannon"].mean())**2)
    r2 = 1 - sse/sst
    aic = len(m)*np.log(sse/len(m)) + 2*3
except Exception as e:
    print(f"Free Hill failed: {e}")
    Vmax, K, n_h, r2, aic = (np.nan,)*5

print(f"\n=== Pooled functional Hill (Shannon ~ S) ===")
print(f"  V={Vmax:.2f}, K={K:.4f}% S, n_Hill={n_h:.2f}, R²={r2:.3f}, AIC={aic:.2f}")
print(f"  Compare to within-EQ 16S Hill: K={EQ_K_PCT}%, n={EQ_N_HILL}")
print(f"  RATIO K (functional/16S): {K/EQ_K_PCT:.2f}")
print(f"  RATIO n (functional/16S): {n_h/EQ_N_HILL:.2f}")

# H2: Within-site Hill (site fixed effect via demean)
# For each site, demean S and Shannon, then fit Hill on demeaned-positive only
m_within = m.copy()
m_within["S_dem"] = m_within.groupby(["site","compartment"])["S"].transform(lambda x: x - x.mean())
m_within["Sh_dem"] = m_within.groupby(["site","compartment"])["shannon"].transform(lambda x: x - x.mean())
# Within-site only useful if we have multiple replicates per (site, compartment)
n_replicated = m_within.groupby(["site","compartment"]).size().gt(1).sum()
print(f"\nReplicated (site,compartment) cells: {n_replicated} of {m_within.groupby(['site','compartment']).ngroups}")

# Within-site OLS (linear): beta_S
import statsmodels.api as sm
m_w = m_within.copy()
m_w["sc"] = m_w["site"].astype(str) + "_" + m_w["compartment"].astype(str)
m_w = m_w.dropna(subset=["shannon","S"])
m_w_dummies = pd.get_dummies(m_w["sc"], drop_first=True).astype(float)
X = sm.add_constant(pd.concat([m_w[["S"]], m_w_dummies], axis=1))
try:
    ols = sm.OLS(m_w["shannon"], X).fit()
    beta_S = ols.params["S"]
    se_S = ols.bse["S"]
    t_S = ols.tvalues["S"]
    p_S_lin = ols.pvalues["S"]
    print(f"Within-site linear: β_S={beta_S:+.3f} (SE={se_S:.3f}, t={t_S:.2f}, p={p_S_lin:.3g})")
except Exception as e:
    beta_S, se_S, t_S, p_S_lin = (np.nan,)*4
    print(f"OLS failed: {e}")

# H3: Functional density as alt mediator
m_alt = m.copy()
m_alt["dens_clip"] = m_alt["K00108_density_per_M_CDS"].clip(lower=1e-3)
try:
    popt2, _ = curve_fit(
        lambda x, V, K2, n2: V*(x**n2)/(K2**n2 + x**n2),  # increasing Hill
        m_alt["dens_clip"], m_alt["shannon"],
        p0=[m_alt["shannon"].max(), m_alt["dens_clip"].median(), 2.0],
        maxfev=10000, bounds=([0,1e-4,0.5],[10,1e3,30]))
    Vmax_d, K_d, n_d = popt2
    pred_d = popt2[0]*(m_alt["dens_clip"]**popt2[2])/(popt2[1]**popt2[2]+m_alt["dens_clip"]**popt2[2])
    sse_d = np.sum((m_alt["shannon"]-pred_d)**2)
    r2_d = 1 - sse_d/np.sum((m_alt["shannon"]-m_alt["shannon"].mean())**2)
    aic_d = len(m_alt)*np.log(sse_d/len(m_alt)) + 2*3
    print(f"\n=== Hill: Shannon ~ K00108_density (positive, increasing) ===")
    print(f"  V={Vmax_d:.2f}, K_density={K_d:.2f}, n={n_d:.2f}, R²={r2_d:.3f}, AIC={aic_d:.2f}")
except Exception as e:
    Vmax_d, K_d, n_d, r2_d, aic_d = (np.nan,)*5
    print(f"Density Hill failed: {e}")

# Save results
out = pd.DataFrame([
    {"specification":"functional_pooled (S→Shannon)",
     "K":K, "n_Hill":n_h, "Vmax":Vmax,
     "R2":r2, "AIC":aic, "n_obs":len(m), "comparator_K":EQ_K_PCT, "comparator_n":EQ_N_HILL},
    {"specification":"functional_within_site_linear",
     "K":np.nan, "n_Hill":np.nan, "Vmax":np.nan,
     "R2":np.nan, "AIC":np.nan, "n_obs":len(m_w),
     "comparator_K":np.nan, "comparator_n":np.nan,
     "linear_beta_S":beta_S, "linear_se_S":se_S, "linear_t_S":t_S, "linear_p_S":p_S_lin},
    {"specification":"functional_density (K00108_dens→Shannon)",
     "K":K_d, "n_Hill":n_d, "Vmax":Vmax_d,
     "R2":r2_d, "AIC":aic_d, "n_obs":len(m_alt),
     "comparator_K":np.nan, "comparator_n":np.nan},
])
out.to_csv(CACHE/"functional_hill_fit.tsv", sep="\t", index=False)
print(f"\nWrote cache/functional_hill_fit.tsv")

with open(CACHE/"functional_hill_fit.txt","w") as fh:
    fh.write("Functional Hill fit (B): K00108 density as functional mediator\n")
    fh.write("="*70 + "\n\n")
    fh.write(f"n samples merged: {len(merged)}\n")
    fh.write(f"sites: {merged['site'].nunique()}\n")
    fh.write(f"\nSpearman K00108_density vs Shannon: rho={rho_Sh:+.3f} (p={p_Sh:.3g})\n")
    fh.write(f"Spearman K00108_density vs S: rho={rho_S:+.3f} (p={p_S:.3g})\n")
    fh.write(f"Spearman S vs Shannon: rho={rho_lin_S:+.3f} (p={p_lin_S:.3g})\n")
    fh.write(out.to_string(index=False))
    fh.write("\n\n")
    fh.write("Pre-registered comparator: within-EQ 16S CSP1-2 Hill\n")
    fh.write(f"  K=0.02% S, n_Hill=9.4 (from within_site_hill_fit.tsv)\n")
print(f"Wrote cache/functional_hill_fit.txt")
