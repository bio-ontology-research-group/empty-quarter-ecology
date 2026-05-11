#!/usr/bin/env python3
"""Decompose the salinity → CSP1-2 Hill threshold into within-site
and between-site components.

If the threshold is identifiable *within* sites (e.g., across
compartments or across trips at a fixed site), it is a much stronger
mechanistic claim than if it is only between-site (where unmeasured
site-level confounders could explain it).

Approach:
  1. Demean S and CSP1-2 by site → within-site variation.
  2. Refit Hill to within-site demeaned data.
  3. Compare fits and report the within-site n_Hill, K, and slope.

Output:
  cache/within_site_hill_fit.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")

# Cell-aggregated panel
cells = (frame.dropna(subset=["S", "csp_relab"])
         .groupby(["trip", "site", "compartment"])
         .agg({"S": "mean", "csp_relab": "mean"})
         .reset_index())
print(f"Cells: n = {len(cells)}, sites = {cells['site'].nunique()}")

# ----------------------------------------------------------------------------
# (1) Pooled Hill fit
# ----------------------------------------------------------------------------
def hill(s, M0, V, K, n):
    return M0 - V * s ** n / (K ** n + s ** n)

def fit_hill(s, m, p0=(0.005, 0.003, 0.1, 3.0)):
    try:
        popt, _ = curve_fit(hill, s, m, p0=p0, maxfev=10000,
                             bounds=([0, 0, 1e-4, 0.5],
                                     [1, 1, 10, 10]))
        return tuple(popt)
    except Exception as e:
        return None

s_all = cells["S"].values
m_all = cells["csp_relab"].values
mask = (s_all > 0) & np.isfinite(s_all) & np.isfinite(m_all)
fit_pool = fit_hill(s_all[mask], m_all[mask])
print(f"\nPooled fit:  M0={fit_pool[0]:.4f} V={fit_pool[1]:.4f} "
      f"K={fit_pool[2]:.4f} n_Hill={fit_pool[3]:.2f}")

# ----------------------------------------------------------------------------
# (2) Within-site demeaning
# ----------------------------------------------------------------------------
mean_S = cells.groupby("site")["S"].transform("mean")
mean_M = cells.groupby("site")["csp_relab"].transform("mean")
cells["S_within"] = cells["S"] - mean_S
cells["M_within"] = cells["csp_relab"] - mean_M
cells["S_between"] = mean_S
cells["M_between"] = mean_M

# Need positive S to fit Hill; demeaned S has negative values, so add
# back the global S median to recover absolute scale while removing
# site-level variation.
S_med = float(cells["S"].median())
cells["S_resid_pos"] = cells["S_within"] + S_med
cells["M_resid_pos"] = cells["M_within"] + cells["csp_relab"].median()

mask_w = (cells["S_resid_pos"] > 0) & np.isfinite(cells["S_resid_pos"]) \
         & np.isfinite(cells["M_resid_pos"])
fit_within = fit_hill(cells.loc[mask_w, "S_resid_pos"].values,
                       cells.loc[mask_w, "M_resid_pos"].values)
if fit_within:
    print(f"Within-site fit: M0={fit_within[0]:.4f} V={fit_within[1]:.4f} "
          f"K={fit_within[2]:.4f} n_Hill={fit_within[3]:.2f}")
else:
    print("Within-site fit: failed to converge")

# ----------------------------------------------------------------------------
# (3) Between-site fit (per-site means)
# ----------------------------------------------------------------------------
site_means = cells.groupby("site").agg(
    {"S": "mean", "csp_relab": "mean"}).reset_index()
mask_b = (site_means["S"] > 0) & np.isfinite(site_means["csp_relab"])
fit_between = fit_hill(site_means.loc[mask_b, "S"].values,
                        site_means.loc[mask_b, "csp_relab"].values)
print(f"Between-site fit: M0={fit_between[0]:.4f} V={fit_between[1]:.4f} "
      f"K={fit_between[2]:.4f} n_Hill={fit_between[3]:.2f}")

# ----------------------------------------------------------------------------
# (4) Linear within-site coefficient (S → CSP1-2 with site fixed effects)
# ----------------------------------------------------------------------------
import statsmodels.formula.api as smf
cells_fe = cells.copy()
cells_fe["site"] = cells_fe["site"].astype(int).astype(str)
cells_fe["trip"] = cells_fe["trip"].astype(int).astype(str)
cells_fe["compartment"] = cells_fe["compartment"].astype(str)
m = smf.ols("csp_relab ~ S + C(site) + C(trip) + C(compartment)",
            data=cells_fe).fit(cov_type="cluster",
                                cov_kwds={"groups": cells_fe["site"].values})
print(f"\nLinear S → CSP1-2 with site FE: "
      f"β_S = {m.params['S']:+.4f}, "
      f"clustered SE = {m.bse['S']:.4f}, "
      f"t = {m.tvalues['S']:+.2f}, "
      f"p = {m.pvalues['S']:.3g}")

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
out = pd.DataFrame([
    {"specification": "pooled",
     "M0": fit_pool[0], "V": fit_pool[1], "K": fit_pool[2],
     "n_Hill": fit_pool[3], "n_obs": int(mask.sum())},
    {"specification": "within-site (demeaned)",
     "M0": fit_within[0] if fit_within else np.nan,
     "V":  fit_within[1] if fit_within else np.nan,
     "K":  fit_within[2] if fit_within else np.nan,
     "n_Hill": fit_within[3] if fit_within else np.nan,
     "n_obs": int(mask_w.sum())},
    {"specification": "between-site (site means)",
     "M0": fit_between[0], "V": fit_between[1],
     "K":  fit_between[2], "n_Hill": fit_between[3],
     "n_obs": int(mask_b.sum())},
    {"specification": "linear S→CSP, site+trip+compartment FE",
     "M0": np.nan, "V": np.nan, "K": np.nan,
     "n_Hill": np.nan, "n_obs": int(len(cells)),
     "linear_beta_S": float(m.params["S"]),
     "linear_se_S": float(m.bse["S"]),
     "linear_t_S": float(m.tvalues["S"]),
     "linear_p_S": float(m.pvalues["S"])},
])
out.to_csv(CACHE / "within_site_hill_fit.tsv", sep="\t", index=False)
print(f"\nwrote {CACHE / 'within_site_hill_fit.tsv'}")
