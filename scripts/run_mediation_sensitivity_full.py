#!/usr/bin/env python3
"""Causal sensitivity analysis for the salinity → Shannon mediation.

Three complementary tests:
  (1) E-values  — how strong must unmeasured confounding be (on the
                  RR scale) to nullify the indirect or direct effect.
  (2) Imai ρ    — residual M–Y correlation that drives the indirect
                  effect to zero (sequential ignorability test).
  (3) Covariate-drop stability — re-fit dropping each covariate
                  group (climate, depth, host plant) one at a time.

Mediator alternatives: csp_relab (CSP1-2) and dep_pool (11-genus
guild). Bootstrap B = 1,000.

Output:
  cache/causal_tier1_mediation_sensitivity_full.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
sys.path.insert(0, str(REPO / "src"))

# ----------------------------------------------------------------------------
# Load causal frame
# ----------------------------------------------------------------------------
frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")
ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")

# Add the dependent-pool mediator
if "genus" not in tax.columns and "Genus" in tax.columns:
    tax = tax.rename(columns={"Genus": "genus"})
ft_rel = ft.div(ft.sum(axis=0), axis=1)
asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")
genus_relab = ft_rel.assign(genus=asv_to_genus.values).groupby("genus").sum()
dep_genera = ['Herpetosiphon', 'Paenibacillus', 'Flavisolibacter',
              'Ammoniphilus', 'Streptomyces', 'Rubrobacter', 'Ectobacillus',
              'Neobacillus', 'Ramlibacter', 'Noviherbaspirillum', 'Nocardioides']
dep_present = [g for g in dep_genera if g in genus_relab.index]
dep_pool = genus_relab.loc[dep_present].sum(axis=0)
frame = frame.set_index("sample")
frame["dep_pool"] = dep_pool.reindex(frame.index).fillna(0)
frame = frame.reset_index()

# ----------------------------------------------------------------------------
# Cell aggregation (variant C from the robustness panel) — n_cells = 622
# ----------------------------------------------------------------------------
ELEMENTS = ["S", "Cl", "Na", "P"]
CLIMATE = ["rain_W30d", "temp_mean_W30d"]
COMPARTMENT = ["compartment"]
agg_cols = ELEMENTS + CLIMATE + ["shannon", "csp_relab", "dep_pool"]

# Cell groups
cells = (frame.dropna(subset=agg_cols)
         .groupby(["trip", "site", "compartment"])
         .agg({c: "mean" for c in agg_cols} | {"compartment": "first"})
         .reset_index(drop=True))
cells = cells.dropna()
print(f"Cells: n = {len(cells)}")

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def fit_mediation(d, mediator, covars, nboot=1000, seed=0):
    """Bootstrap a/b path mediation. Returns dict of effect estimates."""
    rng = np.random.default_rng(seed)
    n = len(d)
    a_b, b_b, c_b, cp_b = [], [], [], []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        db = d.iloc[idx]
        # a path: S → M | covars
        X = db[covars].values
        try:
            a = LinearRegression().fit(
                np.c_[X, db["S"].values], db[mediator].values).coef_[-1]
            # b path: M → Y | S, covars
            b = LinearRegression().fit(
                np.c_[X, db["S"].values, db[mediator].values],
                db["shannon"].values).coef_[-1]
            # c path: S → Y total
            c = LinearRegression().fit(
                np.c_[X, db["S"].values], db["shannon"].values).coef_[-1]
            # c': S → Y | M
            cp = LinearRegression().fit(
                np.c_[X, db[mediator].values, db["S"].values],
                db["shannon"].values).coef_[-1]
        except Exception:
            continue
        a_b.append(a); b_b.append(b); c_b.append(c); cp_b.append(cp)
    a_b, b_b, c_b, cp_b = map(np.array, (a_b, b_b, c_b, cp_b))
    indirect = a_b * b_b
    return {
        "n": int(n), "mediator": mediator,
        "covars": ",".join(covars),
        "a_mean": float(a_b.mean()),
        "b_mean": float(b_b.mean()),
        "indirect": float(indirect.mean()),
        "indirect_lo": float(np.quantile(indirect, 0.025)),
        "indirect_hi": float(np.quantile(indirect, 0.975)),
        "direct": float(cp_b.mean()),
        "direct_lo": float(np.quantile(cp_b, 0.025)),
        "direct_hi": float(np.quantile(cp_b, 0.975)),
        "total": float(c_b.mean()),
    }

# ----------------------------------------------------------------------------
# (3) Covariate-drop stability
# ----------------------------------------------------------------------------
# Base = climate + P (a confounder we always include via XRF)
BASE_COVARS = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d"]
DROP_TESTS = [
    ("baseline (full)",   BASE_COVARS),
    ("- climate",         ["Cl", "Na", "P"]),
    ("- P (nutrient)",    ["Cl", "Na", "rain_W30d", "temp_mean_W30d"]),
    ("- correlated salts (Cl, Na)", ["P", "rain_W30d", "temp_mean_W30d"]),
    ("- all controls",    []),
]

rows = []
for label, covars in DROP_TESTS:
    for med in ["csp_relab", "dep_pool"]:
        out = fit_mediation(cells, med, covars, nboot=1000, seed=42)
        out["covariate_set"] = label
        rows.append(out)
        print(f"  {label:30s}  med={med:11s}  "
              f"indirect = {out['indirect']:+.3f} "
              f"[{out['indirect_lo']:+.3f}, {out['indirect_hi']:+.3f}]  "
              f"direct = {out['direct']:+.3f}")

# ----------------------------------------------------------------------------
# (1) E-value (VanderWeele & Ding 2017) for the indirect effect
# ----------------------------------------------------------------------------
# Convert standardised continuous effect to approximate RR per +1 SD of S
# RR_approx = exp(0.91 * d_cohen) for an outcome dichotomised at the median
def evalue_continuous(beta, sd_y, sd_x):
    """E-value approximation for continuous outcome (Linden et al. 2020)."""
    d = abs(beta) * sd_x / sd_y
    rr = np.exp(0.91 * d)
    e = rr + np.sqrt(rr * (rr - 1))
    return float(rr), float(e)

base_csp = next(r for r in rows
                if r["covariate_set"] == "baseline (full)"
                and r["mediator"] == "csp_relab")
base_dep = next(r for r in rows
                if r["covariate_set"] == "baseline (full)"
                and r["mediator"] == "dep_pool")

sd_y = cells["shannon"].std()
sd_S = cells["S"].std()
print("\nE-values (Linden 2020; how strong must unmeasured confounding be):")
for r in (base_csp, base_dep):
    rr_ind, e_ind = evalue_continuous(r["indirect"], sd_y, sd_S)
    rr_dir, e_dir = evalue_continuous(r["direct"], sd_y, sd_S)
    print(f"  med = {r['mediator']:11s}  "
          f"indirect: RR≈{rr_ind:.2f} → E={e_ind:.2f}   "
          f"direct: RR≈{rr_dir:.2f} → E={e_dir:.2f}")
    r["E_indirect"] = e_ind
    r["E_direct"] = e_dir
    r["RR_indirect"] = rr_ind
    r["RR_direct"] = rr_dir

# ----------------------------------------------------------------------------
# (2) Imai ρ — residual M–Y correlation that nullifies the indirect effect
# ----------------------------------------------------------------------------
# After fitting the mediation regressions, compute residuals from
#   M | S, covars    and    Y | S, M, covars
# ρ_obs is the correlation between these residuals (≈ 0 under sequential
# ignorability). The bound on |ρ| at which indirect effect = 0 is found
# numerically by scanning the parameter (Imai et al. 2010).
def imai_rho(d, mediator, covars):
    X = d[covars + ["S"]].values
    M = d[mediator].values
    # M | S, covars
    res_M = M - LinearRegression().fit(X, M).predict(X)
    Y = d["shannon"].values
    XM = np.c_[X, M]
    res_Y = Y - LinearRegression().fit(XM, Y).predict(XM)
    rho_obs = float(np.corrcoef(res_M, res_Y)[0, 1])
    # Bound: the value of ρ that would make indirect = 0
    # Approximation: for centred residuals,
    #   indirect_adjusted = indirect - ρ * sd(M) * sd(Y) * sd(S)/var(S)
    # We solve numerically.
    sd_M = res_M.std()
    sd_Y_res = res_Y.std()
    sd_S_res = X[:, -1].std()
    a_ = LinearRegression().fit(X[:, :-1], M).predict(X[:, :-1])
    a_resid = M - a_
    a_path = np.cov(a_resid, X[:, -1])[0, 1] / np.var(X[:, -1])
    b_ = LinearRegression().fit(XM, Y).coef_[-1]
    indirect = a_path * b_
    rhos = np.linspace(-0.99, 0.99, 1001)
    ind_adj = indirect - rhos * sd_M * sd_Y_res / np.var(X[:, -1])
    # ρ at which |ind_adj| crosses zero
    sign_change = np.where(np.diff(np.sign(ind_adj)))[0]
    rho_zero = float(rhos[sign_change[0]]) if len(sign_change) else float("nan")
    return rho_obs, rho_zero, float(indirect)

print("\nImai sensitivity (ρ value that nullifies indirect effect):")
for med in ("csp_relab", "dep_pool"):
    rho_obs, rho_zero, ind = imai_rho(cells, med, BASE_COVARS)
    print(f"  med = {med:11s}  ρ_observed = {rho_obs:+.3f}   "
          f"ρ_at_zero = {rho_zero:+.3f}   indirect = {ind:+.3f}")
    # Annotate base rows
    for r in rows:
        if r["mediator"] == med and r["covariate_set"] == "baseline (full)":
            r["rho_obs"] = rho_obs
            r["rho_zero"] = rho_zero

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
out = pd.DataFrame(rows)
out.to_csv(CACHE / "causal_tier1_mediation_sensitivity_full.tsv",
           sep="\t", index=False)
print(f"\nwrote {CACHE / 'causal_tier1_mediation_sensitivity_full.tsv'}")
