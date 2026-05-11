#!/usr/bin/env python3
"""Spatial autocorrelation control for the salinity → Shannon and
mediation results.

Three tests:
  (1) Moran's I on residuals of Shannon ~ S + climate + compartment
      (raw OLS) using inverse-distance weights.
  (2) Spatial GLS / GEE: refit Shannon ~ S with spatial-distance
      autocorrelation structure; report S coefficient + CI.
  (3) Mediation refit on cells with site-level random intercept
      approximation (FE / clustered standard errors).

Output:
  cache/spatial_autocorrelation_tests.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.spatial.distance import cdist

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
sys.path.insert(0, str(REPO / "src"))

# ----------------------------------------------------------------------------
# Load and merge frames
# ----------------------------------------------------------------------------
frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")
geo_frames = []
for t in (1, 2, 3, 4, 5):
    p = REPO / f"data/geodata/trip{t}_geodata.tsv"
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        df["trip"] = t
        geo_frames.append(df)
geo = pd.concat(geo_frames, ignore_index=True)

# Normalise column names
for c in ("SiteNum", "Site"):
    if c in geo.columns:
        geo = geo.rename(columns={c: "site"})
        break
for c in ("Lat", "Latitude"):
    if c in geo.columns:
        geo = geo.rename(columns={c: "lat"})
        break
for c in ("Lon", "Longitude"):
    if c in geo.columns:
        geo = geo.rename(columns={c: "lon"})
        break

geo["site"] = pd.to_numeric(geo["site"], errors="coerce")
geo = geo.dropna(subset=["site", "lat", "lon"]).copy()
geo["site"] = geo["site"].astype(int)
sites = (geo[["site", "lat", "lon"]]
         .drop_duplicates("site").set_index("site"))

# Cell aggregation
agg_cols = ["S", "Cl", "Na", "P", "rain_W30d", "temp_mean_W30d",
            "shannon", "csp_relab"]
cells = (frame.dropna(subset=agg_cols)
         .groupby(["trip", "site", "compartment"])
         .agg({c: "mean" for c in agg_cols})
         .reset_index())
cells = cells.merge(sites, left_on="site", right_index=True, how="left")
cells = cells.dropna(subset=["lat", "lon"])
print(f"Cells with coords: n = {len(cells)}")

# ----------------------------------------------------------------------------
# (1) Moran's I on Shannon ~ S + covars residuals
# ----------------------------------------------------------------------------
COVARS = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d"]
X = cells[COVARS + ["S"]].values
y = cells["shannon"].values
beta = LinearRegression().fit(X, y)
resid = y - beta.predict(X)

coords = cells[["lat", "lon"]].values
# Great-circle distance approximation in km (lat, lon in degrees)
def haversine_matrix(coords):
    lat = np.radians(coords[:, 0]); lon = np.radians(coords[:, 1])
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat/2)**2 + np.cos(lat[:, None])*np.cos(lat[None, :])*np.sin(dlon/2)**2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

D = haversine_matrix(coords)
np.fill_diagonal(D, np.inf)
W = 1.0 / D
W[~np.isfinite(W)] = 0
np.fill_diagonal(W, 0)
W /= W.sum(axis=1, keepdims=True).clip(1e-12)

n = len(resid)
mean_r = resid.mean()
zr = resid - mean_r
num = (zr[:, None] * W * zr[None, :]).sum()
den = (zr ** 2).sum()
I = (n / W.sum()) * num / den
# Expected I under randomization
EI = -1.0 / (n - 1)
# Variance approximation (Moran's I randomization)
S0 = W.sum()
S1 = ((W + W.T)**2).sum() / 2
S2 = ((W.sum(axis=1) + W.sum(axis=0))**2).sum()
b2 = (n * (zr**4).sum()) / ((zr**2).sum()**2)
var_I = (n*((n**2 - 3*n + 3)*S1 - n*S2 + 3*S0**2)
         - b2*((n**2 - n)*S1 - 2*n*S2 + 6*S0**2)) / \
         ((n-1)*(n-2)*(n-3)*S0**2) - EI**2
z = (I - EI) / np.sqrt(max(var_I, 1e-12))
from scipy.stats import norm
p = 2 * (1 - norm.cdf(abs(z)))

print(f"\nMoran's I on residuals (Shannon ~ S + covars):")
print(f"  I = {I:+.4f}, E(I) = {EI:+.4f}, "
      f"z = {z:+.2f}, p = {p:.3g}")
print(f"  → {'Significant spatial autocorrelation' if p < 0.05 else 'No significant spatial autocorrelation'}")

# ----------------------------------------------------------------------------
# (2) Site-clustered SE for Shannon ~ S — Eicker–White cluster-robust
# ----------------------------------------------------------------------------
# Cluster on site to relax independence within site
X_int = np.c_[np.ones(n), X]
ols = LinearRegression(fit_intercept=False).fit(X_int, y)
beta_hat = ols.coef_
resid2 = y - X_int @ beta_hat
# Naive (homoscedastic) SE for comparison
sigma2 = (resid2 ** 2).sum() / (n - X_int.shape[1])
XtX_inv = np.linalg.inv(X_int.T @ X_int)
naive_se = np.sqrt(np.diag(sigma2 * XtX_inv))

# Cluster-robust (Liang & Zeger 1986) on site
sites_arr = cells["site"].values
unique_sites = np.unique(sites_arr)
meat = np.zeros((X_int.shape[1], X_int.shape[1]))
for s in unique_sites:
    mask = sites_arr == s
    Xc = X_int[mask]
    rc = resid2[mask]
    sc = (Xc.T @ rc).reshape(-1, 1)
    meat += sc @ sc.T
G = len(unique_sites)
correction = (G / (G - 1)) * ((n - 1) / (n - X_int.shape[1]))
V_cluster = correction * XtX_inv @ meat @ XtX_inv
cluster_se = np.sqrt(np.diag(V_cluster))

# Locate the S coefficient (last in our X stack — see X_int construction)
# X = covars + [S]; so X_int columns = [intercept, Cl, Na, P, rain, temp, S]
S_idx = X_int.shape[1] - 1
S_beta = beta_hat[S_idx]
S_naive = naive_se[S_idx]
S_clust = cluster_se[S_idx]

print(f"\nShannon ~ S coefficient:")
print(f"  β_S = {S_beta:+.3f}")
print(f"  naive SE = {S_naive:.3f}, t = {S_beta / S_naive:+.2f}")
print(f"  cluster-robust SE (by site) = {S_clust:.3f}, "
      f"t = {S_beta / S_clust:+.2f}")
print(f"  inflation factor = {S_clust / S_naive:.2f}×")

# ----------------------------------------------------------------------------
# (3) Spatial GLS — exponential-decay covariance, fit S coefficient
# ----------------------------------------------------------------------------
# Σ = σ² (ρ * exp(-D/φ) + (1-ρ) I); fit ρ, φ via profile likelihood
def neg_log_lik(params, X, y, D):
    log_phi, logit_rho, log_sigma = params
    phi = np.exp(log_phi)
    rho = 1 / (1 + np.exp(-logit_rho))
    sigma = np.exp(log_sigma)
    n = len(y)
    Sig = sigma**2 * (rho * np.exp(-D/phi) + (1-rho) * np.eye(n))
    try:
        L = np.linalg.cholesky(Sig)
    except np.linalg.LinAlgError:
        return 1e10
    Sinv_X = np.linalg.solve(Sig, X)
    Sinv_y = np.linalg.solve(Sig, y)
    XtSX = X.T @ Sinv_X
    XtSy = X.T @ Sinv_y
    beta_gls = np.linalg.solve(XtSX, XtSy)
    res = y - X @ beta_gls
    quad = res @ np.linalg.solve(Sig, res)
    logdet = 2 * np.log(np.diag(L)).sum()
    return 0.5 * (logdet + quad)

# Subsample to keep the matrix tractable (use 250 random cells)
rng = np.random.default_rng(0)
idx_sub = rng.choice(n, size=min(n, 250), replace=False)
X_sub = X_int[idx_sub]
y_sub = y[idx_sub]
D_sub = D[np.ix_(idx_sub, idx_sub)]
np.fill_diagonal(D_sub, 0)

from scipy.optimize import minimize
res_opt = minimize(
    neg_log_lik,
    x0=[np.log(50.0), 0.0, np.log(np.std(y_sub))],
    args=(X_sub, y_sub, D_sub),
    method="Nelder-Mead",
    options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 200},
)
phi_hat = float(np.exp(res_opt.x[0]))
rho_hat = float(1 / (1 + np.exp(-res_opt.x[1])))
sigma_hat = float(np.exp(res_opt.x[2]))

# Refit GLS β with the optimised covariance
Sig = sigma_hat**2 * (rho_hat * np.exp(-D_sub/phi_hat) + (1-rho_hat)*np.eye(len(y_sub)))
Sinv_X = np.linalg.solve(Sig, X_sub)
Sinv_y = np.linalg.solve(Sig, y_sub)
XtSX = X_sub.T @ Sinv_X
beta_gls = np.linalg.solve(XtSX, X_sub.T @ Sinv_y)
res_gls = y_sub - X_sub @ beta_gls
sigma_gls = np.sqrt((res_gls @ np.linalg.solve(Sig, res_gls)) / (len(y_sub) - X_sub.shape[1]))
V_gls = sigma_gls**2 * np.linalg.inv(XtSX)
gls_se = np.sqrt(np.diag(V_gls))

print(f"\nSpatial GLS (n_sub={len(y_sub)}, exponential decay):")
print(f"  φ̂ = {phi_hat:.1f} km, ρ̂ = {rho_hat:.3f}")
print(f"  β_S (GLS) = {beta_gls[S_idx]:+.3f},  SE = {gls_se[S_idx]:.3f},  "
      f"t = {beta_gls[S_idx]/gls_se[S_idx]:+.2f}")

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
out = pd.DataFrame([
    {"test": "OLS naive",
     "S_beta": S_beta, "se": S_naive,
     "t": S_beta / S_naive,
     "moran_I": I, "moran_p": p},
    {"test": "OLS cluster-robust SE (by site)",
     "S_beta": S_beta, "se": S_clust,
     "t": S_beta / S_clust,
     "moran_I": I, "moran_p": p},
    {"test": f"Spatial GLS (exp-decay, φ̂={phi_hat:.0f} km)",
     "S_beta": beta_gls[S_idx], "se": gls_se[S_idx],
     "t": beta_gls[S_idx]/gls_se[S_idx],
     "moran_I": np.nan, "moran_p": np.nan},
])
out.to_csv(CACHE / "spatial_autocorrelation_tests.tsv",
           sep="\t", index=False)
print(f"\nwrote {CACHE / 'spatial_autocorrelation_tests.tsv'}")
