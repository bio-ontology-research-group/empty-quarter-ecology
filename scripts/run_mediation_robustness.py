#!/usr/bin/env python3
"""Robustness checks for the salinity → CSP1-2 → Shannon mediation.

Produces three side-by-side fits to clarify what changes when the
panel is enlarged from Trip-5 only to all trips:

    A. Trip-5 only (sample level) — reproduces the published 88% number
       using the combined-panel pipeline.
    B. All-trip sample level — same code path, broadcasts XRF cell to
       all 16S replicates of that (trip, site, compartment) cell.
       Pseudo-replicates inflate n; treat as upper bound on power.
    C. All-trip cell level — aggregates Shannon, CSP1-2 and XRF
       within (trip, site, compartment); n = 725 independent cells.
       This is the conservative, properly-de-replicated panel.
    D. Clustered bootstrap (sample level, resamples whole cells) —
       same point estimate as B but standard errors that respect
       within-cell correlation.

Reads the same caches as ``run_causal_tier1.py``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from sklearn.linear_model import LinearRegression

CACHE = REPO / "cache"

# Inputs
ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")
meta = pd.read_parquet(CACHE / "metadata_with_rainfall.parquet").set_index("sample")
xrf = pd.read_csv(REPO / "data" / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                  sep="\t")
xrf["compartment"] = xrf["compartment"].str.lower()

# Per-sample Shannon
def _h(col):
    x = col[col > 0].astype(float)
    if x.empty:
        return float("nan")
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())
meta["shannon"] = ft.apply(_h, axis=0).rename("shannon").reindex(meta.index)

# CSP1-2 relab
ft_rel = ft.div(ft.sum(axis=0), axis=1)
csp_asvs = tax[tax["Taxon"].str.contains("CSP1-2|Dadabacteria",
                                          case=False, regex=True, na=False)].index
meta["csp_relab"] = ft_rel.loc[csp_asvs].sum(axis=0).reindex(meta.index).fillna(0)

# XRF cell panel
ELEMENTS = ["S", "Cl", "Na", "P", "Fe", "Mn", "V", "K", "Ca", "Si"]
xrf_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
              .groupby(["trip", "site", "compartment"])[ELEMENTS]
              .mean().reset_index())

# Build sample-level frame
frame = meta.reset_index().merge(
    xrf_cell, on=["trip", "site", "compartment"], how="left")
frame["compartment_num"] = frame.compartment.map(
    {"surface": 0, "deep": 1, "rhizosphere": 2})

req = ["S", "csp_relab", "shannon", "Cl", "Na", "P",
       "rain_W30d", "temp_mean_W30d", "compartment_num"]


def fit_mediation(d, nboot=2000, seed=0, cluster_col=None):
    """Return point estimates and bootstrap CIs for indirect/direct/total.

    If ``cluster_col`` is given, resample whole clusters (rows sharing
    the same value) on each bootstrap iteration.
    """
    covars = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d", "compartment_num"]
    rng = np.random.default_rng(seed)
    n = len(d)

    a_b, b_b, cp_b = [], [], []
    if cluster_col is None:
        for _ in range(nboot):
            idx = rng.integers(0, n, n)
            db = d.iloc[idx]
            X_a = db[covars].values
            a = LinearRegression().fit(
                np.c_[X_a, db["S"].values], db["csp_relab"].values).coef_[-1]
            X_bT = db[covars + ["S"]].values
            b = LinearRegression().fit(
                np.c_[X_bT, db["csp_relab"].values],
                db["shannon"].values).coef_[-1]
            cp = LinearRegression().fit(
                np.c_[X_a, db["csp_relab"].values, db["S"].values],
                db["shannon"].values).coef_[-1]
            a_b.append(a); b_b.append(b); cp_b.append(cp)
    else:
        clusters = d[cluster_col].unique()
        m = len(clusters)
        idx_by_cluster = {c: d.index[d[cluster_col] == c].tolist() for c in clusters}
        for _ in range(nboot):
            chosen = rng.choice(clusters, size=m, replace=True)
            sample_idx = []
            for c in chosen:
                sample_idx.extend(idx_by_cluster[c])
            db = d.loc[sample_idx]
            X_a = db[covars].values
            a = LinearRegression().fit(
                np.c_[X_a, db["S"].values], db["csp_relab"].values).coef_[-1]
            X_bT = db[covars + ["S"]].values
            b = LinearRegression().fit(
                np.c_[X_bT, db["csp_relab"].values],
                db["shannon"].values).coef_[-1]
            cp = LinearRegression().fit(
                np.c_[X_a, db["csp_relab"].values, db["S"].values],
                db["shannon"].values).coef_[-1]
            a_b.append(a); b_b.append(b); cp_b.append(cp)
    a_b = np.array(a_b); b_b = np.array(b_b); cp_b = np.array(cp_b)
    ind = a_b * b_b
    tot = ind + cp_b
    return {
        "n": n,
        "indirect": ind.mean(),
        "indirect_ci": (np.quantile(ind, 0.025), np.quantile(ind, 0.975)),
        "direct": cp_b.mean(),
        "direct_ci": (np.quantile(cp_b, 0.025), np.quantile(cp_b, 0.975)),
        "total": tot.mean(),
        "total_ci": (np.quantile(tot, 0.025), np.quantile(tot, 0.975)),
        "prop_mediated": float(ind.mean() / tot.mean()) if tot.mean() != 0 else float("nan"),
    }


def fmt(r):
    return (f"  n={r['n']:5d}  indirect={r['indirect']:+.4f} "
            f"[{r['indirect_ci'][0]:+.3f},{r['indirect_ci'][1]:+.3f}]  "
            f"direct={r['direct']:+.4f} "
            f"[{r['direct_ci'][0]:+.3f},{r['direct_ci'][1]:+.3f}]  "
            f"total={r['total']:+.4f}  "
            f"prop_med={r['prop_mediated']:+.3f}")


# A. Trip-5 only (replicate the published 88%)
d_a = (frame[frame.trip == 5].dropna(subset=req)
       .reset_index(drop=True))
A = fit_mediation(d_a)

# B. All trips, sample-level (broadcast)
d_b = frame.dropna(subset=req).reset_index(drop=True)
B = fit_mediation(d_b)

# C. All trips, per-cell aggregation
cell_keys = ["trip", "site", "compartment"]
cell_means = (d_b.groupby(cell_keys)
              .agg({"shannon": "mean", "csp_relab": "mean",
                    "S": "mean", "Cl": "mean", "Na": "mean", "P": "mean",
                    "rain_W30d": "mean", "temp_mean_W30d": "mean",
                    "compartment_num": "first"})
              .reset_index())
C = fit_mediation(cell_means)

# D. Clustered bootstrap on the broadcast frame
d_b = d_b.copy()
d_b["cell"] = (d_b["trip"].astype(str) + "_" + d_b["site"].astype(str)
               + "_" + d_b["compartment"].astype(str))
D = fit_mediation(d_b, cluster_col="cell")


print("\n=== Salinity → CSP1-2 → Shannon mediation, robustness ===")
print(f"A) Trip-5 only (sample-level):     {fmt(A)}")
print(f"B) All trips sample-level (broadcast): {fmt(B)}")
print(f"C) All trips cell-aggregated (n=cells): {fmt(C)}")
print(f"D) All trips, cluster-bootstrap (cells): {fmt(D)}")

rows = [
    {"variant": "A_trip5_sample", **A,
     "indirect_ci_lo": A["indirect_ci"][0], "indirect_ci_hi": A["indirect_ci"][1]},
    {"variant": "B_alltrip_sample", **B,
     "indirect_ci_lo": B["indirect_ci"][0], "indirect_ci_hi": B["indirect_ci"][1]},
    {"variant": "C_alltrip_cell",  **C,
     "indirect_ci_lo": C["indirect_ci"][0], "indirect_ci_hi": C["indirect_ci"][1]},
    {"variant": "D_alltrip_clusterboot", **D,
     "indirect_ci_lo": D["indirect_ci"][0], "indirect_ci_hi": D["indirect_ci"][1]},
]
out_df = pd.DataFrame(rows).drop(columns=["indirect_ci", "direct_ci", "total_ci"])
out_df.to_csv(CACHE / "causal_tier1_mediation_robustness.tsv", sep="\t", index=False)
print(f"\nwrote {CACHE / 'causal_tier1_mediation_robustness.tsv'}")
