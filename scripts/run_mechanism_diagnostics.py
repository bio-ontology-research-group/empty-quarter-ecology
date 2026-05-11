#!/usr/bin/env python3
"""Diagnose where the salinity → diversity mechanism actually flows.

Four discriminating tests on the all-trip XRF panel, all
cell-aggregated (n=622 trip×site×compartment cells, with the same
covariates as the published mediation: Cl, Na, P, rain_W30d,
temp_mean_W30d):

  1. a/b path decomposition  (per compartment + pooled).
       a:    S → CSP1-2  | covars
       b:    CSP1-2 → Shannon | S, covars
       prod: a*b
     A weak a-path means salinity does not hit CSP1-2 hard;
     a weak b-path means CSP1-2 abundance does not load on Shannon.

  2. Alternative mediators (cell-pooled):
       - csp_relab            (the published mediator)
       - dep_pool             (sum of 11 dependent genera)
       - n_cycle (MND1 + Nitrospira)
       - oligo_pool (top-degree network hubs minus CSP1-2)

  3. Non-linear (Hill / saturation) a-path
       Fits S → CSP1-2 as Hill curve M = M0 - V * S^n / (K^n + S^n)
       and reports whether the relationship saturates / has a knee.
       Then re-estimates a*b with the implied non-linear a slope at
       the median S.

  4. Presence/absence keystone
       Binarise CSP1-2 at the median (cell-level) and refit
       mediation with the binary mediator. Tests Power (1996)
       keystone definition: effect via presence rather than mass.

All outputs land in cache/causal_mechanism_diagnostics.tsv plus a
human-readable summary at cache/causal_mechanism_diagnostics.txt.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

CACHE = REPO / "cache"

ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")
meta = pd.read_parquet(CACHE / "metadata_with_rainfall.parquet").set_index("sample")
xrf = pd.read_csv(REPO / "data" / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                  sep="\t")
xrf["compartment"] = xrf["compartment"].str.lower()


def shannon_per_sample(ft):
    def _h(col):
        x = col[col > 0].astype(float)
        if x.empty:
            return float("nan")
        p = x / x.sum()
        return float(-(p * np.log(p)).sum())
    return ft.apply(_h, axis=0)


# ----- Setup ----------------------------------------------------------
meta["shannon"] = shannon_per_sample(ft).reindex(meta.index)
ft_rel = ft.div(ft.sum(axis=0), axis=1)

if "genus" not in tax.columns and "Genus" in tax.columns:
    tax = tax.rename(columns={"Genus": "genus"})
asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")

# CSP1-2 by Taxon string match
csp_asvs = tax[tax["Taxon"].str.contains("CSP1-2|Dadabacteria",
                                          case=False, regex=True, na=False)].index
meta["csp_relab"] = ft_rel.loc[csp_asvs].sum(axis=0).reindex(meta.index).fillna(0)

# Dependent pool (11 facilitated genera identified earlier)
DEP_GENERA = ['Herpetosiphon', 'Paenibacillus', 'Flavisolibacter', 'Ammoniphilus',
              'Streptomyces', 'Rubrobacter', 'Ectobacillus', 'Neobacillus',
              'Ramlibacter', 'Noviherbaspirillum', 'Nocardioides']
genus_relab = (ft_rel.assign(genus=asv_to_genus.values).groupby("genus").sum())
dep_present = [g for g in DEP_GENERA if g in genus_relab.index]
meta["dep_pool"] = genus_relab.loc[dep_present].sum(axis=0).reindex(meta.index).fillna(0)

# N-cycle: MND1 + Nitrospira
n_cycle = [g for g in ("MND1", "Nitrospira") if g in genus_relab.index]
meta["n_cycle"] = (genus_relab.loc[n_cycle].sum(axis=0)
                   if n_cycle else pd.Series(0.0, index=ft.columns)
                   ).reindex(meta.index).fillna(0)

# Oligotroph hubs from network nodes (non-CSP1-2 high-degree taxa)
def load_top_hubs():
    rows = []
    for c in ("surface", "deep", "rhizosphere"):
        nodes = pd.read_csv(CACHE / f"network_nodes_{c}.tsv", sep="\t")
        rows.append(nodes.assign(compartment=c))
    df = pd.concat(rows, ignore_index=True)
    df = df[(df.node != "CSP1-2") & (df.node != "MND1")
            & (df.node != "Nitrospira")]
    # take top 10 by degree across union of compartments
    return (df.groupby("node")["degree"].max()
            .sort_values(ascending=False).head(10).index.tolist())
oligo_hubs = load_top_hubs()
oligo_present = [g for g in oligo_hubs if g in genus_relab.index]
meta["oligo_pool"] = (genus_relab.loc[oligo_present].sum(axis=0)
                      if oligo_present
                      else pd.Series(0.0, index=ft.columns)
                      ).reindex(meta.index).fillna(0)
print(f"oligo hubs: {oligo_present}")

# XRF cell panel + frame
ELEMENTS = ["S", "Cl", "Na", "P", "Fe", "Mn", "V", "K", "Ca", "Si"]
xrf_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
              .groupby(["trip", "site", "compartment"])[ELEMENTS]
              .mean().reset_index())

frame = meta.reset_index().merge(
    xrf_cell, on=["trip", "site", "compartment"], how="left")

REQ = ["S", "Cl", "Na", "P", "rain_W30d", "temp_mean_W30d", "shannon"]
COVARS = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d"]


def cell_agg(df, mediator):
    keys = ["trip", "site", "compartment"]
    return (df.groupby(keys)
            .agg({**{c: "mean" for c in REQ + [mediator]}})
            .reset_index().dropna())


# ----- Test 1: a / b path decomposition --------------------------------
def ab_paths(d, mediator, comp_label, n_boot=2000, seed=0):
    """Return point estimates and 95% CIs for a, b, indirect, direct."""
    rng = np.random.default_rng(seed)
    n = len(d)
    a_b, b_b, c_b, cp_b = [], [], [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        db = d.iloc[idx]
        # a: S → M | covars
        X = db[COVARS].values
        a = LinearRegression().fit(np.c_[X, db["S"].values], db[mediator].values).coef_[-1]
        # b: M → Y | S, covars
        Xb = db[COVARS + ["S"]].values
        b = LinearRegression().fit(np.c_[Xb, db[mediator].values],
                                    db["shannon"].values).coef_[-1]
        # c': S → Y | M, covars
        cp = LinearRegression().fit(np.c_[X, db[mediator].values, db["S"].values],
                                     db["shannon"].values).coef_[-1]
        # c: S → Y | covars (total)
        c = LinearRegression().fit(np.c_[X, db["S"].values],
                                    db["shannon"].values).coef_[-1]
        a_b.append(a); b_b.append(b); cp_b.append(cp); c_b.append(c)
    a_arr, b_arr, cp_arr, c_arr = map(np.array, (a_b, b_b, cp_b, c_b))
    ind = a_arr * b_arr

    def ci(arr): return float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))
    return {
        "compartment": comp_label, "mediator": mediator, "n": n,
        "a": float(a_arr.mean()), "a_lo": ci(a_arr)[0], "a_hi": ci(a_arr)[1],
        "b": float(b_arr.mean()), "b_lo": ci(b_arr)[0], "b_hi": ci(b_arr)[1],
        "indirect": float(ind.mean()),
        "indirect_lo": ci(ind)[0], "indirect_hi": ci(ind)[1],
        "direct": float(cp_arr.mean()),
        "direct_lo": ci(cp_arr)[0], "direct_hi": ci(cp_arr)[1],
        "total": float(c_arr.mean()),
        "prop_med": float(ind.mean() / c_arr.mean()) if c_arr.mean() != 0 else float("nan"),
    }


buf = []
def log(msg=""):
    print(msg)
    buf.append(msg)


log("\n========== TEST 1: a/b path decomposition (cell-level) ==========")
log("Mediator = csp_relab")
log(f"{'compartment':12s} {'n':>4s} {'a (S→M)':>16s} {'b (M→Y|S)':>16s} {'indirect':>16s} {'direct':>16s} {'prop_med':>8s}")
results_t1 = []
for comp in ("ALL", "surface", "deep", "rhizosphere"):
    sub = frame if comp == "ALL" else frame[frame.compartment == comp]
    d = cell_agg(sub, "csp_relab")
    if len(d) < 30:
        continue
    r = ab_paths(d, "csp_relab", comp)
    results_t1.append(r)
    log(f"{comp:12s} {r['n']:4d} "
        f"{r['a']:+.4f}[{r['a_lo']:+.3f},{r['a_hi']:+.3f}] "
        f"{r['b']:+.3f}[{r['b_lo']:+.2f},{r['b_hi']:+.2f}] "
        f"{r['indirect']:+.4f}[{r['indirect_lo']:+.3f},{r['indirect_hi']:+.3f}] "
        f"{r['direct']:+.3f}[{r['direct_lo']:+.2f},{r['direct_hi']:+.2f}] "
        f"{r['prop_med']:+.3f}")

# ----- Test 2: Alternative mediators (pooled) --------------------------
log("\n========== TEST 2: Alternative mediators (pooled) ==========")
log(f"{'mediator':14s} {'n':>4s} {'a (S→M)':>16s} {'b (M→Y|S)':>16s} {'indirect':>16s} {'prop_med':>8s}")
results_t2 = []
for med in ("csp_relab", "dep_pool", "n_cycle", "oligo_pool"):
    d = cell_agg(frame, med)
    if len(d) < 30:
        continue
    r = ab_paths(d, med, "ALL")
    results_t2.append(r)
    log(f"{med:14s} {r['n']:4d} "
        f"{r['a']:+.4f}[{r['a_lo']:+.3f},{r['a_hi']:+.3f}] "
        f"{r['b']:+.3f}[{r['b_lo']:+.2f},{r['b_hi']:+.2f}] "
        f"{r['indirect']:+.4f}[{r['indirect_lo']:+.3f},{r['indirect_hi']:+.3f}] "
        f"{r['prop_med']:+.3f}")


# ----- Test 3: Non-linear (Hill) a-path --------------------------------
def hill(x, M0, V, K, n):
    return M0 - V * (np.power(x, n) / (np.power(K, n) + np.power(x, n)))


log("\n========== TEST 3: Non-linear (Hill saturation) a-path ==========")
log("Fit S → CSP1-2 as Hill saturation curve; report point fit & local slope at median S")
results_t3 = []
for comp in ("ALL", "surface", "deep", "rhizosphere"):
    sub = frame if comp == "ALL" else frame[frame.compartment == comp]
    d = cell_agg(sub, "csp_relab").dropna(subset=["S", "csp_relab"])
    if len(d) < 30 or d["S"].max() <= 0:
        continue
    x = d["S"].values
    y = d["csp_relab"].values
    M0_init = float(np.percentile(y, 80))
    V_init = float(M0_init - np.percentile(y, 20))
    K_init = float(np.median(x[x > 0])) if (x > 0).any() else 0.5
    try:
        popt, _ = curve_fit(hill, x, y, p0=[M0_init, V_init, K_init, 1.0],
                            bounds=([0, 0, 1e-3, 0.3], [1, 1, 50, 5]),
                            maxfev=20_000)
        M0, V, K, n_hill = popt
        # local slope of M w.r.t. S at median S
        S_med = float(np.median(x))
        # d/dS: -V * n*K^n*S^(n-1) / (K^n + S^n)^2
        denom = (K**n_hill + S_med**n_hill)**2
        slope_med = -V * n_hill * (K**n_hill) * (S_med ** (n_hill - 1)) / denom
        # linear a from same data
        Xb = np.c_[d[COVARS].values, x]
        a_lin = LinearRegression().fit(Xb, y).coef_[-1]
        log(f"  {comp:12s} n={len(d):3d}  M0={M0:.4f}  V={V:.4f}  K={K:.2f}  n={n_hill:.2f}  "
            f"d/dS @ med S({S_med:.2f}) = {slope_med:+.4f}  vs linear a = {a_lin:+.4f}")
        results_t3.append({
            "compartment": comp, "n": len(d),
            "M0": M0, "V": V, "K": K, "n_hill": n_hill,
            "slope_at_median_S": slope_med, "linear_a": a_lin,
            "S_median": S_med,
        })
    except Exception as e:
        log(f"  {comp:12s} hill fit failed: {e}")


# ----- Test 4: Presence/absence keystone -------------------------------
log("\n========== TEST 4: Presence/absence keystone ==========")
log("Binarise CSP1-2 at cell-level median; re-fit mediation with binary mediator.")
log(f"{'compartment':12s} {'n':>4s} {'a (S→Mbin)':>16s} {'b (Mbin→Y|S)':>16s} {'indirect':>16s} {'prop_med':>8s}")
results_t4 = []
for comp in ("ALL", "surface", "deep", "rhizosphere"):
    sub = frame if comp == "ALL" else frame[frame.compartment == comp]
    d = cell_agg(sub, "csp_relab")
    if len(d) < 30:
        continue
    cutoff = float(d["csp_relab"].median())
    d = d.copy()
    d["csp_bin"] = (d["csp_relab"] > cutoff).astype(float)
    r = ab_paths(d, "csp_bin", comp)
    results_t4.append(r)
    log(f"{comp:12s} {r['n']:4d} "
        f"{r['a']:+.4f}[{r['a_lo']:+.3f},{r['a_hi']:+.3f}] "
        f"{r['b']:+.3f}[{r['b_lo']:+.2f},{r['b_hi']:+.2f}] "
        f"{r['indirect']:+.4f}[{r['indirect_lo']:+.3f},{r['indirect_hi']:+.3f}] "
        f"{r['prop_med']:+.3f}")


# ----- Save -----------------------------------------------------------
all_rows = (
    [{"test": "T1_ab_path", **r} for r in results_t1]
    + [{"test": "T2_alt_mediator", **r} for r in results_t2]
    + [{"test": "T3_hill_apath", **r} for r in results_t3]
    + [{"test": "T4_presence", **r} for r in results_t4]
)
out_df = pd.DataFrame(all_rows)
out_df.to_csv(CACHE / "causal_mechanism_diagnostics.tsv", sep="\t", index=False)
with open(CACHE / "causal_mechanism_diagnostics.txt", "w") as fh:
    fh.write("\n".join(buf))
log(f"\nwrote {CACHE / 'causal_mechanism_diagnostics.tsv'}")
log(f"wrote {CACHE / 'causal_mechanism_diagnostics.txt'}")
