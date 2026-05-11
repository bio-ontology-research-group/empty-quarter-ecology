#!/usr/bin/env python3
"""Re-run causal tier 1 (DML + mediation) on the all-trip XRF panel.

Mirrors the analytic code in ``notebooks/09_causal_tier1.qmd``.
Outputs:
    cache/causal_tier1_dml_ate.tsv         DML ATE estimates
    cache/causal_tier1_mediation.tsv       bootstrap mediation result
    cache/causal_tier1_panel_fe.tsv        within-transform panel FE
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

CACHE = REPO / "cache"

# 1) Load tables
ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")
meta = pd.read_parquet(CACHE / "metadata_with_rainfall.parquet").set_index("sample")
xrf = pd.read_csv(REPO / "data" / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                  sep="\t")
xrf["compartment"] = xrf["compartment"].str.lower()

# 2) Per-sample Shannon
def _h(col):
    x = col[col > 0].astype(float)
    if x.empty:
        return float("nan")
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())
shan = ft.apply(_h, axis=0).rename("shannon")
meta["shannon"] = shan.reindex(meta.index)

# 3) CSP1-2 / radiation guild relab
ft_rel = ft.div(ft.sum(axis=0), axis=1)
csp_asvs = tax[tax["Taxon"].str.contains("CSP1-2|Dadabacteria",
                                          case=False, regex=True, na=False)].index
meta["csp_relab"] = ft_rel.loc[csp_asvs].sum(axis=0).reindex(meta.index).fillna(0)

if "genus" not in tax.columns and "Genus" in tax.columns:
    tax = tax.rename(columns={"Genus": "genus"})

rad = ["Rubrobacter", "Deinococcus", "Geodermatophilus",
       "Modestobacter", "Blastococcus"]
asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")
rad_asvs = asv_to_genus.index[asv_to_genus.isin(rad)]
meta["rad_relab"] = ft_rel.loc[rad_asvs].sum(axis=0).reindex(meta.index).fillna(0)

# 4) Climate windows
np_path = CACHE / "nasa_power_daily.parquet"
if np_path.exists():
    cx = pd.read_parquet(np_path)
    cx["Date"] = pd.to_datetime(cx["Date"])
    uv_col, rad_col, et_col, sm_col = ("ALLSKY_SFC_UV_INDEX",
                                        "ALLSKY_SFC_SW_DWN",
                                        "EVPTRNS", "GWETTOP")
    meta_r = meta.reset_index()
    if "trip_date" in meta_r.columns:
        meta_r["trip_date"] = pd.to_datetime(meta_r["trip_date"])
        rows = []
        for _, r in meta_r.iterrows():
            s, d = r.site, r.trip_date
            subset = cx[(cx.Site == s) & (cx.Date <= d)]
            w7 = subset[subset.Date > d - pd.Timedelta(days=7)]
            w30 = subset[subset.Date > d - pd.Timedelta(days=30)]
            rows.append({"uv_W7d": w7[uv_col].mean(),
                         "uv_W30d": w30[uv_col].mean(),
                         "rad_W30d": w30[rad_col].mean(),
                         "et0_W30d": w30[et_col].mean(),
                         "sm_W7d": w7[sm_col].mean(),
                         "sm_W30d": w30[sm_col].mean()})
        clim = pd.DataFrame(rows, index=meta_r.index)
        meta = pd.concat([meta_r, clim], axis=1).set_index("sample")
        print(f"  merged NASA-POWER windows: {clim.shape[1]} new vars")

# 5) XRF cell panel (trip x site x compartment) → broadcast to 16S samples
ELEMENTS = ["S", "Cl", "Na", "P", "Fe", "Mn", "V", "K", "Ca", "Si", "Zn"]
xrf_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
              .groupby(["trip", "site", "compartment"])[ELEMENTS]
              .mean().reset_index())
print(f"  XRF cell panel: {xrf_cell.shape}")

frame = meta.reset_index().merge(
    xrf_cell, on=["trip", "site", "compartment"], how="left")
frame["compartment_num"] = frame.compartment.map(
    {"surface": 0, "deep": 1, "rhizosphere": 2})
frame.to_parquet(CACHE / "causal_frame_tier1.parquet")
print(f"  analysis frame: {frame.shape}")
print(f"  with full XRF (S,P,Cl,Na,Fe,V): "
      f"{frame.dropna(subset=['S','P','Cl','Na','Fe','V']).shape[0]} samples")

# 6) DML ATE
def double_ml_ate(Y, T, X, seed=0, K=5):
    n = len(Y)
    kf = KFold(n_splits=K, shuffle=True, random_state=seed)
    Y_res = np.zeros(n); T_res = np.zeros(n)
    for tr, te in kf.split(X):
        my = RandomForestRegressor(n_estimators=200, min_samples_leaf=5,
                                    random_state=seed, n_jobs=-1).fit(X[tr], Y[tr])
        mt = RandomForestRegressor(n_estimators=200, min_samples_leaf=5,
                                    random_state=seed, n_jobs=-1).fit(X[tr], T[tr])
        Y_res[te] = Y[te] - my.predict(X[te])
        T_res[te] = T[te] - mt.predict(X[te])
    theta = (T_res @ Y_res) / (T_res @ T_res)
    resid = Y_res - theta * T_res
    se = np.sqrt(np.mean(resid**2) / np.mean(T_res**2) / n)
    return float(theta), float(se)

xrf_complete = frame.dropna(subset=["S", "P", "Fe", "V", "Cl", "Na"]).copy()
print(f"  All-trip samples with full XRF: n={len(xrf_complete)} "
      f"trips={sorted(xrf_complete['trip'].unique().tolist())}")

elements = ["S", "Cl", "Na", "P", "Fe", "Mn", "V", "K"]
outcomes = ["shannon", "csp_relab", "rad_relab"]

dml_rows = []
np.random.seed(0)
for y in outcomes:
    for e in elements:
        sub = xrf_complete.dropna(
            subset=[y, e] + [x for x in elements if x != e]
            + ["rain_W30d", "temp_mean_W30d"]).copy()
        if len(sub) < 40:
            continue
        X = sub[[x for x in elements if x != e]
                + ["rain_W30d", "temp_mean_W30d", "compartment_num"]].values
        T = sub[e].values
        Y = sub[y].values
        try:
            ate, se = double_ml_ate(Y, T, X, seed=0, K=5)
            ci_lo, ci_hi = ate - 1.96*se, ate + 1.96*se
            dml_rows.append({"outcome": y, "treatment": e, "ate": ate,
                             "se": se, "ci_lo": ci_lo, "ci_hi": ci_hi,
                             "n": len(sub)})
        except Exception as err:
            dml_rows.append({"outcome": y, "treatment": e,
                             "ate": np.nan, "n": len(sub),
                             "err": str(err)[:60]})

dml_df = pd.DataFrame(dml_rows)
dml_df["significant"] = (dml_df.ci_lo * dml_df.ci_hi) > 0
dml_df.to_csv(CACHE / "causal_tier1_dml_ate.tsv", sep="\t", index=False)
print("\nDML ATE estimates (Shannon outcome, all-trip panel):")
print(dml_df[dml_df.outcome == "shannon"][
      ["treatment", "ate", "ci_lo", "ci_hi", "n", "significant"]].to_string(index=False))

# 7) Mediation: S -> CSP1-2 -> Shannon
med = xrf_complete.dropna(subset=["S", "csp_relab", "shannon",
                                  "Cl", "Na", "P", "rain_W30d",
                                  "temp_mean_W30d"]).copy()

def boot_mediation(d, T, M, Y, covars, nboot=2000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(d)
    a_b, b_b, cp_b = [], [], []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        db = d.iloc[idx]
        X_a = db[covars].values
        a = LinearRegression().fit(np.c_[X_a, db[T].values], db[M].values).coef_[-1]
        X_bT = db[covars + [T]].values
        b = LinearRegression().fit(np.c_[X_bT, db[M].values], db[Y].values).coef_[-1]
        cp = LinearRegression().fit(
            np.c_[X_a, db[M].values, db[T].values], db[Y].values).coef_[-1]
        a_b.append(a); b_b.append(b); cp_b.append(cp)
    a_b = np.array(a_b); b_b = np.array(b_b); cp_b = np.array(cp_b)
    ind = a_b * b_b
    tot = ind + cp_b
    return {
        "indirect": ind.mean(),
        "indirect_ci": (np.quantile(ind, 0.025), np.quantile(ind, 0.975)),
        "direct": cp_b.mean(),
        "direct_ci": (np.quantile(cp_b, 0.025), np.quantile(cp_b, 0.975)),
        "total": tot.mean(),
        "total_ci": (np.quantile(tot, 0.025), np.quantile(tot, 0.975)),
        "prop_mediated": (ind.mean() / tot.mean()) if tot.mean() != 0 else np.nan,
    }

covars = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d", "compartment_num"]
out = boot_mediation(med, "S", "csp_relab", "shannon", covars, nboot=2000)
med_df = pd.DataFrame([{
    "treatment": "S", "mediator": "csp_relab", "outcome": "shannon",
    "n": len(med),
    "indirect_effect": out["indirect"],
    "indirect_ci_lo": out["indirect_ci"][0], "indirect_ci_hi": out["indirect_ci"][1],
    "direct_effect": out["direct"],
    "direct_ci_lo": out["direct_ci"][0], "direct_ci_hi": out["direct_ci"][1],
    "total_effect": out["total"],
    "total_ci_lo": out["total_ci"][0], "total_ci_hi": out["total_ci"][1],
    "prop_mediated": out["prop_mediated"]
}])
med_df.to_csv(CACHE / "causal_tier1_mediation.tsv", sep="\t", index=False)
print("\nMediation (S → CSP1-2 → Shannon, all-trip panel):")
print(med_df.iloc[0].to_string())
