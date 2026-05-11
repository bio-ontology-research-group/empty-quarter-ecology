#!/usr/bin/env python3
"""Per-compartment mediation: does the salinity → CSP1-2 → Shannon
indirect share survive in any one compartment when the panel is full?

Mirrors run_mediation_robustness.py variant C (cell-aggregated) but
fits separately within each of {surface, deep, rhizosphere}.
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

ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")
meta = pd.read_parquet(CACHE / "metadata_with_rainfall.parquet").set_index("sample")
xrf = pd.read_csv(REPO / "data" / "geochemistry" / "xrf_lab_table_all_trips.tsv",
                  sep="\t")
xrf["compartment"] = xrf["compartment"].str.lower()


def _h(col):
    x = col[col > 0].astype(float)
    if x.empty:
        return float("nan")
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())
meta["shannon"] = ft.apply(_h, axis=0).rename("shannon").reindex(meta.index)

ft_rel = ft.div(ft.sum(axis=0), axis=1)
csp_asvs = tax[tax["Taxon"].str.contains("CSP1-2|Dadabacteria",
                                          case=False, regex=True, na=False)].index
meta["csp_relab"] = ft_rel.loc[csp_asvs].sum(axis=0).reindex(meta.index).fillna(0)

ELEMENTS = ["S", "Cl", "Na", "P", "Fe", "Mn", "V", "K", "Ca", "Si"]
xrf_cell = (xrf.dropna(subset=["trip", "site", "compartment"])
              .groupby(["trip", "site", "compartment"])[ELEMENTS]
              .mean().reset_index())

frame = meta.reset_index().merge(
    xrf_cell, on=["trip", "site", "compartment"], how="left")

req = ["S", "csp_relab", "shannon", "Cl", "Na", "P",
       "rain_W30d", "temp_mean_W30d"]


def fit_med_cell(d, nboot=2000, seed=0):
    """Cell-aggregated mediation, 95% bootstrap CIs."""
    cell_keys = ["trip", "site", "compartment"]
    agg = (d.groupby(cell_keys)
           .agg({"shannon": "mean", "csp_relab": "mean",
                 "S": "mean", "Cl": "mean", "Na": "mean", "P": "mean",
                 "rain_W30d": "mean", "temp_mean_W30d": "mean"})
           .reset_index())
    agg = agg.dropna()
    covars = ["Cl", "Na", "P", "rain_W30d", "temp_mean_W30d"]
    rng = np.random.default_rng(seed)
    n = len(agg)
    a_b, b_b, cp_b = [], [], []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        db = agg.iloc[idx]
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
    a_b, b_b, cp_b = map(np.array, (a_b, b_b, cp_b))
    ind = a_b * b_b
    tot = ind + cp_b
    return {
        "n": n,
        "indirect": float(ind.mean()),
        "indirect_ci_lo": float(np.quantile(ind, 0.025)),
        "indirect_ci_hi": float(np.quantile(ind, 0.975)),
        "direct": float(cp_b.mean()),
        "direct_ci_lo": float(np.quantile(cp_b, 0.025)),
        "direct_ci_hi": float(np.quantile(cp_b, 0.975)),
        "total": float(tot.mean()),
        "prop_mediated": float(ind.mean() / tot.mean()) if tot.mean() != 0 else float("nan"),
    }


rows = []
for comp in ("surface", "deep", "rhizosphere"):
    sub = frame[frame.compartment == comp].dropna(subset=req)
    if len(sub) < 50:
        print(f"  skip {comp} — n={len(sub)}")
        continue
    out = fit_med_cell(sub)
    print(f"{comp:11s}  n_cells={out['n']:3d}  "
          f"indirect={out['indirect']:+.4f} "
          f"[{out['indirect_ci_lo']:+.3f},{out['indirect_ci_hi']:+.3f}]  "
          f"direct={out['direct']:+.4f} "
          f"[{out['direct_ci_lo']:+.3f},{out['direct_ci_hi']:+.3f}]  "
          f"prop_med={out['prop_mediated']:+.3f}")
    rows.append({"compartment": comp, **out})

# Also: dependent-pool as an alternative mediator
dep_genera = ['Herpetosiphon', 'Paenibacillus', 'Flavisolibacter', 'Ammoniphilus',
              'Streptomyces', 'Rubrobacter', 'Ectobacillus', 'Neobacillus',
              'Ramlibacter', 'Noviherbaspirillum', 'Nocardioides']
if "genus" not in tax.columns and "Genus" in tax.columns:
    tax = tax.rename(columns={"Genus": "genus"})
asv_to_genus = tax["genus"].reindex(ft.index).fillna("Unclassified")
genus_relab = ft_rel.assign(genus=asv_to_genus.values).groupby("genus").sum()
dep_present = [g for g in dep_genera if g in genus_relab.index]
meta["dep_pool"] = genus_relab.loc[dep_present].sum(axis=0).reindex(meta.index).fillna(0)

frame2 = meta.reset_index().merge(
    xrf_cell, on=["trip", "site", "compartment"], how="left")

print("\n--- alt mediator: dependent-pool ---")
for comp in ("surface", "deep", "rhizosphere"):
    sub = frame2[frame2.compartment == comp].dropna(
        subset=["S", "dep_pool", "shannon", "Cl", "Na", "P",
                "rain_W30d", "temp_mean_W30d"])
    if len(sub) < 50:
        continue
    sub2 = sub.rename(columns={"dep_pool": "csp_relab"}).copy()
    out = fit_med_cell(sub2)
    print(f"{comp:11s}  n_cells={out['n']:3d}  "
          f"indirect={out['indirect']:+.4f} "
          f"[{out['indirect_ci_lo']:+.3f},{out['indirect_ci_hi']:+.3f}]  "
          f"direct={out['direct']:+.4f}  "
          f"prop_med={out['prop_mediated']:+.3f}")
    rows.append({"compartment": f"{comp}_dep_pool", **out})

pd.DataFrame(rows).to_csv(CACHE / "causal_tier1_mediation_per_compartment.tsv",
                          sep="\t", index=False)
print(f"\nwrote {CACHE / 'causal_tier1_mediation_per_compartment.tsv'}")
