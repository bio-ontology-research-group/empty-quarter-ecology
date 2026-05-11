#!/usr/bin/env python3
"""Stoichiometric public-good model.

Question: given observed CSP1-2 abundance and per-cell betaine
production rate, can the leak alone supply the dependent guild's
osmoprotection demand under salt stress?

We use literature parameters with explicit ranges (low / mid / high
to bracket plausible biology) and report supply-to-demand ratio
across the salinity range observed in the EQ.

Inputs:
  cache/causal_frame_tier1.parquet  — for cell-level CSP1-2 and
                                       dependent-guild abundances.

Outputs:
  cache/stoichiometric_supply_demand.tsv   — per cell ratio
  Plus a console summary.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

# ----------------------------------------------------------------------------
# Per-cell parameters with low/mid/high bounds
# ----------------------------------------------------------------------------
# Cell density (cells per gram dry soil) — Whitman et al. 1998, soil
CELLS_PER_G = 1e10                       # 10^10 cells / g dry soil

# Per-cell betaine biosynthesis rate (mol betaine s^-1 per cell)
#   Bacterial choline-betaine synthesis ~ 1e-19 to 1e-17 mol/cell/s
#   under salt stress (Köster et al. 2003 OpuA flux; Boch et al. 1996).
BETA_PROD_LOW  = 1e-19
BETA_PROD_MID  = 3e-19
BETA_PROD_HIGH = 1e-18

# Fraction leaked vs retained
LEAK_FRAC_LOW  = 0.05
LEAK_FRAC_MID  = 0.20
LEAK_FRAC_HIGH = 0.50

# Per-cell betaine *demand* under salt stress (mol/cell/s, steady state)
#   B. subtilis with OpuA importing ~ 1e-19 to 5e-19 mol/cell/s under
#   ≥0.4 M NaCl (von Blohn et al. 1997, Höper et al. 2005).
DEM_LOW  = 5e-20
DEM_MID  = 2e-19
DEM_HIGH = 5e-19

# Time over which CSP1-2 supplies the guild (steady-state assumption)
# Dependent guild relative abundance ~ 5--15% across cells.

# ----------------------------------------------------------------------------
# Load observed abundances
# ----------------------------------------------------------------------------
frame = pd.read_parquet(CACHE / "causal_frame_tier1.parquet")
ft = pd.read_parquet(CACHE / "feature_table.parquet")
tax = pd.read_parquet(CACHE / "taxonomy.parquet")
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

cells = (frame.dropna(subset=["S", "csp_relab", "dep_pool", "shannon"])
         .groupby(["trip", "site", "compartment"])
         .agg({"S": "mean", "csp_relab": "mean",
               "dep_pool": "mean", "shannon": "mean"})
         .reset_index())
print(f"Cells: n = {len(cells)}")

# ----------------------------------------------------------------------------
# Compute supply-to-demand ratio per cell, three parameter sets
# ----------------------------------------------------------------------------
def supply_demand(csp_frac, dep_frac, prod_per_cell, leak_frac, demand):
    csp_cells = CELLS_PER_G * csp_frac
    dep_cells = CELLS_PER_G * dep_frac
    supply = csp_cells * prod_per_cell * leak_frac        # mol/g/s
    demand_total = dep_cells * demand                     # mol/g/s
    ratio = supply / demand_total
    return supply, demand_total, ratio

scenarios = [
    ("low",  BETA_PROD_LOW,  LEAK_FRAC_LOW,  DEM_HIGH),
    ("mid",  BETA_PROD_MID,  LEAK_FRAC_MID,  DEM_MID),
    ("high", BETA_PROD_HIGH, LEAK_FRAC_HIGH, DEM_LOW),
]

print("\n--- Per-scenario median supply / demand ratio ---")
rows = []
for label, prod, leak, dem in scenarios:
    sup, dem_t, ratio = supply_demand(
        cells["csp_relab"].values,
        cells["dep_pool"].values,
        prod, leak, dem,
    )
    print(f"  {label:5s}: prod={prod:.0e}, leak={leak:.0%}, demand={dem:.0e}")
    print(f"           median supply = {np.nanmedian(sup):.2e} mol/g/s")
    print(f"           median demand = {np.nanmedian(dem_t):.2e} mol/g/s")
    print(f"           median ratio  = {np.nanmedian(ratio):.2f}  "
          f"(IQR {np.nanquantile(ratio, 0.25):.2f}--"
          f"{np.nanquantile(ratio, 0.75):.2f})")
    rows.append({
        "scenario": label,
        "prod_per_cell": prod, "leak_frac": leak, "demand_per_cell": dem,
        "median_ratio": float(np.nanmedian(ratio)),
        "ratio_q25": float(np.nanquantile(ratio, 0.25)),
        "ratio_q75": float(np.nanquantile(ratio, 0.75)),
        "n_cells": int(len(cells)),
        "frac_cells_ratio_ge_1": float(np.mean(ratio >= 1.0)),
    })

ratio_df = pd.DataFrame(rows)
ratio_df.to_csv(CACHE / "stoichiometric_supply_demand.tsv",
                sep="\t", index=False)
print(f"\nwrote {CACHE / 'stoichiometric_supply_demand.tsv'}")

print("\nInterpretation:")
print("  ratio ≥ 1 → CSP1-2 leak alone can supply dependent-guild demand")
print("  ratio < 1 → CSP1-2 supply is insufficient; another source needed")
for r in rows:
    label = r["scenario"]
    print(f"  {label:5s} scenario: in {100*r['frac_cells_ratio_ge_1']:.0f}% "
          f"of cells the leak alone is sufficient")

# ----------------------------------------------------------------------------
# Sensitivity: what's the breakeven leak fraction?
# ----------------------------------------------------------------------------
# At median CSP and dep_pool, median S, what leak_frac gives ratio = 1?
m_csp = float(cells["csp_relab"].median())
m_dep = float(cells["dep_pool"].median())
breakeven = []
for label, prod, leak, dem in scenarios:
    needed = (m_dep * dem) / (m_csp * prod)
    breakeven.append((label, needed))
print("\nBreakeven leak fraction at median CSP1-2 and dep-pool abundances:")
for label, x in breakeven:
    feas = "feasible (<50%)" if x < 0.5 else "implausible (>50%)"
    print(f"  {label:5s}: required leak = {x:.0%}   ({feas})")
