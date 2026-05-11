#!/usr/bin/env python3
"""Leak-asymmetry test: do bins with K00108 (betA producers) tend to LACK
glycine-betaine uptake genes (consistent with passive leak), while bins
WITHOUT K00108 (potential dependents) tend to CARRY uptake genes?

Reads:
  cache/betA_guild_census.tsv         (per-bin K00108 counts)
  cache/betaine_uptake_census.tsv     (per-bin uptake gene counts)

Outputs:
  cache/leak_asymmetry_test.tsv        (2x2 contingency)
  cache/leak_asymmetry_per_bin.tsv     (joined per-bin)
  cache/leak_asymmetry_test.txt        (text summary)
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency

CACHE = Path(__file__).resolve().parents[1] / "cache"

prod = pd.read_csv(CACHE/"betA_guild_census.tsv", sep="\t")
upt  = pd.read_csv(CACHE/"betaine_uptake_census.tsv", sep="\t")

print(f"Producer census rows: {len(prod)}")
print(f"Uptake census rows: {len(upt)}")

joined = prod.merge(
    upt[["sample_id","bin_id","n_uptake_KEGG","n_uptake_loose"]],
    on=["sample_id","bin_id"], how="inner")
print(f"Joined rows: {len(joined)}")

joined["has_betA"] = (joined["n_K00108_strict"]>0).astype(int)
joined["has_uptake_strict"] = (joined["n_uptake_KEGG"]>0).astype(int)
joined["has_uptake_loose"] = (joined["n_uptake_loose"]>0).astype(int)

joined.to_csv(CACHE/"leak_asymmetry_per_bin.tsv", sep="\t", index=False)

print("\n=== 2x2 contingency: betA × uptake (strict KEGG only) ===")
ct = pd.crosstab(joined["has_betA"], joined["has_uptake_strict"],
                  rownames=["betA"], colnames=["uptake"])
print(ct.to_string())
if ct.shape == (2,2):
    odds, p_fisher = fisher_exact(ct.values, alternative="less")
    print(f"\nFisher's exact (one-sided: betA→fewer uptake): "
          f"OR={odds:.3f}, p={p_fisher:.3g}")
    chi2, p_chi, _, _ = chi2_contingency(ct.values)
    print(f"Chi-square: chi2={chi2:.2f}, p={p_chi:.3g}")
    # Fraction with uptake among producers vs non-producers
    pct_uptake_prod = ct.loc[1,1] / ct.loc[1].sum() * 100
    pct_uptake_nonp = ct.loc[0,1] / ct.loc[0].sum() * 100
    print(f"Uptake-positive among betA producers:    {pct_uptake_prod:.1f}%")
    print(f"Uptake-positive among non-producers:     {pct_uptake_nonp:.1f}%")

print("\n=== 2x2 contingency: betA × uptake (loose product match) ===")
ct2 = pd.crosstab(joined["has_betA"], joined["has_uptake_loose"],
                   rownames=["betA"], colnames=["uptake_loose"])
print(ct2.to_string())
if ct2.shape == (2,2):
    odds2, p2 = fisher_exact(ct2.values, alternative="less")
    print(f"Fisher's exact (one-sided): OR={odds2:.3f}, p={p2:.3g}")

# By compartment
print("\n=== Stratified by compartment ===")
for comp in ["surface","deep","rhizosphere"]:
    sub = joined[joined["compartment"]==comp]
    if len(sub) < 20: continue
    ct = pd.crosstab(sub["has_betA"], sub["has_uptake_strict"],
                      rownames=["betA"], colnames=["uptake"])
    print(f"\n{comp}: n={len(sub)}")
    print(ct.to_string())
    if ct.shape == (2,2) and (ct.values > 0).all():
        odds, p = fisher_exact(ct.values, alternative="less")
        print(f"OR={odds:.3f}, p={p:.3g}")

# Save summary
with open(CACHE/"leak_asymmetry_test.txt","w") as fh:
    fh.write("Leak-asymmetry test: betA producers should LACK uptake (passive leak)\n")
    fh.write("="*70 + "\n\n")
    fh.write(f"Total bins with both producer and uptake annotations: {len(joined)}\n")
    fh.write(f"betA-positive bins: {joined['has_betA'].sum()}\n")
    fh.write(f"Uptake-positive bins (strict KEGG): {joined['has_uptake_strict'].sum()}\n")
    fh.write(f"Uptake-positive bins (loose): {joined['has_uptake_loose'].sum()}\n")
    fh.write("\nSee leak_asymmetry_per_bin.tsv for raw per-bin records.\n")
    fh.write("Compartment-stratified contingency tables in script stdout.\n")
print(f"\nWrote cache/leak_asymmetry_test.txt and cache/leak_asymmetry_per_bin.tsv")
