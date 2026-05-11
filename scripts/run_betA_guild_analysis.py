#!/usr/bin/env python3
"""Analyze the EQ betA guild census across all 312 metagenomes.

Steps:
1. Load census (per-bin K00108 counts) and total CDS per sample.
2. Build per-sample summary: total K00108 hits, K00108 density (per million CDS).
3. Cross-tabulate: bins encoding K00108 across samples, compartments.
4. Compare CSP1-2-only vs guild-level prediction of community Shannon.

Run after `betA_guild_census.tsv` and `total_cds_per_sample.tsv` are
fetched from unimatrix.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"

cen = pd.read_csv(CACHE / "betA_guild_census.tsv", sep="\t")
cds = pd.read_csv(CACHE / "total_cds_per_sample.tsv", sep="\t")
print(f"Census: {len(cen)} bin records across {cen['sample_id'].nunique()} samples")
print(f"CDS counter: {len(cds)} samples, "
      f"{cds['total_cds'].sum():,} total CDS")

# Re-derive compartment (bash regex missed F*/T* prefixes)
import re
def compartment_from_sid(sid):
    s = sid.split("_")[-1]
    m = re.match(r"^[A-Z]*[0-9]+([A-Z]+)r[0-9]+$", s)
    if not m: return "?"
    c = m.group(1)
    return {"PR":"rhizosphere","S":"surface","D":"deep"}.get(c, c)
cen["compartment"] = cen["sample_id"].apply(compartment_from_sid)
cds["compartment"] = cds["sample_id"].apply(compartment_from_sid)

# ---- Per-sample aggregation ----
samp = cen.groupby("sample_id").agg(
    n_bins=("bin_id","count"),
    n_K00108_strict=("n_K00108_strict","sum"),
    n_choline_dehyd_loose=("n_choline_dehyd_loose","sum"),
    n_bins_with_K00108=("n_K00108_strict", lambda x: (x>0).sum()),
).reset_index()
samp = samp.merge(cds[["sample_id","total_cds"]], on="sample_id", how="left")
samp["K00108_density_per_M_CDS"] = samp["n_K00108_strict"] / (samp["total_cds"] / 1e6)
# Re-derive compartment from sample_id (bash script's regex missed F*/T* prefixes)
samp["compartment"] = samp["sample_id"].apply(compartment_from_sid)
samp.to_csv(CACHE / "betA_per_sample_summary.tsv", sep="\t", index=False)

print("\n=== Per-compartment guild summary ===")
g = samp.groupby("compartment").agg(
    n_samples=("sample_id","count"),
    samples_with_K00108=("n_K00108_strict", lambda x: (x>0).sum()),
    median_n_K00108=("n_K00108_strict","median"),
    median_density=("K00108_density_per_M_CDS","median"),
    median_bins_with_K00108=("n_bins_with_K00108","median"),
)
print(g.to_string())

# ---- Per-bin guild members ----
print("\n=== K00108-encoding bins ===")
producers = cen[cen["n_K00108_strict"]>0].copy()
print(f"Total bin-records with K00108: {len(producers)}")
print(f"Unique bins (sample×bin): {producers[['sample_id','bin_id']].drop_duplicates().shape[0]}")
print(f"Distinct bin IDs across samples: {producers['bin_id'].nunique()}")
print(f"Spread across compartments:")
print(producers["compartment"].value_counts().to_string())

# Save per-bin producers
producers.to_csv(CACHE / "betA_producers_per_bin.tsv", sep="\t", index=False)

# ---- UniRef cluster extraction via remote bakta TSVs ----
# Save list of (sample, bin) to fetch UniRef IDs for. Done in companion script.
producers[["trip","site","compartment","sample_id","bin_id"]].to_csv(
    CACHE / "betA_producers_locator.tsv", sep="\t", index=False)

print("\nWrote:")
for p in ["betA_per_sample_summary.tsv", "betA_producers_per_bin.tsv",
          "betA_producers_locator.tsv"]:
    print(f"  cache/{p}")
