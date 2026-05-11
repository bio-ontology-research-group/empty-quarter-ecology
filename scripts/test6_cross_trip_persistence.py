#!/usr/bin/env python3
"""TEST 6: Cross-trip OTU persistence distribution.

For each (OTU, site, compartment), count number of trips it's detected in.
Builds the "trip-presence" distribution per compartment, classifies OTUs as:
  - Persistent core: detected in all 5 trips at the site
  - Recurrent: 3-4 trips
  - Intermittent: 2 trips
  - Episodic: 1 trip only

Also reports the distribution at GLOBAL level (across-sites count).

Inputs:
  cache/feature_table.parquet
  cache/taxonomy.parquet

Outputs:
  cache/test6_persistence/per_OTU_site_persistence.parquet
  cache/test6_persistence/persistence_summary_per_compartment.tsv
  cache/test6_persistence/abundant_otu_persistence.tsv
  cache/test6_persistence/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
OUT = CACHE / "test6_persistence"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    print(f"feature_table: {ft.shape}, parsed: {len(smeta)}", flush=True)

    # Build presence matrix: ASV x (site, comp) -> trip-presence count
    # Sample is detected if relabund > 0 OR raw count > 1; use count > 0 here.
    # First aggregate: for each (site, compartment, trip), is the OTU present?
    # An OTU is "present" at (site, compartment) on trip T if it occurs in
    # ANY rep at that (site, compartment, trip).
    print("Building per-(site, compartment, trip) presence...", flush=True)

    # Reshape feature table to long format only for non-zero entries (memory)
    # Use sparse approach: per (site, comp, trip) bucket, get list of present ASVs.
    smeta_set = smeta.set_index("sample")
    bucket_cols = {}
    for sample, m in smeta_set.iterrows():
        key = (int(m["site"]), m["compartment"], int(m["trip"]))
        bucket_cols.setdefault(key, []).append(sample)

    # For each bucket, OR-presence across replicate samples
    # Build per-bucket presence vector (boolean ASV)
    asv_index = ft.index
    bucket_presence = {}
    for key, cols in bucket_cols.items():
        sub = ft[cols]
        bucket_presence[key] = (sub.sum(axis=1) > 0).values  # ndarray

    # Now per (site, compartment) count trip-presence per OTU
    print("Aggregating trip presence per OTU per (site, compartment)...", flush=True)
    by_sc = {}
    for (site, comp, trip), pres in bucket_presence.items():
        by_sc.setdefault((site, comp), []).append(pres)

    rows = []
    for (site, comp), trips_pres in by_sc.items():
        n_trips_observed = len(trips_pres)  # how many trips covered this (site, comp)
        if n_trips_observed == 0: continue
        stack = np.vstack(trips_pres)  # n_trips x n_ASV
        trip_count_per_otu = stack.sum(axis=0)  # ASV-level int
        for asv_idx, k in enumerate(trip_count_per_otu):
            if k == 0: continue
            rows.append({"site": site, "compartment": comp,
                          "ASV": asv_index[asv_idx],
                          "n_trips_present": int(k),
                          "n_trips_observed": int(n_trips_observed)})
    persist_df = pd.DataFrame(rows)
    print(f"per-OTU-site rows (where trips_present > 0): {len(persist_df)}", flush=True)
    persist_df.to_parquet(OUT / "per_OTU_site_persistence.parquet")

    # Per-compartment distribution: of OTUs that are observed at >=1 trip,
    # what fraction is observed at all 5 trips?
    summary_rows = []
    for comp in ["rhizosphere", "surface", "deep"]:
        sub = persist_df[persist_df["compartment"] == comp]
        if len(sub) == 0: continue
        # restrict to (site, comp) cells with all 5 trips covered for fair comparison
        sub5 = sub[sub["n_trips_observed"] == 5]
        n5 = len(sub5)
        for k in range(1, 6):
            cnt = int((sub5["n_trips_present"] == k).sum())
            summary_rows.append({"compartment": comp, "trips_present": k,
                                  "n_OTU_site_records": cnt,
                                  "frac_of_records": cnt / n5 if n5 else np.nan})
    sumdf = pd.DataFrame(summary_rows)
    sumdf.to_csv(OUT / "persistence_summary_per_compartment.tsv",
                 sep="\t", index=False)

    # Abundant OTU subset: top 1% by total abundance
    totals = ft.sum(axis=1)
    n_top = max(int(len(ft) * 0.01), 200)
    top_asvs = totals.sort_values(ascending=False).head(n_top).index
    print(f"\nTop-{n_top} ASVs by total abundance — their persistence distribution:",
          flush=True)
    top_persist = persist_df[persist_df["ASV"].isin(top_asvs)
                              & (persist_df["n_trips_observed"] == 5)]
    top_persist["genus"] = top_persist["ASV"].map(tax["genus"].fillna("Unclassified"))
    top_summary = (top_persist.groupby("compartment")["n_trips_present"]
                              .value_counts().unstack(fill_value=0))
    print(top_summary.to_string())
    top_persist.to_csv(OUT / "abundant_otu_persistence.tsv", sep="\t", index=False)

    # Compute fraction of TOTAL READS attributable to persistent vs episodic OTUs
    print("\nFraction of total reads in each persistence class (per compartment):",
          flush=True)
    rd_summary_rows = []
    for comp in ["rhizosphere", "surface", "deep"]:
        comp_samples = [s for s in ft.columns
                        if s in smeta_set.index
                        and smeta_set.loc[s, "compartment"] == comp]
        if len(comp_samples) == 0: continue
        sub_ft = ft[comp_samples]
        total_reads = sub_ft.values.sum()
        sub_persist = persist_df[(persist_df["compartment"] == comp)
                                  & (persist_df["n_trips_observed"] == 5)]
        for k in range(1, 6):
            asvs_k = set(sub_persist[sub_persist["n_trips_present"] == k]["ASV"])
            asv_idx_k = sub_ft.index.intersection(asvs_k)
            reads_k = sub_ft.loc[asv_idx_k].values.sum() if len(asv_idx_k) else 0
            rd_summary_rows.append({"compartment": comp, "trips_present": k,
                                     "n_OTU": len(asvs_k),
                                     "total_reads": int(reads_k),
                                     "frac_of_compartment_reads":
                                         reads_k / total_reads if total_reads else np.nan})
    rd_df = pd.DataFrame(rd_summary_rows)
    print(rd_df.to_string(index=False))
    rd_df.to_csv(OUT / "reads_by_persistence_class.tsv", sep="\t", index=False)

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 6: Cross-trip OTU persistence\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("Trip-presence distribution per compartment "
                 "(restricted to (site, comp) cells observed at all 5 trips):\n\n")
        fh.write(sumdf.to_string(index=False))
        fh.write("\n\nFraction of compartment reads carried by each "
                 "persistence class:\n")
        fh.write(rd_df.to_string(index=False))
        fh.write("\n\nTop-1% abundant OTUs persistence (counts of OTU-site records):\n")
        fh.write(top_summary.to_string())
        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  if frac_5trips dominates (>0.4): communities are stable\n"
                 "    -> dispersal-mixing operates on existing populations\n")
        fh.write("  if frac_1trip dominates (>0.4): communities are episodic blooms\n"
                 "    -> seed-bank / pulse-driven assembly\n")
        fh.write("  intermediate -> mixed (typical soil)\n")
        fh.write("\nReads-fraction comparison: how much of the SIGNAL is from\n"
                 "  persistent vs episodic taxa? Often persistent taxa carry the\n"
                 "  bulk of reads even if they're a minority of taxa.\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
