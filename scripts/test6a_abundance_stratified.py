#!/usr/bin/env python3
"""TEST 6A: Abundance-stratified persistence.

For each persistence class (n_trips_present in 1..5), report the
distribution of read counts at detection events. If 1-trip OTUs
are detected at much lower counts than 5-trip OTUs, the 67%
ephemeral signal is largely a sub-detection sampling artifact.

Inputs:
  cache/feature_table.parquet
  cache/test6_persistence/per_OTU_site_persistence.parquet  (from Test 6)

Output:
  cache/test6_disconfirmation/abundance_per_persistence_class.tsv
  cache/test6_disconfirmation/test6a_summary.txt
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
OUT = CACHE / "test6_disconfirmation"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    persist = pd.read_parquet(CACHE / "test6_persistence" /
                               "per_OTU_site_persistence.parquet")
    print(f"feature_table: {ft.shape}", flush=True)
    print(f"persistence rows: {len(persist)}", flush=True)
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)

    # Build per (ASV, site, comp, trip) read count
    # — use SUM across the 3 reps per (site, comp, trip)
    sm_set = smeta.set_index("sample")
    bucket_cols = {}
    for s, m in sm_set.iterrows():
        key = (int(m["site"]), m["compartment"], int(m["trip"]))
        bucket_cols.setdefault(key, []).append(s)

    # Restrict to (site, comp) cells observed at all 5 trips
    cells_5trip = {}
    sc = {}
    for (site, comp, trip), cols in bucket_cols.items():
        sc.setdefault((site, comp), set()).add(trip)
    cells_5trip = {(s, c) for (s, c), tset in sc.items() if tset == {1, 2, 3, 4, 5}}
    print(f"(site, comp) cells with all 5 trips: {len(cells_5trip)}", flush=True)

    # For each OTU-site-comp where in cells_5trip, build per-trip read counts
    rec = []
    persist_lkp = persist.set_index(["ASV", "site", "compartment"])["n_trips_present"]
    for (site, comp) in cells_5trip:
        comp_cells = [(site, comp, t) for t in (1, 2, 3, 4, 5)]
        per_trip_sums = {}
        for cell in comp_cells:
            cols = bucket_cols.get(cell, [])
            if not cols: continue
            per_trip_sums[cell[2]] = ft[cols].sum(axis=1)
        # ASVs present at this (site, comp): index into persist
        try:
            asvs_here = persist[(persist["site"] == site)
                                & (persist["compartment"] == comp)
                                & (persist["n_trips_observed"] == 5)]
        except Exception:
            continue
        for _, prow in asvs_here.iterrows():
            asv = prow["ASV"]
            n_trips = int(prow["n_trips_present"])
            for trip, sums in per_trip_sums.items():
                if asv in sums.index and sums[asv] > 0:
                    rec.append({"ASV": asv, "site": site, "compartment": comp,
                                 "trip": trip, "reads": int(sums[asv]),
                                 "n_trips_present": n_trips})
    df = pd.DataFrame(rec)
    print(f"detection events: {len(df)}", flush=True)
    df.to_parquet(OUT / "abundance_per_persistence_class.parquet")

    # Summary: per (compartment, n_trips_present), distribution of reads
    print("\n=== Reads at detection per persistence class (median, p25, p75) ===")
    agg = df.groupby(["compartment", "n_trips_present"])["reads"].agg(
        ["count", "median",
         lambda x: float(np.percentile(x, 25)),
         lambda x: float(np.percentile(x, 75)),
         "max"]).rename(columns={"<lambda_0>": "p25", "<lambda_1>": "p75"})
    print(agg.to_string())
    agg.to_csv(OUT / "abundance_per_persistence_class.tsv", sep="\t")

    # How does median abundance scale with persistence class?
    print("\n=== Pivot: median read count per (compartment, persistence class) ===")
    piv = df.groupby(["compartment", "n_trips_present"])["reads"].median().unstack()
    print(piv.to_string())

    with open(OUT / "test6a_summary.txt", "w") as fh:
        fh.write("Test 6A: Abundance-stratified persistence\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Detection events analyzed: {len(df):,}\n")
        fh.write(f"(site, comp) cells with all 5 trips: {len(cells_5trip)}\n\n")
        fh.write("Per (compartment, n_trips_present) reads-at-detection:\n")
        fh.write(agg.to_string())
        fh.write("\n\nMedian reads per (compartment, persistence class):\n")
        fh.write(piv.to_string())
        fh.write("\n\nDISCONFIRMATION KEY:\n")
        fh.write("  If median reads for n_trips=1 is < 10 AND median for n_trips=5\n"
                 "    is > 50, the '67% ephemeral' is largely sub-detection sampling.\n")
        fh.write("  If medians are comparable across persistence classes, the\n"
                 "    ephemerality signal is real (not just abundance/detection).\n")
    print(f"\nWrote {OUT}/test6a_summary.txt")


if __name__ == "__main__":
    main()
