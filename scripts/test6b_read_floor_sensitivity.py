#!/usr/bin/env python3
"""TEST 6B: Read-floor sensitivity for persistence.

Re-do Test 6 with progressively higher detection floors:
  floor in [1, 5, 10, 50, 100, 500] reads per sample.

If at floor=10 the 1-trip ephemeral fraction drops from 67% to ~20%,
those 47 percentage points were sampling noise.

Inputs:
  cache/feature_table.parquet

Output:
  cache/test6_disconfirmation/floor_sensitivity.tsv
  cache/test6_disconfirmation/test6b_summary.txt
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

FLOORS = [1, 5, 10, 50, 100, 500]


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    print(f"feature_table: {ft.shape}", flush=True)

    sm_set = smeta.set_index("sample")
    bucket_cols = {}
    for s, m in sm_set.iterrows():
        key = (int(m["site"]), m["compartment"], int(m["trip"]))
        bucket_cols.setdefault(key, []).append(s)

    # restrict to (site, comp) cells with all 5 trips
    sc = {}
    for (site, comp, trip) in bucket_cols:
        sc.setdefault((site, comp), set()).add(trip)
    cells_5trip = {(s, c) for (s, c), tset in sc.items() if tset == {1, 2, 3, 4, 5}}
    print(f"(site, comp) cells with all 5 trips: {len(cells_5trip)}", flush=True)

    # Pre-compute per (site, comp, trip) sum across reps
    per_cell_sum = {}
    for cell, cols in bucket_cols.items():
        per_cell_sum[cell] = ft[cols].sum(axis=1).values

    asv_index = ft.index
    rows = []
    for floor in FLOORS:
        # For each (site, comp) in cells_5trip, count trips where each ASV
        # has >= floor reads (summed across reps).
        # Then aggregate the persistence-class distribution per compartment.
        per_comp = {"rhizosphere": np.zeros(6), "surface": np.zeros(6),
                     "deep": np.zeros(6)}  # index 0..5
        per_comp_total_rec = {"rhizosphere": 0, "surface": 0, "deep": 0}
        for (site, comp) in cells_5trip:
            stack = []
            for trip in (1, 2, 3, 4, 5):
                arr = per_cell_sum.get((site, comp, trip))
                if arr is None: continue
                stack.append(arr >= floor)
            if len(stack) != 5: continue
            S = np.vstack(stack)
            counts = S.sum(axis=0)  # ASV-level: in how many trips present
            # Distribution of ASVs by trip-count (excluding 0)
            for k in range(1, 6):
                per_comp[comp][k] += int((counts == k).sum())
            per_comp_total_rec[comp] += int((counts >= 1).sum())
        for comp in ("rhizosphere", "surface", "deep"):
            tot = per_comp_total_rec[comp]
            for k in range(1, 6):
                rows.append({"floor_reads": floor, "compartment": comp,
                              "trips_present": k,
                              "n_records": int(per_comp[comp][k]),
                              "frac_of_records": (per_comp[comp][k] / tot
                                                   if tot > 0 else 0)})
        print(f"  floor={floor:>4d}: rhizo n={per_comp_total_rec['rhizosphere']:,} "
              f"surface n={per_comp_total_rec['surface']:,} "
              f"deep n={per_comp_total_rec['deep']:,}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "floor_sensitivity.tsv", sep="\t", index=False)

    # Make a clear table: 1-trip-only fraction per (compartment, floor)
    one_trip = df[df["trips_present"] == 1]
    five_trip = df[df["trips_present"] == 5]
    pivot1 = one_trip.pivot_table(index="compartment", columns="floor_reads",
                                    values="frac_of_records").round(3)
    pivot5 = five_trip.pivot_table(index="compartment", columns="floor_reads",
                                    values="frac_of_records").round(3)

    print("\n=== fraction 1-trip-only at increasing read floor ===")
    print(pivot1.to_string())
    print("\n=== fraction 5-trip-persistent at increasing read floor ===")
    print(pivot5.to_string())

    with open(OUT / "test6b_summary.txt", "w") as fh:
        fh.write("Test 6B: Read-floor sensitivity for OTU persistence\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Floor levels tested: {FLOORS}\n")
        fh.write(f"(site, comp) cells with all 5 trips: {len(cells_5trip)}\n\n")
        fh.write("Fraction of OTU-records appearing in ONLY 1 of 5 trips,\n"
                 "by detection floor:\n")
        fh.write(pivot1.to_string())
        fh.write("\n\nFraction in ALL 5 of 5 trips:\n")
        fh.write(pivot5.to_string())
        fh.write("\n\nDISCONFIRMATION KEY:\n")
        fh.write("  - if 1-trip fraction at floor=10 drops to <0.4 (from 0.67),\n"
                 "    >40% of the original ephemerality signal is sub-detection\n"
                 "    sampling noise.\n")
        fh.write("  - if 5-trip fraction at floor=10 RISES to >0.1 (from 0.01-0.03),\n"
                 "    the apparent 'tiny core' was deflated by background.\n")
    print(f"\nWrote {OUT}/test6b_summary.txt")


if __name__ == "__main__":
    main()
