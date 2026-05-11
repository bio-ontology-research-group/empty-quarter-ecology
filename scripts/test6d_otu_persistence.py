#!/usr/bin/env python3
"""TEST 6D: Re-do persistence at 99%-OTU level (collapsed from ASVs).

If 1-trip ephemerality fraction drops substantially after clustering,
sequencing/PCR error inflation explained part of the original signal.

Inputs:
  cache/feature_table.parquet
  cache/test6_disconfirmation/asv_to_otu_99.tsv  (from unimatrix vsearch)

Output:
  cache/test6_disconfirmation/test6d_persistence_at_otu_level.tsv
  cache/test6_disconfirmation/test6d_summary.txt
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


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)

    # Load ASV -> OTU map
    map_df = pd.read_csv(OUT / "asv_to_otu_99.tsv", sep="\t",
                          header=None, names=["asv", "otu"])
    map_dict = dict(zip(map_df["asv"], map_df["otu"]))
    print(f"ASV->OTU map: {len(map_dict)} entries", flush=True)
    print(f"unique OTUs: {len(set(map_dict.values()))}", flush=True)

    # Aggregate feature_table to OTU level
    asv_to_otu = pd.Series(map_dict)
    # Restrict to ASVs present in both
    common_asvs = ft.index.intersection(asv_to_otu.index)
    print(f"common ASVs: {len(common_asvs)}", flush=True)
    ft_sub = ft.loc[common_asvs]
    otu_assignment = asv_to_otu.loc[common_asvs]
    # GroupBy OTU sum
    print("Aggregating ASV counts to OTU level...", flush=True)
    otu_table = ft_sub.groupby(otu_assignment.values).sum()
    print(f"OTU table: {otu_table.shape}", flush=True)

    # Now redo persistence at OTU level
    sm_set = smeta.set_index("sample")
    bucket_cols = {}
    for s, m in sm_set.iterrows():
        key = (int(m["site"]), m["compartment"], int(m["trip"]))
        bucket_cols.setdefault(key, []).append(s)
    sc = {}
    for (site, comp, trip) in bucket_cols:
        sc.setdefault((site, comp), set()).add(trip)
    cells_5trip = {(s, c) for (s, c), tset in sc.items() if tset == {1, 2, 3, 4, 5}}
    print(f"(site, comp) cells with all 5 trips: {len(cells_5trip)}", flush=True)

    # Per (site, comp, trip), get OTU presence (>0 reads aggregated across reps)
    per_cell_pres = {}
    for cell, cols in bucket_cols.items():
        per_cell_pres[cell] = (otu_table[cols].sum(axis=1) > 0).values

    rows = []
    per_comp = {"rhizosphere": np.zeros(6),
                "surface": np.zeros(6),
                "deep": np.zeros(6)}
    per_comp_total = {"rhizosphere": 0, "surface": 0, "deep": 0}
    for (site, comp) in cells_5trip:
        stack = []
        for trip in (1, 2, 3, 4, 5):
            arr = per_cell_pres.get((site, comp, trip))
            if arr is None: continue
            stack.append(arr)
        if len(stack) != 5: continue
        S = np.vstack(stack)
        counts = S.sum(axis=0)
        for k in range(1, 6):
            per_comp[comp][k] += int((counts == k).sum())
        per_comp_total[comp] += int((counts >= 1).sum())

    for comp in ("rhizosphere", "surface", "deep"):
        tot = per_comp_total[comp]
        for k in range(1, 6):
            rows.append({"compartment": comp, "trips_present": k,
                          "n_OTU_records": int(per_comp[comp][k]),
                          "frac_of_records": (per_comp[comp][k] / tot
                                                if tot > 0 else 0)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "test6d_persistence_at_otu_level.tsv", sep="\t", index=False)

    print("\n=== Persistence distribution at 99%-OTU level ===")
    pivot = df.pivot_table(index="compartment", columns="trips_present",
                            values="frac_of_records").round(3)
    print(pivot.to_string())

    # Compare to ASV-level (for reference, from previous run)
    asv_level = pd.read_csv(CACHE / "test6_persistence" /
                              "persistence_summary_per_compartment.tsv",
                              sep="\t")
    asv_pivot = asv_level.pivot_table(index="compartment",
                                        columns="trips_present",
                                        values="frac_of_records").round(3)
    print("\n=== ASV-level reference (from Test 6) ===")
    print(asv_pivot.to_string())

    print("\n=== DELTA: OTU - ASV (positive = increased fraction at OTU level) ===")
    delta = pivot - asv_pivot
    print(delta.round(3).to_string())

    with open(OUT / "test6d_summary.txt", "w") as fh:
        fh.write("Test 6D: Persistence at 99%-OTU level\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"ASVs collapsed from 75,469 to 37,518 99%-OTUs\n")
        fh.write(f"(20-50% reduction varies by compartment)\n\n")
        fh.write("Persistence distribution at OTU level:\n")
        fh.write(pivot.to_string())
        fh.write("\n\nFor comparison, ASV-level (Test 6):\n")
        fh.write(asv_pivot.to_string())
        fh.write("\n\nDELTA (OTU - ASV):\n")
        fh.write(delta.round(3).to_string())
        fh.write("\n\nDISCONFIRMATION KEY:\n")
        fh.write("  If 1-trip fraction at OTU drops by >0.10 (= >10 percentage\n"
                 "    points), sequencing/PCR error inflation explains a notable\n"
                 "    fraction of the apparent ephemerality.\n")
        fh.write("  If unchanged (delta < 0.05), error inflation is not the\n"
                 "    explanation -- the ephemerality is biological signal.\n")
    print(f"\nWrote {OUT}/test6d_summary.txt")


if __name__ == "__main__":
    main()
