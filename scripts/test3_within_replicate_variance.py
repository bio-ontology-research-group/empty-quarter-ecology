#!/usr/bin/env python3
"""TEST 3: Within-replicate vs between-site variance.

For each (compartment, trip) cell:
  - Pairwise BC across all sample-pairs
  - Split pairs into:
      * within-replicate: same (site, compartment, trip) -- different rep cores
      * between-site: same (compartment, trip) -- different sites
  - Compute median BC per category, ratio
  - PERMANOVA-style R^2: how much variance does site explain within (comp, trip)?

If within-replicate variance is comparable to (or > 30% of) between-site variance,
the distance-decay framework needs revising — micro-spatial heterogeneity
dominates over kilometer-scale geography.

Inputs:
  cache/feature_table.parquet   75469 ASVs x 1227 samples

Outputs:
  cache/test3_variance/within_vs_between_BC.tsv
  cache/test3_variance/permanova_per_compartment_trip.tsv
  cache/test3_variance/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
OUT = CACHE / "test3_variance"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260509)
N_PERM = 999  # PERMANOVA permutations


def permanova_r2(D: np.ndarray, groups: np.ndarray) -> tuple:
    """Compute PERMANOVA pseudo-F and R^2 + permutation p-value."""
    n = len(groups)
    if n < 4:
        return (np.nan, np.nan, np.nan)
    ss_total = float(np.sum(D ** 2)) / (2 * n)
    unique_g = np.unique(groups)
    ss_within = 0.0
    for g in unique_g:
        idx = np.where(groups == g)[0]
        if len(idx) < 2:
            continue
        sub = D[np.ix_(idx, idx)]
        ss_within += float(np.sum(sub ** 2)) / (2 * len(idx))
    ss_between = ss_total - ss_within
    a = len(unique_g)
    if a < 2 or n - a < 1:
        return (np.nan, np.nan, np.nan)
    pseudo_F = (ss_between / (a - 1)) / (ss_within / (n - a)) if ss_within > 0 else np.nan
    r2 = ss_between / ss_total if ss_total > 0 else np.nan
    # permutation p
    cnt = 0
    for _ in range(N_PERM):
        gp = RNG.permutation(groups)
        sw = 0.0
        for g in unique_g:
            idx = np.where(gp == g)[0]
            if len(idx) < 2: continue
            sub = D[np.ix_(idx, idx)]
            sw += float(np.sum(sub ** 2)) / (2 * len(idx))
        sb = ss_total - sw
        f = (sb / (a - 1)) / (sw / (n - a)) if sw > 0 else np.nan
        if f >= pseudo_F:
            cnt += 1
    p = (cnt + 1) / (N_PERM + 1)
    return (float(pseudo_F), float(r2), float(p))


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    smeta = parse_samples_to_df(ft.columns)
    smeta["site"] = smeta["site"].astype(int)
    print(f"feature_table: {ft.shape}, parsed: {len(smeta)}", flush=True)

    rel = ft.div(ft.sum(axis=0).replace(0, 1), axis=1)

    rows = []
    perm_rows = []
    for (comp, trip), grp in smeta.groupby(["compartment", "trip"]):
        if len(grp) < 6: continue
        cols = grp["sample"].tolist()
        sub_rel = rel[cols].T.values
        D = squareform(pdist(sub_rel, metric="braycurtis"))
        site_arr = grp["site"].values
        # Find within-rep pairs (same site, different rep) and between-site pairs
        n = len(cols)
        within = []; between = []
        for i in range(n):
            for j in range(i + 1, n):
                if site_arr[i] == site_arr[j]:
                    within.append(D[i, j])
                else:
                    between.append(D[i, j])
        if len(within) < 5 or len(between) < 5:
            continue
        rows.append({"compartment": comp, "trip": int(trip),
                     "n_samples": n, "n_sites": int(grp["site"].nunique()),
                     "n_within_rep_pairs": len(within),
                     "n_between_site_pairs": len(between),
                     "median_BC_within": float(np.median(within)),
                     "median_BC_between": float(np.median(between)),
                     "mean_BC_within": float(np.mean(within)),
                     "mean_BC_between": float(np.mean(between)),
                     "ratio_within_over_between": float(np.median(within)
                                                          / np.median(between))})
        # PERMANOVA: site explains how much variance?
        F, r2, p = permanova_r2(D, site_arr)
        perm_rows.append({"compartment": comp, "trip": int(trip),
                          "n_samples": n, "n_sites": int(grp["site"].nunique()),
                          "pseudo_F": F, "R2_site": r2, "p_value": p})
        print(f"  {comp:>11s} trip={trip} n={n} sites={grp['site'].nunique()} "
              f"within_med={np.median(within):.3f} between_med={np.median(between):.3f} "
              f"ratio={np.median(within)/np.median(between):.3f} "
              f"R2_site={r2:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "within_vs_between_BC.tsv", sep="\t", index=False)
    pdf = pd.DataFrame(perm_rows)
    pdf.to_csv(OUT / "permanova_per_compartment_trip.tsv", sep="\t", index=False)

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 3: Within-replicate vs between-site variance\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"PERMANOVA permutations: {N_PERM}\n\n")
        fh.write("Per (compartment, trip) BC and site-R^2:\n")
        joined = df.merge(pdf[["compartment", "trip", "R2_site", "p_value"]],
                          on=["compartment", "trip"])
        fh.write(joined.to_string(index=False))
        fh.write("\n\nMEDIAN of within/between ratio (across trips), per compartment:\n")
        agg = df.groupby("compartment")["ratio_within_over_between"].median()
        fh.write(agg.to_string())
        fh.write("\n\nMEDIAN R^2(site) across trips, per compartment:\n")
        agg2 = pdf.groupby("compartment")["R2_site"].median()
        fh.write(agg2.to_string())
        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  ratio_within/between < 0.3 -> distance-decay holds (sites differ)\n")
        fh.write("  ratio 0.5-0.7 -> sub-meter heterogeneity is comparable to inter-site\n"
                 "                 (distance-decay framework needs revising)\n")
        fh.write("  ratio > 0.8 -> SURPRISING: micro-spatial dominates over geography\n")
        fh.write("  R^2(site) < 0.3 -> sites do not explain much community variance\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
