#!/usr/bin/env python3
"""TEST 4: Functional redundancy via Allison-Martiny slope.

Pairwise functional BC ~ taxonomic BC. Per (compartment), per trip.
Slope close to 1 = no redundancy (every taxonomic substitution carries
a functional cost). Slope << 1 = redundancy buffers function against
taxonomic turnover.

Inputs:
  cache/feature_table.parquet
  data/functional/picrust2/path_abun_unstrat.tsv

Outputs:
  cache/test4_allison_martiny/per_compartment_trip_slope.tsv
  cache/test4_allison_martiny/all_pair_data.parquet
  cache/test4_allison_martiny/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import linregress, spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "test4_allison_martiny"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    path = pd.read_csv(DATA / "functional" / "picrust2" / "path_abun_unstrat.tsv",
                       sep="\t", index_col=0)
    print(f"feature_table: {ft.shape}", flush=True)
    print(f"PICRUSt2 paths: {path.shape}", flush=True)

    common = sorted(set(ft.columns) & set(path.columns))
    print(f"shared samples: {len(common)}", flush=True)
    smeta = parse_samples_to_df(common)
    smeta["site"] = smeta["site"].astype(int)

    rel_t = ft[common].div(ft[common].sum(axis=0).replace(0, 1), axis=1)
    rel_f = path[common].div(path[common].sum(axis=0).replace(0, 1), axis=1)

    rows = []
    pair_records = []
    for (comp, trip), grp in smeta.groupby(["compartment", "trip"]):
        if len(grp) < 10: continue
        cols = grp["sample"].tolist()
        Tt = rel_t[cols].T.values
        Tf = rel_f[cols].T.values
        Dt = squareform(pdist(Tt, metric="braycurtis"))
        Df = squareform(pdist(Tf, metric="braycurtis"))
        iu = np.triu_indices(len(cols), k=1)
        x = Dt[iu]; y = Df[iu]
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if valid.sum() < 30: continue
        x = x[valid]; y = y[valid]
        slope, intercept, r, p, _ = linregress(x, y)
        rho, prho = spearmanr(x, y)
        rows.append({
            "compartment": comp, "trip": int(trip),
            "n_samples": len(cols), "n_pairs": int(valid.sum()),
            "slope_func_vs_tax": float(slope),
            "intercept": float(intercept),
            "pearson_r": float(r), "p_lin": float(p),
            "spearman_rho": float(rho), "p_spearman": float(prho),
            "median_taxonomic_BC": float(np.median(x)),
            "median_functional_BC": float(np.median(y)),
            "ratio_functional_to_taxonomic_BC":
                float(np.median(y) / np.median(x)) if np.median(x) > 0 else np.nan,
        })
        # Sample 200 pairs for plotting
        n_save = min(200, len(x))
        idx = np.random.default_rng(20260509).choice(len(x), n_save, replace=False)
        for i in idx:
            pair_records.append({"compartment": comp, "trip": int(trip),
                                  "tax_BC": float(x[i]), "func_BC": float(y[i])})
        print(f"  {comp:>11s} trip={trip}: slope={slope:+.3f}, r={r:.3f}, "
              f"med_tax={np.median(x):.3f}, med_func={np.median(y):.3f}, "
              f"ratio_func/tax={np.median(y)/np.median(x):.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_compartment_trip_slope.tsv", sep="\t", index=False)
    pdf = pd.DataFrame(pair_records)
    pdf.to_parquet(OUT / "all_pair_data.parquet")

    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 4: Allison-Martiny functional redundancy slope\n")
        fh.write("=" * 70 + "\n\n")
        fh.write("Per (compartment, trip) regression: functional_BC = a + b * taxonomic_BC\n\n")
        fh.write(df.to_string(index=False))
        fh.write("\n\nMEDIAN slope per compartment:\n")
        agg = df.groupby("compartment")[["slope_func_vs_tax",
                                            "ratio_functional_to_taxonomic_BC"]].median()
        fh.write(agg.to_string())
        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  slope ~ 1.0  -> NO REDUNDANCY (every taxonomic change -> functional change)\n")
        fh.write("  slope ~ 0.5  -> moderate redundancy\n")
        fh.write("  slope ~ 0.2  -> high redundancy (Allison-Martiny null hypothesis)\n")
        fh.write("\n  ratio_func/tax > 0.8 -> functions track taxa closely (no redundancy)\n")
        fh.write("  ratio_func/tax < 0.4 -> functions much less variable than taxa (redundant)\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
