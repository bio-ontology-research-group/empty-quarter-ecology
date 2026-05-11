#!/usr/bin/env python3
"""Population-level analyses using the relic-likelihood indicator.

Stratifies the EQ ASV pool by relic_score_full_gb (Track C; AUC 0.785) into
deciles. For each (site, comp, trip), computes:

  1. Per-sample relic vs alive read fraction
  2. Alpha diversity (Shannon, ASV count) of relic-only and alive-only subsets
  3. Beta diversity (BC) within and between sites for relic vs alive subsets
  4. Compartment + trip patterns of relic fraction
  5. Cosmopolitan-vs-endemic enrichment in relic vs alive
  6. Geographic patterns: per-site relic load
  7. Taxonomic decomposition: which taxa dominate past vs present
  8. CSP1-2 / Rubellimicrobium classification (their relic_score)

Outputs:
  cache/relic_population/
    per_sample_relic_fraction.tsv
    alpha_div_by_pool.tsv
    beta_div_within_between.tsv
    relic_fraction_by_compartment_trip.tsv
    cosmopolitan_enrichment.tsv
    per_site_relic_load.tsv
    taxonomy_relic_vs_alive.tsv
    csp_rubelli_status.tsv
    summary.txt
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
OUT = CACHE / "relic_population"
OUT.mkdir(parents=True, exist_ok=True)

# Decile cuts
RELIC_HIGH = 0.7
ALIVE_HIGH = 0.3


def main():
    print("Loading inputs ...", flush=True)
    ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                       "relic_indicator_with_damage_per_asv.tsv", sep="\t")
    print(f"  indicator: {ind.shape}", flush=True)
    ind = ind[["asv_id", "relic_score_full_gb", "log_mean_abund",
                 "frac_deep", "frac_surface", "frac_rhizosphere",
                 "emp_cosmo_90", "weighted_median_ratio"]].copy()
    ind["pool"] = np.where(ind["relic_score_full_gb"] >= RELIC_HIGH, "relic",
                              np.where(ind["relic_score_full_gb"] <= ALIVE_HIGH,
                                        "alive", "ambiguous"))
    print(f"  pool counts: relic={int((ind['pool']=='relic').sum())}, "
          f"alive={int((ind['pool']=='alive').sum())}, "
          f"ambig={int((ind['pool']=='ambiguous').sum())}",
          flush=True)

    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    print(f"  feature table: {ft.shape}", flush=True)
    smeta = parse_samples_to_df(ft.columns).set_index("sample")
    smeta["site"] = smeta["site"].astype(int)

    # 1. Per-sample relic vs alive READ fraction
    print("\n[1] per-sample relic / alive READ fractions ...", flush=True)
    relic_asvs = set(ind.loc[ind["pool"] == "relic", "asv_id"])
    alive_asvs = set(ind.loc[ind["pool"] == "alive", "asv_id"])
    ambig_asvs = set(ind.loc[ind["pool"] == "ambiguous", "asv_id"])
    rel_idx = ft.index.isin(relic_asvs)
    ali_idx = ft.index.isin(alive_asvs)
    amb_idx = ft.index.isin(ambig_asvs)
    sample_total = ft.sum(axis=0)
    relic_reads = ft.loc[rel_idx].sum(axis=0)
    alive_reads = ft.loc[ali_idx].sum(axis=0)
    ambig_reads = ft.loc[amb_idx].sum(axis=0)
    rec = pd.DataFrame({"sample": ft.columns,
                          "total_reads": sample_total.values,
                          "relic_reads": relic_reads.values,
                          "alive_reads": alive_reads.values,
                          "ambig_reads": ambig_reads.values})
    rec["relic_frac"] = rec["relic_reads"] / rec["total_reads"].replace(0, 1)
    rec["alive_frac"] = rec["alive_reads"] / rec["total_reads"].replace(0, 1)
    rec = rec.merge(smeta.reset_index(), on="sample", how="left")
    rec.to_csv(OUT / "per_sample_relic_fraction.tsv", sep="\t", index=False)

    # 4. Compartment+trip stratification
    print("\n[4] relic fraction by compartment x trip ...", flush=True)
    by_ct = (rec.groupby(["compartment", "trip"])
             .agg(n_samples=("sample", "count"),
                  median_relic_frac=("relic_frac", "median"),
                  mean_relic_frac=("relic_frac", "mean"),
                  p25_relic=("relic_frac", lambda x: float(np.percentile(x, 25))),
                  p75_relic=("relic_frac", lambda x: float(np.percentile(x, 75))),
                  median_alive_frac=("alive_frac", "median"))
             .reset_index())
    by_ct.to_csv(OUT / "relic_fraction_by_compartment_trip.tsv",
                  sep="\t", index=False)
    print("  per-compartment (across all trips):")
    print(rec.groupby("compartment")[["relic_frac", "alive_frac"]]
          .median().round(3).to_string())

    # 6. Per-site relic load
    print("\n[6] per-site relic load ...", flush=True)
    per_site = (rec.groupby("site")
                .agg(n_samples=("sample", "count"),
                     median_relic_frac=("relic_frac", "median"),
                     median_alive_frac=("alive_frac", "median"))
                .reset_index())
    per_site.to_csv(OUT / "per_site_relic_load.tsv", sep="\t", index=False)
    print(f"  median relic_frac across sites: "
          f"min={per_site['median_relic_frac'].min():.3f} "
          f"max={per_site['median_relic_frac'].max():.3f} "
          f"median={per_site['median_relic_frac'].median():.3f}", flush=True)

    # 2. Alpha diversity per pool (Shannon, ASV count)
    print("\n[2] alpha diversity per pool ...", flush=True)
    def shannon(x):
        x = x[x > 0]
        if len(x) == 0: return 0.0
        p = x / x.sum()
        return float(-(p * np.log(p)).sum())

    alpha = []
    for sample in ft.columns:
        col = ft[sample]
        for pool, mask in (("relic", rel_idx), ("alive", ali_idx),
                             ("all", np.ones(len(ft), dtype=bool))):
            sub = col.values[mask]
            alpha.append({"sample": sample, "pool": pool,
                            "asv_count": int((sub > 0).sum()),
                            "shannon": shannon(sub),
                            "total_reads": int(sub.sum())})
    alpha_df = pd.DataFrame(alpha).merge(smeta.reset_index(),
                                            on="sample", how="left")
    alpha_df.to_csv(OUT / "alpha_div_by_pool.tsv", sep="\t", index=False)
    pivot_alpha = (alpha_df.groupby(["pool", "compartment"])
                   [["shannon", "asv_count"]].median().round(2))
    print(pivot_alpha.to_string(), flush=True)

    # 5. Cosmopolitan enrichment
    print("\n[5] cosmopolitan enrichment ...", flush=True)
    ind_for_cosmo = ind[["asv_id", "pool", "emp_cosmo_90"]]
    cosmo = (ind_for_cosmo.groupby("pool")["emp_cosmo_90"]
             .agg(["count", "sum", "mean"])
             .rename(columns={"sum": "n_cosmo", "mean": "frac_cosmo"}))
    cosmo.to_csv(OUT / "cosmopolitan_enrichment.tsv", sep="\t")
    print(cosmo.round(3).to_string())

    # 3. Beta diversity within vs between sites for each pool
    # Use compartment-stratified, on rarefied to min n_reads per sample.
    print("\n[3] beta diversity (BC) within vs between sites per pool ...",
          flush=True)
    # Restrict to one compartment at a time and one trip (use trip with most
    # samples). For computational cost, downsample to top 200 ASVs per pool by
    # variance for the BC.
    bc_results = []
    for comp in ("rhizosphere", "surface", "deep"):
        for trip in (1, 2, 3, 4, 5):
            sub_meta = smeta[(smeta["compartment"] == comp) &
                              (smeta["trip"] == trip)].copy()
            sub_samples = list(sub_meta.index.intersection(ft.columns))
            if len(sub_samples) < 6: continue
            for pool, mask in (("relic", rel_idx), ("alive", ali_idx)):
                M = ft.loc[mask, sub_samples].T  # samples x ASVs
                # Drop empty samples
                row_sum = M.sum(axis=1)
                M = M[row_sum > 0]
                if len(M) < 6: continue
                Mn = M.div(M.sum(axis=1), axis=0).fillna(0)
                # BC distance
                D = pdist(Mn.values, metric="braycurtis")
                Ds = squareform(D)
                # within-site: same site_int -> within
                site_arr = sub_meta.loc[M.index, "site"].values
                within, between = [], []
                for i in range(len(M)):
                    for j in range(i + 1, len(M)):
                        if site_arr[i] == site_arr[j]:
                            within.append(Ds[i, j])
                        else:
                            between.append(Ds[i, j])
                bc_results.append({
                    "compartment": comp, "trip": trip, "pool": pool,
                    "n_samples": len(M),
                    "within_mean": float(np.mean(within)) if within else np.nan,
                    "within_n": len(within),
                    "between_mean": float(np.mean(between)) if between else np.nan,
                    "between_n": len(between),
                    "ratio_within_between": (np.mean(within) / np.mean(between)
                                               if (within and between and
                                                    np.mean(between) > 0) else np.nan),
                })
    bc_df = pd.DataFrame(bc_results)
    bc_df.to_csv(OUT / "beta_div_within_between.tsv", sep="\t", index=False)
    print(bc_df.round(3).head(20).to_string(index=False))

    # 7. Taxonomic decomposition
    print("\n[7] taxonomic decomposition ...", flush=True)
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index()
    tax = tax.rename(columns={"ASV": "asv_id"})
    m = ind.merge(tax, on="asv_id", how="left")
    for level in ("phylum", "class", "family", "genus"):
        if level not in m.columns: continue
        ct = (m.groupby([level, "pool"]).size().unstack(fill_value=0)
               .reindex(columns=["alive", "ambiguous", "relic"],
                          fill_value=0))
        ct["total"] = ct.sum(axis=1)
        ct = ct[ct["total"] >= 50]
        ct["frac_relic"] = ct["relic"] / ct["total"]
        ct["frac_alive"] = ct["alive"] / ct["total"]
        ct["log2_relic_alive"] = np.log2(
            (ct["relic"] + 1) / (ct["alive"] + 1))
        top_relic = ct.sort_values("log2_relic_alive", ascending=False).head(10)
        top_alive = ct.sort_values("log2_relic_alive").head(10)
        print(f"\n  Top 10 {level} enriched in RELIC:")
        print(top_relic[["alive", "relic", "total",
                            "frac_relic", "log2_relic_alive"]]
              .round(3).to_string())
        print(f"\n  Top 10 {level} enriched in ALIVE:")
        print(top_alive[["alive", "relic", "total",
                            "frac_alive", "log2_relic_alive"]]
              .round(3).to_string())
        ct.to_csv(OUT / f"taxonomy_{level}_relic_vs_alive.tsv",
                   sep="\t")

    # 8. CSP1-2 / Rubellimicrobium status
    print("\n[8] CSP1-2 / Rubellimicrobium status ...", flush=True)
    csp_fasta = CACHE / "csp1-2_asvs.fasta"
    if csp_fasta.exists():
        csp_ids = set()
        with open(csp_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    csp_ids.add(line[1:].strip().split()[0])
        csp_status = ind[ind["asv_id"].isin(csp_ids)].copy()
        csp_status.to_csv(OUT / "csp_rubelli_status.tsv",
                           sep="\t", index=False)
        print(f"  CSP1-2 ASVs found: {len(csp_status)}", flush=True)
        if len(csp_status):
            print(f"  pool breakdown:")
            print(csp_status["pool"].value_counts().to_string())
            print(f"  relic_score percentiles:")
            for q in (10, 25, 50, 75, 90):
                print(f"    p{q}: "
                      f"{np.percentile(csp_status['relic_score_full_gb'], q):.3f}")

    # Summary
    print("\nWriting summary ...", flush=True)
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Population-level relic analyses\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"ASV pool definitions:\n")
        fh.write(f"  relic_score >= {RELIC_HIGH}  -> 'relic'  "
                  f"(n={int((ind['pool']=='relic').sum())})\n")
        fh.write(f"  relic_score <= {ALIVE_HIGH}  -> 'alive'  "
                  f"(n={int((ind['pool']=='alive').sum())})\n")
        fh.write(f"  middle band             -> 'ambiguous'  "
                  f"(n={int((ind['pool']=='ambiguous').sum())})\n\n")

        fh.write("--- Per-sample fractions (medians) ---\n")
        fh.write(rec.groupby("compartment")[["relic_frac",
                                                  "alive_frac"]]
                  .median().round(3).to_string())
        fh.write("\n\n--- Per (compartment, trip) ---\n")
        fh.write(by_ct.round(3).to_string(index=False))

        fh.write("\n\n--- Per-site relic load summary ---\n")
        fh.write(per_site["median_relic_frac"].describe().round(3).to_string())

        fh.write("\n\n--- Cosmopolitan enrichment (EMP min25 90bp) ---\n")
        fh.write(cosmo.round(3).to_string())

        fh.write("\n\n--- Alpha diversity (per pool x compartment, medians) "
                  "---\n")
        fh.write(pivot_alpha.to_string())

        fh.write("\n\n--- Beta diversity (BC) within/between sites ---\n")
        fh.write(bc_df.round(3).to_string(index=False))
    print(f"Wrote {OUT}/summary.txt", flush=True)


if __name__ == "__main__":
    main()
