#!/usr/bin/env python3
"""Reconstruct PICRUSt2 metagenome predictions per sample restricted to the
alive pool (relic_score < 0.3) and compare to the standard all-ASV PICRUSt2.

Steps:
  1. Load per-ASV KO predictions (PICRUSt2 KO_predicted.tsv).
  2. Load feature_table (ASV x sample).
  3. Define pool membership from relic_indicator.
  4. For each sample: aggregate KO copies x ASV reads over (alive | relic | all).
  5. Compare:
       - per-sample functional richness (n KOs detected) per pool
       - per-sample functional Shannon per pool
       - top KOs differentially enriched between alive vs all
  6. Repeat Test 5 osmolyte uptake/biosynth on alive-only.

Outputs:
  cache/relic_population/alive_only_metagenome_pred.parquet
  cache/relic_population/alive_relic_func_diversity.tsv
  cache/relic_population/alive_vs_all_top_KOs.tsv
  cache/relic_population/test5_osmolyte_alive.tsv
  cache/relic_population/alive_only_functions_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
OUT = CACHE / "relic_population"
PIC = Path("/home/leechuck/Public/software/empty-quarter/data/processed/"
            "functional/picrust2")

RELIC_HIGH = 0.7
ALIVE_HIGH = 0.3


def main():
    print("Loading indicator and pool labels ...", flush=True)
    ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                       "relic_indicator_with_damage_per_asv.tsv", sep="\t",
                       usecols=["asv_id", "relic_score_full_gb"])
    ind["pool"] = np.where(ind["relic_score_full_gb"] >= RELIC_HIGH, "relic",
                              np.where(ind["relic_score_full_gb"] <= ALIVE_HIGH,
                                        "alive", "ambiguous"))
    relic_set = set(ind.loc[ind["pool"] == "relic", "asv_id"])
    alive_set = set(ind.loc[ind["pool"] == "alive", "asv_id"])

    print("Loading feature table ...", flush=True)
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    print(f"  {ft.shape}", flush=True)

    print("Loading per-ASV KO predictions (large) ...", flush=True)
    # KO_predicted is 4.6GB; load with float32 to fit
    ko_pred = pd.read_csv(PIC / "KO_predicted.tsv", sep="\t",
                            index_col=0, dtype={"sequence": str})
    print(f"  {ko_pred.shape}  dtype={ko_pred.dtypes.iloc[0]}", flush=True)
    # Coerce to float32 to save mem
    ko_pred = ko_pred.astype(np.float32)

    # Restrict to ASVs present in both
    common = ft.index.intersection(ko_pred.index)
    print(f"  ASVs in BOTH ft and KO_pred: {len(common)}", flush=True)
    ft = ft.loc[common]
    ko_pred = ko_pred.loc[common]

    # Build alive/relic feature subtables
    ft_alive = ft.loc[ft.index.isin(alive_set)]
    ft_relic = ft.loc[ft.index.isin(relic_set)]
    print(f"  alive ASVs in ft: {len(ft_alive)}", flush=True)
    print(f"  relic ASVs in ft: {len(ft_relic)}", flush=True)

    # KO predictions for the pools
    ko_alive = ko_pred.loc[ko_pred.index.isin(alive_set)]
    ko_relic = ko_pred.loc[ko_pred.index.isin(relic_set)]

    # Per-sample metagenome prediction: KO copies x reads, summed over ASVs
    print("\nComputing per-sample metagenome predictions per pool ...", flush=True)
    # All-ASV
    print("  all ...", flush=True)
    M_all = ko_pred.T.dot(ft.values)  # (n_KO, n_samples)? actually .dot wants
    # Better: align ASVs along same axis
    # ft: ASV x sample. ko_pred: ASV x KO. We want KO x sample = ko_pred.T x ft
    M_all = ko_pred.T.values.dot(ft.values)  # (n_KO, n_samples)
    print(f"    M_all shape: {M_all.shape}", flush=True)

    print("  alive ...", flush=True)
    M_alive = ko_alive.T.values.dot(ft_alive.values)
    print(f"    M_alive shape: {M_alive.shape}", flush=True)

    print("  relic ...", flush=True)
    M_relic = ko_relic.T.values.dot(ft_relic.values)
    print(f"    M_relic shape: {M_relic.shape}", flush=True)

    sample_cols = list(ft.columns)
    ko_index = list(ko_pred.columns)

    M_all_df = pd.DataFrame(M_all, index=ko_index, columns=sample_cols)
    M_alive_df = pd.DataFrame(M_alive, index=ko_index, columns=sample_cols)
    M_relic_df = pd.DataFrame(M_relic, index=ko_index, columns=sample_cols)
    M_alive_df.to_parquet(OUT / "alive_only_metagenome_pred.parquet")
    M_relic_df.to_parquet(OUT / "relic_only_metagenome_pred.parquet")

    # Per-sample functional diversity per pool
    print("\nFunctional diversity per pool per sample ...", flush=True)
    def shannon(col):
        x = col[col > 0]
        if len(x) == 0: return 0.0
        p = x / x.sum()
        return float(-(p * np.log(p)).sum())

    rec = []
    for s in sample_cols:
        for label, M in (("all", M_all_df[s]), ("alive", M_alive_df[s]),
                            ("relic", M_relic_df[s])):
            rec.append({"sample": s, "pool": label,
                          "ko_richness": int((M > 0).sum()),
                          "ko_shannon": shannon(M.values),
                          "ko_total": float(M.sum())})
    div_df = pd.DataFrame(rec)
    smeta = parse_samples_to_df(sample_cols).set_index("sample")
    div_df = div_df.merge(smeta.reset_index(), on="sample", how="left")
    div_df.to_csv(OUT / "alive_relic_func_diversity.tsv", sep="\t", index=False)

    print("\n  per-pool x compartment medians:")
    print(div_df.groupby(["pool", "compartment"])
          [["ko_richness", "ko_shannon"]].median().round(2).to_string())

    # Per-sample relabund per pool, then compute log2(alive_relabund / all_relabund)
    print("\nDifferentially abundant KOs (alive vs all) ...", flush=True)
    M_all_rel = M_all_df.div(M_all_df.sum(axis=0).replace(0, 1), axis=1)
    M_alive_rel = M_alive_df.div(M_alive_df.sum(axis=0).replace(0, 1), axis=1)
    M_relic_rel = M_relic_df.div(M_relic_df.sum(axis=0).replace(0, 1), axis=1)

    # Mean relabund per pool across samples
    mean_all = M_all_rel.mean(axis=1)
    mean_alive = M_alive_rel.mean(axis=1)
    mean_relic = M_relic_rel.mean(axis=1)
    diff = pd.DataFrame({
        "ko": ko_index,
        "mean_relabund_all": mean_all.values,
        "mean_relabund_alive": mean_alive.values,
        "mean_relabund_relic": mean_relic.values,
    })
    diff["log2_alive_over_all"] = np.log2(
        (diff["mean_relabund_alive"] + 1e-12) /
        (diff["mean_relabund_all"] + 1e-12))
    diff["log2_alive_over_relic"] = np.log2(
        (diff["mean_relabund_alive"] + 1e-12) /
        (diff["mean_relabund_relic"] + 1e-12))
    diff = diff.sort_values("log2_alive_over_relic", ascending=False)
    diff.to_csv(OUT / "alive_vs_all_top_KOs.tsv", sep="\t", index=False)

    # Show top 15 enriched in alive vs relic
    top_a = diff.head(15)
    top_r = diff.tail(15)
    print("\nTop 15 KOs ENRICHED in alive vs relic:")
    print(top_a[["ko", "mean_relabund_alive", "mean_relabund_relic",
                  "log2_alive_over_relic"]].round(5).to_string(index=False))
    print("\nTop 15 KOs ENRICHED in relic vs alive:")
    print(top_r[["ko", "mean_relabund_alive", "mean_relabund_relic",
                  "log2_alive_over_relic"]].round(5).to_string(index=False))

    # Test 5: osmolyte uptake (BetA, ProP, OpuA-D, etc.) vs biosynthesis
    # Use known KO IDs. From Test 5 memory: uptake/biosynth ratio was ~230x in
    # METAGENOMES. Re-test with the alive-only PICRUSt2-style aggregation.
    print("\n=== Test 5 reanalysis on alive subset ===", flush=True)
    # Glycine betaine uptake: K02000-3 (proU/proV), K05845, K05846, K05847,
    #                          K05848 (opuA-D)
    # Glycine betaine biosynth (from choline): K00108 (betA), K00130 (betB)
    # Simplified set from Test 5 (verify by re-using the same KOs)
    uptake_kos = ["K02000", "K02001", "K02002", "K02003",
                    "K05845", "K05846", "K05847", "K05848",
                    "K05874", "K05875", "K05876", "K05877"]
    biosynth_kos = ["K00108", "K00130", "K17755"]
    upt_present = [k for k in uptake_kos if k in M_all_df.index]
    bio_present = [k for k in biosynth_kos if k in M_all_df.index]
    print(f"  uptake KOs found: {len(upt_present)}/{len(uptake_kos)}")
    print(f"  biosynth KOs found: {len(bio_present)}/{len(biosynth_kos)}")

    rows = []
    for label, M in (("all", M_all_df), ("alive", M_alive_df),
                       ("relic", M_relic_df)):
        upt = M.loc[upt_present].sum(axis=0)
        bio = M.loc[bio_present].sum(axis=0)
        ratio = (upt / bio.replace(0, np.nan))
        rows.append({"pool": label,
                      "median_uptake_total": float(upt.median()),
                      "median_biosynth_total": float(bio.median()),
                      "median_ratio": float(ratio.median()),
                      "mean_ratio": float(ratio.mean()),
                      "n_samples_with_ratio": int(ratio.notna().sum())})
    test5 = pd.DataFrame(rows)
    test5.to_csv(OUT / "test5_osmolyte_alive.tsv", sep="\t", index=False)
    print(test5.round(2).to_string(index=False))

    # Write summary
    with open(OUT / "alive_only_functions_summary.txt", "w") as fh:
        fh.write("Alive-only functional reanalysis\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"ASVs in feature table & KO predictions: {len(common)}\n")
        fh.write(f"Alive ASVs: {len(ft_alive)}\n")
        fh.write(f"Relic ASVs: {len(ft_relic)}\n\n")

        fh.write("--- Per-pool x compartment functional diversity (medians) ---\n")
        fh.write(div_df.groupby(["pool", "compartment"])
                  [["ko_richness", "ko_shannon"]].median().round(2)
                  .to_string())

        fh.write("\n\n--- Top 15 KOs enriched in alive vs relic ---\n")
        fh.write(top_a[["ko", "mean_relabund_alive", "mean_relabund_relic",
                          "log2_alive_over_relic"]].round(5)
                  .to_string(index=False))

        fh.write("\n\n--- Top 15 KOs enriched in relic vs alive ---\n")
        fh.write(top_r[["ko", "mean_relabund_alive", "mean_relabund_relic",
                          "log2_alive_over_relic"]].round(5)
                  .to_string(index=False))

        fh.write("\n\n--- Test 5 osmolyte uptake/biosynth ratio reanalysis ---\n")
        fh.write(test5.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
