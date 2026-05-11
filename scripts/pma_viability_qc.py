#!/usr/bin/env python3
"""PMA per-pair viability QC.

Assesses whether within-ASV PMA signal is consistent enough across the 9 T/UT
pairs to support an ASV-level relic-likelihood indicator. If signal is mostly
noise, Path B is not viable.

Computes:
  - Per-(ASV, pair) T/UT ratios
  - Within-ASV consistency metrics (CV, range, IQR-of-log-ratio)
  - Between- vs within-ASV variance decomposition
  - Pair-pair correlation matrix on log ratios
  - Strata: rhizosphere (6 pairs at C1+C2) vs surface (3 pairs at C2)
  - Stability under thresholding (T/UT<0.1 = relic, >0.5 = alive)
  - Bimodality / mode count of viability distribution

Inputs:
  /home/leechuck/Public/software/empty-quarter/relic-dna/ASV_table_rel_abundance.tsv
  cache/test6_disconfirmation/pma_to_eq_match.tsv

Outputs:
  cache/test6_disconfirmation/qc_per_asv_pair_ratios.parquet
  cache/test6_disconfirmation/qc_per_asv_consistency.tsv
  cache/test6_disconfirmation/qc_pair_pair_corr.tsv
  cache/test6_disconfirmation/qc_summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
RELIC = Path("/home/leechuck/Public/software/empty-quarter/relic-dna")
OUT = CACHE / "test6_disconfirmation"

# Threshold for "detected" in a UT sample -- relative abundance > 1e-5
DETECT_FLOOR = 1e-5
# Min n pairs an ASV must be detected in to be assessed
MIN_PAIRS = 3
# eps for ratio computation
EPS = 1e-7
# ratio clip
RATIO_CLIP = 10.0


def parse_pma_sample(s: str):
    m = re.match(r"^(C[12])([RS])(\d+)(T|UT)$", s)
    if not m:
        return None
    return {"site": m.group(1), "comp_code": m.group(2),
            "rep": int(m.group(3)), "treatment": m.group(4)}


def main():
    print("Loading PMA relative abundance table ...", flush=True)
    rel = pd.read_csv(RELIC / "ASV_table_rel_abundance.tsv", sep="\t",
                      index_col=0)
    print(f"  PMA ASVs x samples: {rel.shape}", flush=True)

    samp = []
    for s in rel.columns:
        m = parse_pma_sample(s)
        if m:
            m["sample_orig"] = s
            samp.append(m)
    sm = pd.DataFrame(samp)
    pairs = (sm.pivot(index=["site", "comp_code", "rep"],
                      columns="treatment", values="sample_orig")
             .dropna(subset=["T", "UT"])
             .reset_index())
    pairs["pair_id"] = pairs.apply(
        lambda r: f"{r['site']}{r['comp_code']}{r['rep']}", axis=1)
    print(f"  T/UT pairs: {len(pairs)}", flush=True)
    print(f"  pair IDs: {list(pairs['pair_id'])}", flush=True)

    # Per (ASV, pair): T, UT, ratio
    print("\nComputing per-(ASV, pair) ratios ...", flush=True)
    rec = []
    for _, p in pairs.iterrows():
        T = rel[p["T"]]; UT = rel[p["UT"]]
        det_mask = UT > DETECT_FLOOR
        for asv_id in rel.index[det_mask]:
            t_val = float(T[asv_id]); ut_val = float(UT[asv_id])
            ratio = (t_val + EPS) / (ut_val + EPS)
            rec.append({
                "asv_id": asv_id,
                "pair_id": p["pair_id"],
                "site": p["site"],
                "comp_code": p["comp_code"],
                "T": t_val,
                "UT": ut_val,
                "ratio": min(ratio, RATIO_CLIP),
                "log_ratio": np.log10(min(ratio, RATIO_CLIP) + EPS),
            })
    pr = pd.DataFrame(rec)
    print(f"  total per-ASV-pair detections: {len(pr)}", flush=True)
    pr.to_parquet(OUT / "qc_per_asv_pair_ratios.parquet")

    # Per-ASV consistency metrics (only ASVs detected in >=MIN_PAIRS)
    print(f"\nComputing per-ASV consistency (min {MIN_PAIRS} pairs) ...",
          flush=True)
    g = pr.groupby("asv_id")
    cons = g.agg(n_pairs=("pair_id", "count"),
                 mean_ratio=("ratio", "mean"),
                 median_ratio=("ratio", "median"),
                 std_ratio=("ratio", "std"),
                 mean_log=("log_ratio", "mean"),
                 std_log=("log_ratio", "std"),
                 mean_UT=("UT", "mean")).reset_index()
    cons["cv_ratio"] = cons["std_ratio"] / cons["mean_ratio"].replace(0, np.nan)
    # IQR of log_ratio
    iqr = g["log_ratio"].agg(lambda x: float(np.percentile(x, 75) -
                                              np.percentile(x, 25)))
    cons = cons.merge(iqr.rename("iqr_log").reset_index(), on="asv_id")
    cons_q = cons[cons["n_pairs"] >= MIN_PAIRS].copy()
    print(f"  ASVs in >=3 pairs: {len(cons_q)}", flush=True)
    print(f"  ASVs in >=6 pairs: {(cons_q['n_pairs']>=6).sum()}", flush=True)
    print(f"  ASVs in 9 pairs:    {(cons_q['n_pairs']==9).sum()}", flush=True)
    cons_q.to_csv(OUT / "qc_per_asv_consistency.tsv", sep="\t", index=False)

    # Variance decomposition
    print("\nVariance decomposition (log_ratio): ...", flush=True)
    pr_q = pr[pr["asv_id"].isin(cons_q["asv_id"])]
    grand_mean = pr_q["log_ratio"].mean()
    asv_means = pr_q.groupby("asv_id")["log_ratio"].mean()
    n_per_asv = pr_q.groupby("asv_id")["log_ratio"].count()
    ss_between = float(((asv_means - grand_mean) ** 2 * n_per_asv).sum())
    ss_within = float((pr_q.merge(
        asv_means.rename("asv_mean").reset_index(), on="asv_id")
        .assign(d=lambda x: (x["log_ratio"] - x["asv_mean"]) ** 2)["d"]
        .sum()))
    n_total = len(pr_q); k = len(asv_means)
    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n_total - k)
    f_ratio = ms_between / ms_within if ms_within > 0 else float("inf")
    icc = (ms_between - ms_within) / (
        ms_between + (n_per_asv.mean() - 1) * ms_within)
    print(f"  SS_between={ss_between:.2f}  SS_within={ss_within:.2f}",
          flush=True)
    print(f"  MS_between={ms_between:.4f}  MS_within={ms_within:.4f}",
          flush=True)
    print(f"  F = {f_ratio:.2f}", flush=True)
    print(f"  ICC(1,k) ≈ {icc:.3f}    (1=perfect; 0=noise)", flush=True)

    # Pair-pair correlation matrix
    print("\nPair-pair correlation on log_ratio (ASVs detected in both) ...",
          flush=True)
    pivot = (pr.pivot_table(index="asv_id", columns="pair_id",
                              values="log_ratio")
             .dropna(thresh=2))
    pair_ids = list(pivot.columns)
    rho_rows = []
    for i, pi in enumerate(pair_ids):
        for j, pj in enumerate(pair_ids):
            if j <= i: continue
            both = pivot[[pi, pj]].dropna()
            if len(both) < 30:
                continue
            r, p = spearmanr(both[pi], both[pj])
            rho_rows.append({"pair_a": pi, "pair_b": pj, "n": len(both),
                              "rho": float(r), "p": float(p)})
    rho_df = pd.DataFrame(rho_rows)
    rho_df.to_csv(OUT / "qc_pair_pair_corr.tsv", sep="\t", index=False)
    print(f"  median pair-pair Spearman rho: {rho_df['rho'].median():.3f}",
          flush=True)
    print(f"  range: [{rho_df['rho'].min():.3f}, "
          f"{rho_df['rho'].max():.3f}]", flush=True)

    # Strata: rhizosphere vs surface
    print("\nStrata: rhizosphere (6 pairs) vs surface (3 pairs) ...",
          flush=True)
    rhizo_pairs = [p for p in pair_ids if p[2] == "R"]
    surf_pairs = [p for p in pair_ids if p[2] == "S"]
    print(f"  rhizo pair IDs: {rhizo_pairs}")
    print(f"  surf  pair IDs: {surf_pairs}")
    # Per-ASV: median in rhizosphere vs median in surface
    rhizo_med = (pr[pr["pair_id"].isin(rhizo_pairs)]
                 .groupby("asv_id")["log_ratio"].median())
    surf_med = (pr[pr["pair_id"].isin(surf_pairs)]
                .groupby("asv_id")["log_ratio"].median())
    common = rhizo_med.index.intersection(surf_med.index)
    if len(common) >= 30:
        r, p = spearmanr(rhizo_med.loc[common], surf_med.loc[common])
        print(f"  rhizo vs surf median log_ratio:  rho={r:.3f}  p={p:.3g}  "
              f"n={len(common)}", flush=True)

    # Site stratum: C1 vs C2 (rhizosphere only -- both sites have R)
    c1_pairs = [p for p in rhizo_pairs if p[:2] == "C1"]
    c2_pairs = [p for p in rhizo_pairs if p[:2] == "C2"]
    c1_med = (pr[pr["pair_id"].isin(c1_pairs)]
              .groupby("asv_id")["log_ratio"].median())
    c2_med = (pr[pr["pair_id"].isin(c2_pairs)]
              .groupby("asv_id")["log_ratio"].median())
    common = c1_med.index.intersection(c2_med.index)
    if len(common) >= 30:
        r, p = spearmanr(c1_med.loc[common], c2_med.loc[common])
        print(f"  C1R vs C2R median log_ratio:     rho={r:.3f}  p={p:.3g}  "
              f"n={len(common)}", flush=True)
    c1_c2_rho = float(r) if len(common) >= 30 else None

    # Stability of binary classification under thresholding
    print("\nClassification stability "
          "(relic if T/UT<0.1, alive if T/UT>0.5 per pair) ...", flush=True)
    pr["cls"] = np.where(pr["ratio"] < 0.1, "relic",
                          np.where(pr["ratio"] > 0.5, "alive", "ambig"))
    cls = (pr.groupby(["asv_id", "cls"]).size().unstack(fill_value=0)
           .reindex(columns=["alive", "ambig", "relic"], fill_value=0))
    cls["n_total"] = cls.sum(axis=1)
    cls = cls[cls["n_total"] >= 3]
    cls["majority"] = cls[["alive", "ambig", "relic"]].idxmax(axis=1)
    cls["majority_frac"] = cls[["alive", "ambig", "relic"]].max(axis=1) \
        / cls["n_total"]
    cls["stable"] = cls["majority_frac"] >= 0.66
    print(f"  ASVs assessed: {len(cls)}", flush=True)
    print(f"  fraction stably classified (>=66% of pairs agree): "
          f"{cls['stable'].mean():.3f}", flush=True)
    print(f"  majority class breakdown:")
    print(cls["majority"].value_counts().to_string())
    print(f"  among stable ASVs:")
    print(cls.loc[cls["stable"], "majority"].value_counts().to_string())

    # Bimodality / distribution shape of mean log_ratio
    mean_log = cons_q["mean_log"]
    p10 = float(np.percentile(mean_log, 10))
    p25 = float(np.percentile(mean_log, 25))
    p50 = float(np.percentile(mean_log, 50))
    p75 = float(np.percentile(mean_log, 75))
    p90 = float(np.percentile(mean_log, 90))
    print("\nMean log_ratio distribution per ASV (n>=3 pairs):")
    print(f"  p10={p10:.2f}  p25={p25:.2f}  p50={p50:.2f}  "
          f"p75={p75:.2f}  p90={p90:.2f}", flush=True)
    print(f"  fraction with mean_log < -1 (T/UT<0.1, relic-like): "
          f"{(mean_log < -1).mean():.3f}", flush=True)
    print(f"  fraction with mean_log > -0.3 (T/UT>0.5, alive-like): "
          f"{(mean_log > -0.3).mean():.3f}", flush=True)
    print(f"  fraction in middle band [-1, -0.3]: "
          f"{((mean_log >= -1) & (mean_log <= -0.3)).mean():.3f}", flush=True)

    # Write summary
    with open(OUT / "qc_summary.txt", "w") as fh:
        fh.write("PMA per-pair viability QC\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Pairs: {len(pairs)} ({sum(p[2]=='R' for p in pair_ids)} "
                  f"rhizosphere + {sum(p[2]=='S' for p in pair_ids)} surface)\n")
        fh.write(f"Pair IDs: {list(pairs['pair_id'])}\n\n")
        fh.write(f"Total per-ASV-pair detections (UT > {DETECT_FLOOR}): "
                  f"{len(pr)}\n")
        fh.write(f"ASVs detected in >=3 pairs: {len(cons_q)}\n")
        fh.write(f"ASVs detected in all 9 pairs: "
                  f"{(cons_q['n_pairs']==9).sum()}\n\n")

        fh.write("--- Variance decomposition (on log10 T/UT ratio) ---\n")
        fh.write(f"  SS_between (ASV identity)  = {ss_between:.2f}\n")
        fh.write(f"  SS_within  (pair-pair noise) = {ss_within:.2f}\n")
        fh.write(f"  MS_between = {ms_between:.4f}\n")
        fh.write(f"  MS_within  = {ms_within:.4f}\n")
        fh.write(f"  F = MS_b / MS_w = {f_ratio:.2f}\n")
        fh.write(f"  ICC(1,k) ≈ {icc:.3f}  "
                  "[1=perfect; 0=noise; >0.5 = workable]\n\n")

        fh.write("--- Pair-pair correlation (Spearman, log_ratio) ---\n")
        fh.write(f"  median rho across {len(rho_df)} pair-pair comparisons: "
                  f"{rho_df['rho'].median():.3f}\n")
        fh.write(f"  range:  [{rho_df['rho'].min():.3f}, "
                  f"{rho_df['rho'].max():.3f}]\n")
        if c1_c2_rho is not None:
            fh.write(f"  C1R vs C2R median log_ratio rho: {c1_c2_rho:.3f}\n")
        fh.write("\n")

        fh.write("--- Classification stability ---\n")
        fh.write(f"  ASVs in >=3 pairs: {len(cls)}\n")
        fh.write(f"  Stably classified (>=66% same label): "
                  f"{cls['stable'].mean():.3f}\n")
        fh.write(f"  Majority class breakdown:\n"
                  f"{cls['majority'].value_counts().to_string()}\n")
        fh.write(f"  Stable subset:\n"
                  f"{cls.loc[cls['stable'], 'majority'].value_counts().to_string()}\n\n")

        fh.write("--- Distribution shape ---\n")
        fh.write(f"  mean_log_ratio percentiles: p10={p10:.2f}  p25={p25:.2f}  "
                  f"p50={p50:.2f}  p75={p75:.2f}  p90={p90:.2f}\n")
        fh.write(f"  relic-like (T/UT<0.1):  "
                  f"{(mean_log < -1).mean():.3f}\n")
        fh.write(f"  alive-like (T/UT>0.5):  "
                  f"{(mean_log > -0.3).mean():.3f}\n")
        fh.write(f"  ambiguous middle band: "
                  f"{((mean_log >= -1) & (mean_log <= -0.3)).mean():.3f}\n\n")

        fh.write("--- VERDICT ---\n")
        if icc > 0.5 and rho_df['rho'].median() > 0.4:
            fh.write("  Signal is workable: between-ASV variance dominates "
                      "within-ASV noise; pair-pair concordance is moderate-"
                      "to-strong. PMA can serve as one axis of a composite "
                      "indicator. Path B is feasible.\n")
        elif icc > 0.3 and rho_df['rho'].median() > 0.25:
            fh.write("  Signal is BORDERLINE workable: between-ASV variance "
                      "exceeds noise but only modestly. PMA can contribute as "
                      "a noisy axis in a composite, NOT as a standalone "
                      "ASV-level classifier.\n")
        else:
            fh.write("  Signal is too noisy at ASV level: within-pair noise "
                      "dominates between-ASV variance. PMA can ONLY support "
                      "population-level fractional claims (e.g., 'X% of "
                      "DNA at C1+C2 is relic-like') and NOT per-ASV "
                      "indicator construction.\n")
    print(f"\nWrote {OUT}/qc_summary.txt", flush=True)


if __name__ == "__main__":
    main()
