#!/usr/bin/env python3
"""Add MAG-presence prior to the relic-likelihood indicator.

Logic: an ASV whose 16S sequence matches a MAG-derived 16S at >=97%
identity is mechanistically incompatible with being purely relic — relic
DNA fragments don't co-assemble into 4-Mb genomes with consistent
coverage. So MAG match -> strong alive prior.

Combines with prior-augmented score (relic_indicator_with_priors.py
output) via Bayesian log-odds:
  logit(P(relic | data, MAG)) = logit(P_post) + delta_MAG

where delta_MAG depends on MAG-match strength:
  has_match=False              -> delta_MAG = 0      (no info)
  has_match (>=97% pid, 1 MAG) -> delta_MAG = -1.5   (alive bias)
  has_match (>=99% pid, multi) -> delta_MAG = -2.5   (strong alive bias)

Outputs:
  cache/relic_priors/asv_mag_matches.tsv
  cache/relic_priors/relic_score_with_mag_prior.tsv
  cache/relic_priors/csp_with_mag_prior.tsv
  cache/relic_priors/mag_prior_summary.txt
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "relic_priors"


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def main():
    print("Loading vsearch matches ...", flush=True)
    cols = ["asv_id", "mag_id", "pid", "alen", "mm", "gaps", "qs", "qe",
            "ts", "te", "evalue", "bs"]
    m = pd.read_csv(CACHE / "mag_16s" / "asv_to_mag_match.tsv",
                     sep="\t", header=None, names=cols)
    print(f"  raw match records: {len(m)}", flush=True)
    print(f"  unique ASVs matching: {m['asv_id'].nunique()}", flush=True)

    # Per-ASV: n_MAG_matches, max_pid, n_high_pid (>=99%), n_strong_match (>=99% AND alen>=200)
    agg = (m.groupby("asv_id")
           .agg(n_mag_matches=("mag_id", "nunique"),
                max_pid=("pid", "max"),
                mean_pid=("pid", "mean"),
                max_alen=("alen", "max"))
           .reset_index())
    agg["n_high_pid"] = (m[m["pid"] >= 99]
                          .groupby("asv_id")["mag_id"].nunique()
                          .reindex(agg["asv_id"]).fillna(0).astype(int).values)
    agg["has_match"] = 1
    agg["has_strong_match"] = ((agg["max_pid"] >= 99) &
                                  (agg["max_alen"] >= 200)).astype(int)
    agg.to_csv(OUT / "asv_mag_matches.tsv", sep="\t", index=False)
    print(f"  ASVs with high-pid match (>=99%): "
          f"{(agg['n_high_pid'] > 0).sum()}", flush=True)
    print(f"  ASVs with strong match (>=99%, >=200bp): "
          f"{agg['has_strong_match'].sum()}", flush=True)

    # Load posterior + add MAG features
    print("\nMerging with prior-augmented relic score ...", flush=True)
    p = pd.read_csv(OUT / "relic_score_with_priors.tsv", sep="\t")
    p = p.merge(agg, on="asv_id", how="left")
    p["has_match"] = p["has_match"].fillna(0).astype(int)
    p["has_strong_match"] = p["has_strong_match"].fillna(0).astype(int)
    p["n_mag_matches"] = p["n_mag_matches"].fillna(0).astype(int)
    p["max_pid"] = p["max_pid"].fillna(0)

    # Define delta_MAG (log-odds shift toward alive)
    def mag_delta(row):
        if row["has_strong_match"]:
            # Strong evidence ASV is from a successfully-assembled MAG
            return -2.5
        elif row["has_match"]:
            # Weak match (>=97%, possibly short alignment)
            return -1.5
        return 0.0

    p["delta_mag"] = p.apply(mag_delta, axis=1)
    p["log_post_with_mag"] = (logit(p["relic_score_posterior"])
                                  + p["delta_mag"])
    p["relic_score_with_mag"] = sigmoid(p["log_post_with_mag"])

    # Save
    p.to_csv(OUT / "relic_score_with_mag_prior.tsv", sep="\t",
              index=False)

    # Summary
    print("\n=== Pool counts after MAG prior ===", flush=True)
    for label, col in (("model_only", "relic_score_full_gb"),
                            ("posterior_no_MAG", "relic_score_posterior"),
                            ("posterior_+_MAG", "relic_score_with_mag")):
        n_alive = (p[col] <= 0.3).sum()
        n_amb = ((p[col] > 0.3) & (p[col] < 0.7)).sum()
        n_relic = (p[col] >= 0.7).sum()
        print(f"  {label:<20}  alive={n_alive:>6}  ambig={n_amb:>6}  "
              f"relic={n_relic:>6}", flush=True)

    print("\n=== Score distribution comparison ===", flush=True)
    for q in (10, 25, 50, 75, 90):
        print(f"  p{q:>2}  model: "
              f"{np.percentile(p['relic_score_full_gb'], q):.3f}    "
              f"post(prior): "
              f"{np.percentile(p['relic_score_posterior'], q):.3f}    "
              f"post(+MAG): "
              f"{np.percentile(p['relic_score_with_mag'], q):.3f}",
              flush=True)

    # CSP1-2 check
    csp_fasta = CACHE / "csp1-2_asvs.fasta"
    csp_ids = set()
    with open(csp_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                csp_ids.add(line[1:].strip().split()[0])
    csp = p[p["asv_id"].isin(csp_ids)].copy()
    print(f"\n=== CSP1-2 status with MAG prior ===", flush=True)
    print(f"  CSP1-2 ASVs: {len(csp)}", flush=True)
    print(f"  CSP1-2 ASVs with ANY MAG match (>=97% pid): "
          f"{int(csp['has_match'].sum())}", flush=True)
    print(f"  CSP1-2 ASVs with STRONG MAG match (>=99%, >=200bp): "
          f"{int(csp['has_strong_match'].sum())}", flush=True)
    print(f"  CSP1-2 max_pid distribution:")
    print(f"    median: {csp['max_pid'].median():.2f}")
    print(f"    range:  [{csp['max_pid'].min():.2f}, "
          f"{csp['max_pid'].max():.2f}]")
    print(csp[["asv_id", "genus", "has_match", "max_pid", "n_mag_matches",
                  "relic_score_full_gb", "relic_score_posterior",
                  "relic_score_with_mag"]].round(3).to_string(index=False),
          flush=True)
    print(f"\n  CSP1-2 median relic_score:")
    print(f"    model only:        "
          f"{csp['relic_score_full_gb'].median():.3f}", flush=True)
    print(f"    +taxonomy prior:   "
          f"{csp['relic_score_posterior'].median():.3f}", flush=True)
    print(f"    +MAG prior:        "
          f"{csp['relic_score_with_mag'].median():.3f}", flush=True)
    csp.to_csv(OUT / "csp_with_mag_prior.tsv", sep="\t", index=False)

    # Top 10 alive genera in MAG-augmented pool
    print(f"\n=== Top 15 genera in alive subset (MAG-augmented) ===",
          flush=True)
    alive_mag = p[p["relic_score_with_mag"] <= 0.3]
    print(alive_mag["genus"].value_counts().head(15).to_string())

    # MAG-only "definitely alive" set
    print(f"\n=== ASVs with strong MAG match (definitively alive) ===",
          flush=True)
    strong = p[p["has_strong_match"] == 1]
    print(f"  Total: {len(strong)}", flush=True)
    print(f"  Top 15 genera:")
    print(strong["genus"].value_counts().head(15).to_string())
    print(f"  Phyla distribution:")
    print(strong["phylum"].value_counts().head(10).to_string())

    # How does the model rate them?
    print(f"\n  Of these {len(strong)} 'definitively alive' ASVs:")
    print(f"    model relic_score median: "
          f"{strong['relic_score_full_gb'].median():.3f}", flush=True)
    print(f"    n with model_relic_score >= 0.7 (would have been called "
          f"relic without MAG): {(strong['relic_score_full_gb'] >= 0.7).sum()}",
          flush=True)
    print(f"    n with model_relic_score <= 0.3 (model agrees alive): "
          f"{(strong['relic_score_full_gb'] <= 0.3).sum()}", flush=True)

    # Pool transitions
    p["pool_post"] = np.where(p["relic_score_posterior"] <= 0.3, "alive",
                                np.where(p["relic_score_posterior"] >= 0.7,
                                          "relic", "ambig"))
    p["pool_with_mag"] = np.where(p["relic_score_with_mag"] <= 0.3, "alive",
                                       np.where(p["relic_score_with_mag"] >= 0.7,
                                                 "relic", "ambig"))
    print(f"\n=== Pool transitions (posterior_prior -> +MAG) ===",
          flush=True)
    print(pd.crosstab(p["pool_post"], p["pool_with_mag"], margins=True)
          .to_string())

    # Save summary
    with open(OUT / "mag_prior_summary.txt", "w") as fh:
        fh.write("MAG-presence prior summary\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Total ASVs: {len(p)}\n")
        fh.write(f"ASVs with any MAG match (>=97%): "
                  f"{int(p['has_match'].sum())} "
                  f"({p['has_match'].mean()*100:.1f}%)\n")
        fh.write(f"ASVs with strong MAG match (>=99%, >=200bp): "
                  f"{int(p['has_strong_match'].sum())} "
                  f"({p['has_strong_match'].mean()*100:.1f}%)\n\n")

        fh.write("CSP1-2 (24 ASVs):\n")
        fh.write(f"  with MAG match: {int(csp['has_match'].sum())}\n")
        fh.write(f"  with strong MAG match: "
                  f"{int(csp['has_strong_match'].sum())}\n")
        fh.write(f"  max_pid range: [{csp['max_pid'].min():.2f}, "
                  f"{csp['max_pid'].max():.2f}]\n")
        fh.write(f"  median relic_score (model): "
                  f"{csp['relic_score_full_gb'].median():.3f}\n")
        fh.write(f"  median relic_score (post): "
                  f"{csp['relic_score_posterior'].median():.3f}\n")
        fh.write(f"  median relic_score (+MAG): "
                  f"{csp['relic_score_with_mag'].median():.3f}\n\n")

        fh.write("Pool counts:\n")
        for label, col in (("model_only", "relic_score_full_gb"),
                                ("post_prior", "relic_score_posterior"),
                                ("post_+_MAG", "relic_score_with_mag")):
            n_alive = (p[col] <= 0.3).sum()
            n_relic = (p[col] >= 0.7).sum()
            fh.write(f"  {label}: alive={n_alive}, relic={n_relic}\n")
    print(f"\nWrote {OUT}/mag_prior_summary.txt", flush=True)


if __name__ == "__main__":
    main()
