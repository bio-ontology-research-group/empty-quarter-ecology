#!/usr/bin/env python3
"""Bayesian augmentation of the relic-likelihood indicator with
fundamental biological priors.

Combines:
  log_odds(P(relic|data)) = log_odds_prior(taxonomy + ecology)
                           + log_odds_evidence(model output)

Priors (ordered by confidence):
  Tier 1: Habitat exclusions    (strong)
  Tier 2: Extremophile / arid specialists (alive)
  Tier 3: PMA cell-biology corrections (down-weight Bacilli alive bias)
  Tier 5: Phylum baseline       (weak)

Outputs:
  cache/relic_priors/per_asv_priors.tsv
  cache/relic_priors/relic_score_with_priors.tsv
  cache/relic_priors/csp_with_priors.tsv
  cache/relic_priors/comparison_with_no_prior.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "relic_priors"
OUT.mkdir(parents=True, exist_ok=True)

# Tier 1: Habitat-exclusion priors (strong relic prior)
# Genera physiologically impossible in oxic, hyperarid, low-organic
# soil environment.
T1_RELIC_PRIORS = {
    # Obligate marine
    "Pelagibacter": 0.97, "Pelagibacterales": 0.97,
    "Prochlorococcus": 0.97, "Synechococcus": 0.85,
    # marine clades not really in soil
    "Pseudoalteromonas": 0.92, "Vibrio": 0.92,
    "Alteromonas": 0.92, "Marinobacterium": 0.92,
    # Obligate methanogens
    "Methanosarcina": 0.95, "Methanobrevibacter": 0.95,
    "Methanobacterium": 0.95, "Methanocaldococcus": 0.95,
    "Methanocorpusculum": 0.95, "Methanocella": 0.95,
    # Obligate anaerobes typically
    "Clostridioides": 0.88,  # spore-formers but obligate anaerobes
    # Mammalian gut obligates
    "Bacteroides": 0.85,  # in soil = unusual unless atypical
    "Bifidobacterium": 0.85,
    "Lactobacillus": 0.80,
    # Animal pathogens
    "Mycoplasma": 0.90,
    # Strict obligate intracellular
    "Wolbachia": 0.95, "Rickettsia": 0.92,
}

# Tier 2: Extremophile / arid specialist alive priors (strong alive)
T2_ALIVE_PRIORS = {
    # Halophilic Bacilli & Pseudomonadota
    "Halomonas": 0.15, "Halothermothrix": 0.18,
    "Halobacillus": 0.18, "Halothiobacillus": 0.20,
    "Salimicrobium": 0.20, "Salinimicrobium": 0.18,
    "Salinibacter": 0.20, "Salinibacillus": 0.20,
    "Salinispora": 0.20, "Salinicoccus": 0.20,
    # Halophilic archaea
    "Halobacterium": 0.15, "Halococcus": 0.18,
    "Haloarcula": 0.18, "Haloferax": 0.18,
    # Drought / desiccation tolerant Bacilli
    "Bacillus": 0.30,  # weak alive (spore-formers ubiquitous in soil)
    "Aureibacillus": 0.25, "Lysinibacillus": 0.25,
    "Sediminibacillus": 0.25, "Halalkalibacter": 0.20,
    "Aquibacillus": 0.25, "Litchfieldia": 0.25,
    "Ornithinibacillus": 0.25, "Domibacillus": 0.25,
    "Neobacillus": 0.25, "Ammoniphilus": 0.25,
    # Pseudomonads (often alive in arid)
    "Pseudomonas": 0.30, "Acinetobacter": 0.35,
    # UV/desiccation-resistant
    "Deinococcus": 0.15, "Truepera": 0.18,
    # Desiccation-tolerant cyanobacteria
    "Chroococcidiopsis": 0.20, "Acaryochloris": 0.25,
    # Thermo/halo-tolerant actinos
    "Rubrobacter": 0.30, "Conexibacter": 0.35,
    # Bacteroidota saline/halotolerant
    "Cytophaga": 0.40, "Pontibacter": 0.35,
    "Arcticibacter": 0.40, "Nibribacter": 0.35,
}

# Tier 5: Phylum-level priors (default for taxa not in T1/T2)
PHYLUM_PRIORS = {
    "Bacillota": 0.55,            # weak alive bias  (Bacilli + Clostridia mix)
    "Bacteroidota": 0.55,         # weak alive (halotolerant common)
    "Pseudomonadota": 0.60,       # neutral
    "Halobacterota": 0.30,        # alive in saline
    "Deinococcota": 0.30,         # alive (extremophiles)
    "Cyanobacteriota": 0.65,      # mixed - alive in biocrust, relic at depth
    "Actinomycetota": 0.65,       # mixed
    "Acidobacteriota": 0.75,     # weak relic (oligotroph background)
    "Gemmatimonadota": 0.75,     # weak relic
    "Chloroflexota": 0.75,        # weak relic
    "Planctomycetota": 0.75,      # weak relic (more typical of mesic)
    "Verrucomicrobiota": 0.65,
    "Myxococcota": 0.70,
    "Methylomirabilota": 0.85,    # methanotroph, no methane source -> relic
    "Spirochaetota": 0.80,
    "Fibrobacterota": 0.85,
    "Chlamydiota": 0.85,
    "Nitrospirota": 0.65,
    "Nitrospinota": 0.70,
    "Thermodesulfobacteriota": 0.70,
    "Desulfobacterota": 0.65,
    "Desulfobacterota_D": 0.50,   # NEUTRAL - CSP1-2 home in GTDB
    "Dadabacteria": 0.50,         # NEUTRAL - CSP1-2 home in SILVA
    "Patescibacteria": 0.75,      # CPR, often dormant
    "Crenarchaeota": 0.70,
    "Thermoplasmatota": 0.70,
    "Chloroflexi": 0.75,
}

# Default if no taxonomy info
DEFAULT_PRIOR = 0.65   # weakly relic (matches background biome composition)

# PMA correction: cell-biology adjustments
# Spore-formers (Bacillota, certain Clostridia, Actinomycetota with spores):
#   PMA may over-classify spores as alive. Apply DOWN-weight when model says
#   alive. We don't apply to model RELIC predictions.
SPORE_FORMER_GENERA = {
    "Bacillus", "Clostridium", "Geobacillus", "Aureibacillus", "Domibacillus",
    "Lysinibacillus", "Litchfieldia", "Halobacillus", "Halothermothrix",
    "Sediminibacillus", "Aquibacillus", "Streptomyces", "Salinispora",
    "Sporosarcina", "Sporomusa", "Anoxybacillus",
}


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def get_prior(genus, phylum):
    """Return P(relic | genus, phylum) using priors hierarchy."""
    if isinstance(genus, str):
        if genus in T1_RELIC_PRIORS: return T1_RELIC_PRIORS[genus]
        if genus in T2_ALIVE_PRIORS: return T2_ALIVE_PRIORS[genus]
        # Substring matches (e.g., "Halomonas sp." contains "Halomonas")
        for g, v in T1_RELIC_PRIORS.items():
            if g in genus: return v
        for g, v in T2_ALIVE_PRIORS.items():
            if g in genus: return v
    if isinstance(phylum, str) and phylum in PHYLUM_PRIORS:
        return PHYLUM_PRIORS[phylum]
    return DEFAULT_PRIOR


def is_spore_former(genus):
    if not isinstance(genus, str): return False
    if genus in SPORE_FORMER_GENERA: return True
    for g in SPORE_FORMER_GENERA:
        if g in genus: return True
    return False


def main():
    print("Loading inputs ...", flush=True)
    ind = pd.read_csv(CACHE / "test6_disconfirmation" /
                       "relic_indicator_with_damage_per_asv.tsv", sep="\t")
    tax = pd.read_parquet(CACHE / "taxonomy.parquet").reset_index().rename(
        columns={"ASV": "asv_id"})

    df = ind.merge(tax[["asv_id", "phylum", "genus"]], on="asv_id",
                     how="left")
    print(f"  ASVs: {len(df)}", flush=True)
    print(f"  with genus: "
          f"{df['genus'].notna().sum()}", flush=True)
    print(f"  with phylum: {df['phylum'].notna().sum()}", flush=True)

    # Compute prior per ASV
    df["prior_relic"] = df.apply(lambda r: get_prior(r["genus"], r["phylum"]),
                                     axis=1)
    df["spore_former"] = df["genus"].apply(is_spore_former)

    # Bayesian combination
    # logit(post) = logit(prior) + (logit(model) - logit(0.5))
    # This subtracts the "neutral" contribution of 0.5 from model so that
    # the model effectively contributes its evidence above/below neutral.
    df["log_prior"] = logit(df["prior_relic"])
    df["log_model"] = logit(df["relic_score_full_gb"])
    df["log_post"] = df["log_prior"] + (df["log_model"] - logit(0.5))

    # Spore-former correction: down-weight model "alive" predictions for
    # spore-formers (PMA over-counts spores). I.e., add +0.5 logit to push
    # toward relic.
    spore_correction = 0.5
    df.loc[df["spore_former"] &
              (df["relic_score_full_gb"] < 0.5), "log_post"] += spore_correction

    df["relic_score_posterior"] = sigmoid(df["log_post"])

    print("\n=== Score distribution comparison ===", flush=True)
    for q in (10, 25, 50, 75, 90):
        print(f"  p{q:>2}  model_only: "
              f"{np.percentile(df['relic_score_full_gb'], q):.3f}    "
              f"prior_only: {np.percentile(df['prior_relic'], q):.3f}    "
              f"posterior:  "
              f"{np.percentile(df['relic_score_posterior'], q):.3f}",
              flush=True)

    print("\n=== Pool counts ===", flush=True)
    for label, score_col in (("model_only", "relic_score_full_gb"),
                                  ("prior_only", "prior_relic"),
                                  ("posterior", "relic_score_posterior")):
        n_alive = (df[score_col] <= 0.3).sum()
        n_amb = ((df[score_col] > 0.3) & (df[score_col] < 0.7)).sum()
        n_relic = (df[score_col] >= 0.7).sum()
        print(f"  {label:<12}  alive={n_alive:>5}  ambig={n_amb:>5}  "
              f"relic={n_relic:>6}", flush=True)

    # Save
    df.to_csv(OUT / "relic_score_with_priors.tsv", sep="\t", index=False)

    # CSP1-2 specifically
    print("\n=== CSP1-2 status with priors ===", flush=True)
    csp_fasta = CACHE / "csp1-2_asvs.fasta"
    csp_ids = set()
    with open(csp_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                csp_ids.add(line[1:].strip().split()[0])
    csp = df[df["asv_id"].isin(csp_ids)].copy()
    print(f"  CSP1-2 ASVs in scoring: {len(csp)}", flush=True)
    print(csp[["asv_id", "phylum", "genus", "prior_relic",
                  "relic_score_full_gb", "relic_score_posterior"]].round(3)
          .head(24).to_string(index=False), flush=True)
    print(f"\n  CSP1-2 model-only  median score: "
          f"{csp['relic_score_full_gb'].median():.3f}", flush=True)
    print(f"  CSP1-2 prior-only  median score: "
          f"{csp['prior_relic'].median():.3f}", flush=True)
    print(f"  CSP1-2 posterior   median score: "
          f"{csp['relic_score_posterior'].median():.3f}", flush=True)
    csp.to_csv(OUT / "csp_with_priors.tsv", sep="\t", index=False)

    # Compare alive subsets
    print("\n=== Top 10 genera in alive subset (model_only) ===", flush=True)
    alive_m = df[df["relic_score_full_gb"] <= 0.3]
    print(alive_m["genus"].value_counts().head(10).to_string())

    print("\n=== Top 10 genera in alive subset (posterior) ===", flush=True)
    alive_p = df[df["relic_score_posterior"] <= 0.3]
    print(alive_p["genus"].value_counts().head(10).to_string())

    # How many ASVs change pool with priors?
    df["pool_model"] = np.where(df["relic_score_full_gb"] <= 0.3, "alive",
                                    np.where(df["relic_score_full_gb"] >= 0.7,
                                              "relic", "ambig"))
    df["pool_post"] = np.where(df["relic_score_posterior"] <= 0.3, "alive",
                                    np.where(df["relic_score_posterior"] >= 0.7,
                                              "relic", "ambig"))
    print("\n=== Pool transitions (model -> posterior) ===", flush=True)
    print(pd.crosstab(df["pool_model"], df["pool_post"], margins=True)
          .to_string())

    # Save
    summary = {
        "n_total": len(df),
        "n_with_T1_T2_prior": int((df["genus"].isin(set(T1_RELIC_PRIORS)) |
                                       df["genus"].isin(set(T2_ALIVE_PRIORS))
                                       ).sum()),
        "n_spore_former_corrected": int(df["spore_former"].sum()),
        "n_alive_model": int((df["relic_score_full_gb"] <= 0.3).sum()),
        "n_alive_posterior": int((df["relic_score_posterior"] <= 0.3).sum()),
        "csp_model_median": float(csp["relic_score_full_gb"].median()),
        "csp_posterior_median": float(csp["relic_score_posterior"].median()),
    }
    pd.DataFrame([summary]).to_csv(OUT / "comparison_with_no_prior.tsv",
                                        sep="\t", index=False)
    print(f"\nWrote outputs to {OUT}/", flush=True)


if __name__ == "__main__":
    main()
