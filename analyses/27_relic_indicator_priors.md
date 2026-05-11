# 27. Relic indicator with Bayesian + MAG-presence priors

**Question.** Can we further improve the relic indicator by adding fundamental biological priors (background knowledge about what is biologically alive vs dead)?

**Method.**
- **Bayesian taxonomic priors** (`scripts/relic_indicator_with_priors.py`):
  - T1 habitat exclusions (e.g. aquatic-only genera can't be alive in hyperarid soil).
  - T2 extremophiles (radiation-resistant, halotolerant → upweight as alive).
  - T5 phylum baselines.
  - Bayesian augmentation: log_odds(post) = log_odds(prior) + log_odds(evidence).
- **MAG-presence prior** (`scripts/relic_indicator_with_mag_prior.py`):
  - Match each EQ ASV against 16S extracted from 9,229 metagenome-assembled genomes (barrnap, then vsearch ≥99% identity).
  - Strong match (Δ = −2.5 logit shift toward alive), weak (Δ = −1.5), no match (Δ = 0).

**Inputs.**
- Track-C composite score (#25)
- MAG 16S library (`cache/mag_16s/`)

**Key results.**
- MAG-presence prior **overturns the CSP1-2 "relic" verdict** (median 0.744 → 0.241 — i.e. CSP1-2 has many MAG matches, so it's alive).
- 16/24 CSP1-2 ASVs match MAGs at ≥99% — clear alive evidence.
- Final relic-likelihood score: `cache/relic_priors/relic_score_with_mag_prior.tsv` — canonical per-ASV indicator.

**Status.** solid — canonical model. Used to define `alive` (score ≤ 0.3) and `relic` subsets for all downstream re-analyses (#29–33).

**Cross-refs.** 25, 28, 29–33, 32 (specific CSP1-2 correction).
