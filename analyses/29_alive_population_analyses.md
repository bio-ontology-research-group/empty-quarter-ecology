# 29. Alive-population analyses (post-relic-filter)

**Question.** What does the EQ community look like once relic ASVs are excluded?

**Method.** Apply canonical relic score (#27, `relic_score_with_mag_prior.tsv`) to partition ASVs at threshold 0.3 (alive ≤ 0.3 vs relic > 0.3). Produce alive-only feature table; recompute diversity, abundance, taxonomic composition. `scripts/relic_population_analyses.py`, output `cache/relic_population/`.

**Inputs.**
- `cache/feature_table.parquet`
- `cache/relic_priors/relic_score_with_mag_prior.tsv`

**Outputs.**
- `cache/feature_table_alive.parquet`
- `cache/feature_table_relic.parquet`

**Key results.**
- **~84% of ASVs** classified as relic; **60–80% of reads** are relic.
- Alive community dominated by:
  - **Halotolerant Bacilli** (Aquibacillus, Oceanobacillus, Halobacillus, …)
  - **Halomonas**
  - **Pseudomonas**
- Cosmopolitan ASVs are **more alive**, not more relic (counter to one prior hypothesis).
- CSP1-2 classifies as relic by 16S alone, but rescued by MAG prior (#32).

**Interpretation.** EQ "is mostly dead." The alive fraction is the well-known **halotolerant Bacilli + γ-proteobacteria** crew. The rare diversity is in the relic / wind-deposited pool.

**Status.** solid (with #28 calibration caveat for rare taxa).

**Cross-refs.** 27, 28, 32, every alive-only re-analysis (30–33).
