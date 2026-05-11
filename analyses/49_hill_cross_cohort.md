# 49. Hill diversity — cross-cohort transferability test

**Question.** Does the EQ-specific Hill diversity profile (#18, q-exponent ~2.4) transfer to Atacama?

**Method.** Fit Hill profile per site in Atacama dataset under same protocol as EQ within-site fits. `scripts/run_hill_cross_cohort.py`, output `cache/hill_cross_cohort_atacama.tsv`, `cache/hill_cross_cohort_atacama_fits.txt`.

**Key results.**
- EQ within-site Hill q-exponent: **~2.4** (steep, dominance-skewed).
- Atacama Hill q-exponent: **~1.16** (much flatter, more uniformly diverse).
- Atacama profile is approximately **linear**, EQ is sharply non-linear.

**Interpretation.** The EQ Hill profile **does not transfer** — Atacama is a fundamentally less-dominance-skewed system. Important for pre-registration: claiming "deserts follow Hill q~2" based on EQ data alone would be a false-generalisation.

**Status.** solid. Constrains generality claims.

**Cross-refs.** 18, 47, 48.
