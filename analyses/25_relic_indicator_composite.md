# 25. Relic-likelihood indicator (Path B composite)

**Question.** Without PMA on every sample, can we infer per-ASV relic-likelihood from amplicon + auxiliary features?

**Method (Path B).** Train a composite classifier (logistic regression + gradient boosting) on the 18 PMA pairs as ground truth, with features per ASV:
1. ASV sequence features (length, GC content, dinucleotide ratios) — for **Track C amplicon damage proxies**.
2. Detection patterns across compartments / trips.
3. Per-sample abundance distribution.

Scripts: `scripts/relic_indicator.py` (initial), `scripts/relic_indicator_with_damage_proxies.py` (Track C augmented).

**Inputs.**
- PMA T/UT pairs (`/home/leechuck/Public/software/empty-quarter/relic-dna/`)
- `cache/feature_table.parquet`
- `cache/taxonomy.parquet`

**Key results.**
- Initial logistic model: **AUC = 0.703** (cross-validated).
- + Track C damage proxies: **AUC = 0.785** (substantial gain).
- ASV length, GC content, and a TC/CC ratio (UV-damage signature analog) all moved in the expected direction.

**Outputs.**
- `cache/relic_priors/relic_score.tsv` (initial)
- Track-C augmented score in same dir

**Status.** solid foundation; superseded by **prior-augmented** indicator #27.

**Cross-refs.** 24 (PMA ground truth), 26 (mapDamage), 27 (Bayesian + MAG prior), 28 (sensitivity).
