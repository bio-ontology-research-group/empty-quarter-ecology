# 39. Two-strategy temporal analysis

**Question.** Does the alive EQ community partition into two coherent, anti-correlated strategies — and how does each strategy track over time?

**Method.** Define:
- **Strategy A** (Bacteroidota DOM-cyclers + Massilia): Nibribacter, Flavisolibacter, Solirubrobacter, Telluribacter, Rubellimicrobium, Massilia.
- **Strategy B** (halotolerant): Aquibacillus, Oceanobacillus, Halobacillus, Halomonas, Pseudomonas, …

Per (sample, trip), compute log₂(ΣA / ΣB); dominant = A if log₂ > 0 else B. Aggregate per-trip per-site. `scripts/two_strategy_temporal.py`, `cache/two_strategy_temporal/`.

**Key results.**
- **50 sites are A-dominant**, **10 sites B-dominant** in baseline (T1 + T3 combined).
- **51 (site, comp) cells switch dominance across trips.**
- Strategy A correlates with d7 precip (ρ = +0.21).
- **Strategy A declines T1 → T5** (more wet → more A; system dries on average across the study window).
- Surface compartment has the highest A-dominance fraction (75–92%).

**Interpretation.** Two anti-correlated strategies with a precipitation-pulse trigger. Replaces the obscured single-keystone story (#32) with a clean **alternate stable states** narrative.

**Status.** solid (with scrutiny round #40 confirming effect size is moderate).

**Outputs.**
- `cache/two_strategy_temporal/per_sample_strategy_with_precip.tsv`

**Cross-refs.** 27, 38, 40, 41, 42, 45.
