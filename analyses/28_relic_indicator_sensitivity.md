# 28. Relic indicator — critical sensitivity tests

**Question.** How robust is the relic-likelihood model? Where does it fail?

**Method.** A battery of stress tests on the prior-augmented indicator (#27). `scripts/relic_indicator_sensitivity.py`, outputs in `cache/relic_sensitivity/`.

1. **Cross-site validation** — train on N−1 sites, test on held-out site.
2. **Random null** — permute labels, refit; compute null AUC distribution.
3. **Abundance stratified** — does indicator work for rare ASVs as well as abundant?
4. **Per-prior contribution** — what does each Bayesian prior add?

**Key results.**
- Cross-site **AUC drops 0.79 → 0.69** (substantial generalization gap, but still significantly above null).
- Random null AUC distribution centered at 0.50 (clean).
- **False-relic rate ≈ 80% for low-abundance specialists** — model loses calibration in the long tail.
- Each prior layer adds incrementally; MAG-presence prior is the largest single contribution.

**Interpretation.** Indicator is **trustworthy for moderately-to-highly abundant taxa** but **unreliable for rare specialists**. All downstream alive-only re-analyses (#29–33) are dominated by the well-calibrated abundant fraction, so this is OK for the strategy A/B story — but rare-taxon-specific claims would be unsafe.

**Status.** solid. Bounds the scope of confident inference.

**Cross-refs.** 25, 26, 27, 29–33.
