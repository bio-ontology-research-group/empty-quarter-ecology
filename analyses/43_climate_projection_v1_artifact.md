# 43. Climate projection v1 — cross-sectional (ARTIFACT)

**Question.** Under CMIP6 SSP scenarios, what fraction of EQ sites flip from A-dominant to B-dominant?

**Method (v1).** Per-sample logistic regression: `dominant ~ AnnualMeanTemp + AnnualTotalPrecip + sabkha_score + Latitude + Longitude`. Apply CMIP6 ΔT, ΔP_pct deltas (SSP1-2.6, SSP2-4.5, SSP3-7.0 at 2050 + 2100). `scripts/two_strategy_climate_projection.py`, `cache/two_strategy_projection/`.

**Original headline.** Under SSP3-7.0_2100: fraction-B collapses **28.9% → 1.9%** — i.e. warming dramatically reduces B.

**Why it's an artifact.** The cross-sectional logit captures **geographic clustering** of sabkha sites (which are coastal, cooler, AND happen to be B-dominant). AnnualMeanTemp coefficient = **−0.91** — but this is identifying the sabkha-coastal cluster, NOT a causal "warming kills B" relationship.

When you then apply +4 °C to every sample, the model interprets this as "everywhere now looks less like the cold-coastal sabkha cluster" → predicts collapse of B. **This is exactly the cross-sectional confounding trap.**

**Lesson.** Cross-sectional regression of (state) on (climate + chemistry) cannot disentangle climate-driven shifts from geographic-cluster spurious associations. Must use **longitudinal trip-to-trip transitions** (#44–45).

**Status.** **retracted (artifact).** Preserved for documentation; do not cite.

**Outputs (for record).**
- `cache/two_strategy_projection/scenario_summary.tsv`
- `cache/two_strategy_projection/per_site_projection.tsv`
- `cache/two_strategy_projection/per_site_scenarios.tsv`

**Cross-refs.** 44, 45.
