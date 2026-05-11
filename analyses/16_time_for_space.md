# 16. Time-for-space substitution

**Question.** Can we use trip-to-trip variation at the same site as a proxy for cross-site climate variation?

**Method.** Build per-(site, trip) climate vector; compute within-site temporal vs across-site spatial dissimilarities. `scripts/run_time_for_space.py`, output `cache/tfs/`.

**Key results.**
- Within-site temporal climate variation spans a meaningful fraction of across-site climate space.
- Supports using transition-based modeling (#41, #44, #45) rather than relying solely on cross-sectional comparisons.

**Interpretation.** This is the methodological foundation for the longitudinal climate projection (v2, v3) — temporal trip-to-trip variation is informative about climate-driven shifts.

**Status.** solid.

**Cross-refs.** 41, 44, 45.
