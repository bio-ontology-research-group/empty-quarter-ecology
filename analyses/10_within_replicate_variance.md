# 10. Within-replicate variance vs between-site (Test 3)

**Question.** Are the field replicates meaningfully tight relative to between-site variance? (Quality check on whether site is a real factor.)

**Method.** Per-(site, comp, trip) within-replicate Bray–Curtis variance vs between-site BC variance at fixed trip × compartment. `scripts/test3_within_replicate_variance.py`, output `cache/test3_variance/`.

**Key results.**
- Within-replicate variance is consistently **lower than between-site variance** (~1.5–3× lower).
- Replicates cluster tightly; site is a meaningful factor.

**Interpretation.** Sampling design is sound — site-level effects in downstream analyses are not noise.

**Status.** solid.

**Cross-refs.** Background QC for everything.
