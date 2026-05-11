# 18. Within-site Hill profile fits

**Question.** What is the Hill diversity profile (q = 0, 1, 2) per site, and does it fit a power-law form?

**Method.** Fit Hill profile per site, per compartment. `scripts/run_within_site_hill.py`, output `cache/within_site_hill_fit.tsv`.

**Key results.**
- EQ within-site Hill profile decays steeply (typical of dominance-skewed systems).
- Best-fit q-exponent **~2.4** for EQ.
- Sets up the cross-cohort comparison (#49): Atacama gives ~1.16 (much flatter, more evenly diverse), so EQ's pattern doesn't transfer.

**Status.** solid.

**Cross-refs.** 49.
