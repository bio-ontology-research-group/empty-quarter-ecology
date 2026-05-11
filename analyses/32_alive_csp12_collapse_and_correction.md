# 32. CSP1-2 keystone — collapse, then correction

**Question.** Does CSP1-2 retain its keystone status when we restrict to alive ASVs?

**Method.**
1. Build co-occurrence network on alive subset (compositional Spearman + BH q < 0.05). `scripts/relic_keystone_hunt.py`, `cache/keystone_hunt/`.
2. Re-rank keystones in alive networks per compartment.
3. Apply MAG-presence prior correction (`scripts/relic_indicator_with_mag_prior.py`) and re-classify CSP1-2 ASVs.

**Initial result.** CSP1-2 **drops out** of top keystones in all three compartments on the naïve alive subset.

**Correction.** MAG-presence prior (#27) **rescues** CSP1-2:
- 16 of 24 CSP1-2 ASVs match MAGs at ≥99% identity → CSP1-2 IS alive.
- After MAG-prior re-classification, CSP1-2 ranks **moderate (degree ≈ 3), not dominant (degree 22–41)** — i.e. it's a real alive node, but **not** the topological hub the original paper claimed.

**Interpretation.** Original "CSP1-2 keystone" claim was **two errors compounded**:
1. Inflated by relic-DNA inclusion.
2. Overstated by topological vs functional confusion (see #38, the discriminating-knockout test).

True status: CSP1-2 is one node within a **DOM-cycling guild** (Strategy A; #39), but not the single keystone.

**Status.** solid. Paper-rewriting finding.

**Cross-refs.** 27, 34, 38, 39, 53.
