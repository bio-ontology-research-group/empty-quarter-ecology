# 24. PMA viability QC

**Question.** Are the PMA-treated / untreated paired samples consistent enough to use as ground truth for a relic-DNA detection model?

**Method.** For 18 paired (T = PMA-treated, UT = untreated) samples from Trip 5, compute per-pair correlation, ICC, and presence/absence consistency. `scripts/pma_viability_qc.py`.

**Data location.** Direct experiment files at `/home/leechuck/Public/software/empty-quarter/relic-dna/` (Trip 5 paired T/UT samples).

**Key results.**
- **Per-pair ICC = 0.53** — moderate.
- PMA treatment removes a substantial fraction (median ~70%) of reads, consistent with most DNA being extracellular/dead.
- Workable for ground truth, though not perfect — drives our **composite-indicator** approach (#25) rather than direct PMA classification.

**Status.** solid (workable but moderate-quality ground truth).

**Cross-refs.** 25, 26, 27, 28.
