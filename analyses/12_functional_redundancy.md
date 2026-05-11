# 12. Functional redundancy — Allison–Martiny slope (Test 4)

**Question.** How redundant is the EQ microbiome — i.e. how much taxonomic turnover is needed for a unit of functional change?

**Method.** Allison & Martiny (2008) framework: regress functional dissimilarity (PICRUSt2 pathway) against taxonomic dissimilarity (ASV) across all site-pairs. Slope < 1 = redundant. `scripts/test4_allison_martiny.py`, output `cache/test4_allison_martiny/`.

**Key results.**
- **Slope = 0.21** — strongly sub-linear → **high functional redundancy**.
- Note: alive-only re-analysis (#33) initially showed a doubling of this slope, but that was later confirmed as an artifact of the relic filter (see relic indicator corrections #28).

**Interpretation.** Heavy redundancy is consistent with depth-as-refugium framing and with the dispersal-driven story (any of many taxa can do the same job).

**Status.** solid (with note that alive-only doubling was an artifact).

**Cross-refs.** 4, 13, 14, 33.
