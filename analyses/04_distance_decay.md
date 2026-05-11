# 04. Distance-decay — taxonomic vs functional

**Question.** Does taxonomic similarity decay with geographic distance faster than functional similarity?

**Method.** Pairwise BC-similarity vs km-distance for ASV-level (taxonomic) and PICRUSt2 pathway-level (functional) data, per compartment. `cache/distance_decay.tsv`, `cache/distance_decay_tax_vs_func.tsv`.

**Key results.**
- Taxonomic decay steeper than functional decay in all compartments → **functional convergence** despite taxonomic turnover.
- Decay slope steepest in rhizosphere, shallowest in deep — supports depth-as-refugium framing (RQ04 supplement).

**Interpretation.** Classic redundancy result: many taxa, similar functional repertoires. Bookends with the Allison-Martiny redundancy test (#12) and Black Queen ratio (#13).

**Status.** solid.

**Cross-refs.** 12, 13, 14, 33.
