# 03. PERMANOVA — habitat partitioning

**Question.** Which factors (compartment, trip, site) explain the most community variation?

**Method.** PERMANOVA on Bray–Curtis and Aitchison distances, with terms compartment, trip, site, and their interactions. `cache/permanova_by_term.tsv`, `cache/permanova_habitat.tsv`, `cache/distance_bray.parquet`, `cache/distance_aitchison.parquet`.

**Key results.**
- **Compartment** is the dominant axis (R² ≈ 0.10–0.15 depending on metric).
- **Site** (geography) is the next largest.
- **Trip** is smaller but non-negligible — supports treating trips as repeated measurements.
- LMM (`cache/lmm_compartment_trip.txt`, `cache/lmm_interaction.txt`) confirms compartment×trip interaction.

**Status.** solid (referenced in main text PERMANOVA table).

**Cross-refs.** 04 (distance-decay), 06 (iCAMP), 31 (alive temporal stability).
