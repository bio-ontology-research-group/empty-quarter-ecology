# 50. betA guild — metagenomic census

**Question.** Is the predicted "betaine-producing keystone guild" (from PICRUSt2) confirmed at the metagenome level?

**Method.** Census **9,229 MAG bins** from **296 metagenomes** (EQ + ancillary). Search for K00108 (betA, choline dehydrogenase) using HMM with strict thresholds; per-bin and per-sample summaries. `scripts/run_betA_guild_analysis.py`, `scripts/run_osmoprotectant_transporters.py`.

**Inputs.**
- 9,229 MAG bins (rhizosphere + deep + surface)
- 296 metagenomes

**Key results.**
- **152 producer bin IDs** carry K00108 with confident HMM scores.
- Producers are **6× enriched in rhizosphere** vs deep (rhizosphere ~6%, deep ~1%).
- Per-sample summary: producer count varies 10× across (site, comp); maximum producer fraction ~12% in stressed rhizosphere.

**Outputs.**
- `cache/betA_guild_census.tsv`
- `cache/betA_producers_per_bin.tsv`
- `cache/betA_producers_locator.tsv`
- `cache/betA_per_sample_summary.tsv`
- `cache/betA_field_coverage.tsv`
- `cache/betaine_uptake_census.tsv`
- Per-genus betA HMM tables: `cache/Flavisolibacter_K00108.tbl`, `cache/Rubellimicrobium_K00108.tbl`, `cache/Solirubrobacter_K00108.tbl`, `cache/Telluribacter_K00108.tbl`.

**Interpretation.** The producer-guild exists at the metagenomic level — but **#51 immediately complicates the story** (leak-asymmetry fails). And **#52** further shows multiple genera (not just CSP1-2 / Nibribacter) host betA.

**Status.** solid census; mechanism reinterpreted by #51, #52.

**Cross-refs.** 13, 51, 52.
