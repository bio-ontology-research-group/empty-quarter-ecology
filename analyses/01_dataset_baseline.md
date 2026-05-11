# 01. Dataset baseline

**Question.** What is the sampling design and pipeline that all downstream analyses build on?

**Design.** 1,237 16S rRNA amplicon samples from 60 sites spanning the Rub' al-Khali (Empty Quarter, Saudi Arabia), 5 seasonal campaigns (Trips 1–5, 2023–2024), 3 compartments per site (surface 0 cm, deep 3–5 cm, rhizosphere).

**Trips & sample-name prefixes**
- T1 (Mar 2023, n≈60, no prefix)
- T2 (Jul 2023, n=8 sites, prefix `T`)
- T3 (Feb 2024, n=60, prefix `F`)
- T4 (Aug 2024, n=60, prefix `S`)
- T5 (Oct 2025, n=60, prefix `V`)

**Pipeline.** nf-core/Ampliseq v2.14.0 → QIIME 2 v2024.10 → PICRUSt2 v2.4.1 for functional prediction. Taxonomy via SILVA. XRF geochemistry was run on 158 samples from T5 plus 1 rep per (site, comp) for T1–T4 = 725 measurements.

**Key inputs.**
- `cache/feature_table.parquet` — ASV × sample counts (canonical)
- `cache/taxonomy.parquet` — ASV-level SILVA assignments
- `cache/metadata.parquet` — sample-level metadata
- `data/geodata/trip*_geodata.tsv` — site coords + CenterDate per trip
- `cache/qc_report.json`

**Helper.** `scripts/_sample_parse.py` parses sample names with trip prefixes correctly — use it everywhere.

**Status.** solid (canonical).

**Cross-refs.** Every other entry consumes these tables.
