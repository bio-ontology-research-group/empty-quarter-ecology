# 35. Nibribacter MAG retrieval + functional annotation

**Question.** What does Nibribacter actually do? Pull its MAGs, annotate them.

**Method.**
- Identify all MAGs assigned to Nibribacter genus in the metagenomic guild census.
- ORF prediction with prodigal.
- Kofam HMM search (no kofamscan binary available; used plain hmmsearch + ko_list.txt thresholds).
- `scripts/nibribacter_kegg_analysis.py` (naïve, regex-based — superseded by #36).

**Outputs.**
- `cache/nibribacter_mags/per_mag_ko_assignments.tsv`
- Raw HMM tables `cache/nibribacter_mags/*.hmm.tbl`
- 15PRr3_SemiBin_102.faa (representative ORFs)

**Initial findings (naïve).** Lots of betaine-related hits, polysaccharide degradation, glycoside hydrolases, sigma factors, …

**Critical bug.** The regex `betA` matched **anything containing "beta-"** in KO definitions — including K01918 (pantoate-β-alanine ligase). Reported "96 betaine ORFs" was inflated; real count is **1**. Fixed in #36.

**Status.** superseded by #36 — preserved here only as historical record.

**Cross-refs.** 36, 50.
