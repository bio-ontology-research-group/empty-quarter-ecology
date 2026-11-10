# Rub’ al-Khali Ecology Manuscript

## Manuscript Sources & Build

The active manuscript draft is `alt_stats.tex`, with companion supplementary materials in `supp_altstat.tex`.
(The older `main.tex` and `supplement.tex` are preserved for historical record).

To compile:
```bash
pdflatex alt_stats.tex && bibtex alt_stats && pdflatex alt_stats.tex && pdflatex alt_stats.tex
pdflatex supp_altstat.tex && pdflatex supp_altstat.tex
```

## Scientific Narrative & Data Status

- **Focal Findings**: Spatial lineage replacement across a 1,043 km west--east hyperarid transect; surface enrichment of compatible solute uptake (BetL) and oxidative defence (AhpCF); selective depletion of high-affinity group 1l hydrogenase genes beside roots.
- **Amplicon Survey**: 1,186 deduplicated 16S rRNA gene profiles across 60 sites and five sampling campaigns.
- **Metagenome Cohort**: 554 paired-end shotgun sequencing runs from 551 soil samples (with three duplicate-run samples), yielding 549 MEGAHIT assemblies.
- **MAG Catalogue**: 5,438 species representatives (4,847 bacterial, 591 archaeal) dereplicated by dRep at 95% ANI, classified with GTDB-Tk v2.7.2 against GTDB release 232 (`gtdb_r232_classification`).
- **Gene Calling**: Pyrodigal v3.7.1 in metagenomic mode (`EQ-PROT-PYRODIGAL-REP-v1.0.0`), distinct from the companion Bakta v1.12.0 release annotation product (`EQ-PROT-BAKTA-v1.0.0`).
- **Matched Cohort**: 442 sample profiles across 58 sites matched between the 1,186-sample amplicon cohort and CoverM-mapped metagenomic abundances.
- **Pending Reconciliations**:
  - `[PENDING: ECOL-READ-CUTOFF-RECONCILIATION]`: Tracking the recent manuscript text citing a $<10{,}000$ read cutoff versus canonical data descriptor and companion records showing ten profiles with $<1{,}000$ reads (including `T1Dr1` at 934 reads).
  - `[PENDING: TAXONOMY-SNAPSHOT-SENSITIVITY]`: Tracking order- and genus-adjusted `gRodon2` regressions under canonical GTDB release 232 versus analytical snapshot.
  - Full microTrait ecology integration paragraph awaits root statistics.
