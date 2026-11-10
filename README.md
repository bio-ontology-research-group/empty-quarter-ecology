# Rub' al-Khali microbial biogeography and functional differentiation

This public BORG repository is the reproducible companion to the manuscript
*Microbial biogeography and functional differentiation across the Rub' al-Khali*
(historical amplicon-only snapshot title: *Landscape-scale bacterial biogeography across the Rub' al-Khali*).
It contains the active integrated paper (`empty-quarter-amplicon/alt_stats.tex`)
and supplement (`empty-quarter-amplicon/supp_altstat.tex`), the preserved historical
amplicon snapshot (`main.tex` and `supplement.tex`), the analysis programs used
for their claims, the canonical machine-readable results, the submitted figures,
regression tests, and checksums and custody records for key methodological sources.

The active manuscript integrates the 16S rRNA survey (1,186 deduplicated primary-transect
profiles across 60 sites) with genome-resolved shotgun metagenomics (546 core-transect
libraries across 59 sites; 545 read-depth-modeled libraries; 442 exact specimen-matched
amplicon--metagenome libraries across 58 sites; 441 matched read-depth-modeled libraries)
across the 5,438-representative genome catalogue (`EQ-MAG-REFINED-v1.0.0`). Its main
results concern landscape distance decay and eastern megabarchan taxonomic versus functional
diversity contrasts, compartment turnover and surface compatible-solute/oxidative-stress
enrichment (`BetL`, `AhpCF`, and the complete 41-rule microTrait HMM and KO panel),
root-adjacent trace-gas oxidation depletion (`group 1l` [NiFe]-hydrogenase), subsurface
Thermoplasmatota enrichment (bridged to the companion archaeal proteome study), environmental
and rainfall associations, relic-DNA checks, and assay-aware low-biomass controls.

The active manuscript preserves Rund Tawfiq's 2 October 2026 Overleaf revision
and her 9 October 2026 figure revision (`handover/figure-revision-20261009-rund.md`),
with factual corrections documented in the coauthor review package. The paper
and supplement in `empty-quarter-amplicon/` are synchronized with Overleaf,
including the included tables and generated figures. The analysis inventory
and validation records identify the current result generation for each claim.

`DATA_REPOSITORY.lock` pins the shared data release containing the Site 52
coordinate correction. The daily Open-Meteo correction in
`analysis/v3/open_meteo_site52_corrected_20261004/` replaces that release's
remaining stale Site 52 daily series; its original response, request and hashes
are packaged alongside the corrected input. The common-denominator pH
partition is in `analysis/v3/ph_partition_20261004/`.

## Repository relationship

Large biological inputs, sample and control metadata, climate products,
geochemistry, pH measurements, and the knowledge graph live in the separate
public repository
[`empty-quarter-data-paper`](https://github.com/bio-ontology-research-group/empty-quarter-data-paper).
[`DATA_REPOSITORY.lock`](DATA_REPOSITORY.lock) pins its exact commit and the
SHA-256 digests of its release, bulk-input, environment, and workflow
manifests. The two repositories therefore form one auditable release without
duplicating multi-gigabyte inputs in Git.

The companion archaeal proteome and mechanistic investigation lives in
`empty-quarter-proteome` (local checkout `/home/leechuck/Public/software/empty-quarter-proteome`).
The common genome cohort, GTDB release 232 taxonomy reconciliation, and the pending
status of mechanistic results are documented in
[`docs/SHARED_FOUNDATION.md`](docs/SHARED_FOUNDATION.md).

## Quick verification

Clone the data repository beside this repository, install its bulk inputs, and
then run:

```bash
make bootstrap DATA_REPO=../empty-quarter-data-paper
make verify
make test DATA_REPO=../empty-quarter-data-paper
make figures PYTHON=../empty-quarter-data-paper/.conda-env/bin/python
make paper
```

`make figures` renders the five main-text figures and three supplementary or
archived figures from their current result tables in a temporary directory and requires byte-identical PDFs and review
manifests. Byte-level rendering uses the exact Linux environment in the data
repository's `environment/conda-linux-64.lock`; the renderer fails before
writing output if Python, Matplotlib, or FreeType differs. `make paper` builds
the active manuscript (`supp_altstat.tex` followed by `alt_stats.tex`, which
inputs `mag_lifeonair_methods.tex`, while `supp_altstat.tex` inputs
`mag_lifeonair_supplement.tex`). The preserved historical amplicon-only targets
(`main.tex` and `supplement.tex`) can be built separately via `make paper-legacy`.
Both targets fix `SOURCE_DATE_EPOCH=1785888000` and `FORCE_SOURCE_DATE=1`, so
repeated builds in the pinned TeX environment produce deterministic PDFs.

## Genome-resolved integration and Stage 4 synthesis matrices

`analysis/genome_integration_20261108/` contains the Stage 4 synthesis code
(`code/stage4_ecology_synthesis.py`, `code/common_science.py`), verification tests
(`tests/test_science_20261108.py`), and machine-readable result matrices (`results/`):

- `task_d_h6_baseline_vs_protection_surface_deep.tsv`: 49 surface-minus-subsurface
  rows across stress-protection, compatible-solute, and baseline DNA-repair panels
  (paired site-block bootstrap intervals, site-level Wilcoxon signed-rank BH
  $q$-values across all 41 rules and within exploratory families, INLA spatial-model
  posterior summaries across 7 specifications, and three-component genus turnover
  decomposition under `EQ-TAX-LEGACY-ANALYTICAL-v0`).
- `task_d_h6_h3_systematic_matrix.tsv` and `task_d_all_41_microtrait_and_curated_matrix.tsv`:
  complete 250-row and 136-row panel-by-contrast matrices covering all 41 microTrait
  HMM rules, 39 available KO projections (`rep_ko_v2.tsv`, `EQ-KO-LEGACY-ANALYTICAL-v0`),
  and curated/legacy trace-gas and stress panels across `Surface-minus-Deep` and
  `Root-minus-Deep`. Machine-readable pathway interpretation flags in column
  `pathway_interpretation_flag` distinguish `HMM_KO_ROUTE_DISAGREEMENT_FLAGGED`
  (single-step HMM vs KO/two-step strict trehalose synthesis disagreement) and
  `SINGLE_MARKER_NOT_FULL_SPORULATION_PATHWAY` (`spoIIGA` protease marker alone).
- `task_d_eastern_joint_diversity_contrast.tsv`: conditional spatial-block bootstrap
  sensitivity for the eastern megabarchan versus uruq diversity slope difference
  ($\beta_{\mathrm{tax}} - \beta_{\mathrm{KO}}$; $685/1{,}000$ valid sample-level
  block resamples and $315/1{,}000$ rank-deficient draws omitting the eastern
  megabarchan group).
- `turnover_summary.tsv`, `spatial_fixed_effects.tsv`, `PICRUSt_MAG_effect_concordance.tsv`,
  and `task_d_turnover_allread_vs_mapped_comparison.tsv`: source tables for the
  three-component genus turnover decomposition, INLA spatial fits, and PICRUSt2 vs
  MAG KO concordance.

Every real knowledge-graph generation or semantic-validation run must execute
on `ws` or Ontolinator. The repository deliberately has no local KG target.
Use [`REPRODUCE.md`](REPRODUCE.md) and
[`workflow/run_on_remote.sh`](workflow/run_on_remote.sh) for the complete
cross-paper rebuild.

## Contents

- `empty-quarter-amplicon/`: active manuscript, supplement, bibliography,
  reviewed PDFs, generated pH constants, and submitted figures;
- `analysis/`: analysis programs and canonical result bundles;
- `tests/`: statistical, provenance, claim, and manuscript regression tests;
- `metadata/`: the orientation boundary needed for the landscape figure;
- `literature/CITATION_SOURCES.tsv`: DOI, retrieval, and local-custody
  checksums for methodological sources, without redistributing source bytes;
- `environment/`: the same pinned environment specification used by the data
  repository; and
- `archive/` and `ecology-paper/`: explicit guards for retired manuscript
  material.

This is a submission candidate, not a public data deposit. Public
accessions, final licences, and a DOI remain author-controlled release gates.
