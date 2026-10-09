# Rub' al-Khali bacterial biogeography

This public BORG repository is the reproducible companion to the manuscript
*Landscape-scale bacterial biogeography across the Rub' al-Khali*. It contains the active paper
and supplement, the analysis programs used for their claims, the canonical
machine-readable results, the submitted figures, regression tests, and
checksums and custody records for key methodological sources.

The manuscript reports a bacterial survey across the Rub' al-Khali.
Its main results concern geographic organization, paired soil-position
differences, environmental associations, predicted functional profiles, a
short observational association between rain and richness, relic-DNA checks,
and assay-aware low-biomass controls. Regression tests cover specified
calculations, figure semantics and selected result/prose checks; they do not
regenerate every number or the upstream raw-read processing.

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
only `main.tex` and `supplement.tex`; the main
manuscript contains no included prose fragments. It fixes
`SOURCE_DATE_EPOCH=1785888000` and `FORCE_SOURCE_DATE=1`, so repeated builds in
the pinned TeX environment produce byte-identical PDFs.

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
