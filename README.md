# Rub' al-Khali bacterial biogeography

This private BORG repository is the reproducible companion to the manuscript
*Landscape-scale bacterial biogeography across the Rub' al-Khali reveals
recurring spatial and soil-position structure*. It contains the active paper
and supplement, the analysis programs used for their claims, the canonical
machine-readable results, the four submitted figures, regression tests, and
byte-verifiable copies of key methodological sources.

The manuscript reports the first broad bacterial survey of the Rub' al-Khali.
Its main results concern geographic organization, paired soil-position
differences, environmental associations, predicted functional profiles, a
short observational association between rain and richness, relic-DNA checks,
and assay-aware low-biomass controls. The wording tests connect each numerical
claim to the corresponding result file and enforce the stated limits on
interpretation.

## Repository relationship

Large biological inputs, sample and control metadata, climate products,
geochemistry, pH measurements, and the knowledge graph live in the separate
private repository
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
make figures
make paper
```

`make figures` renders the four figures from the committed canonical result
tables in a temporary directory and requires byte-identical PDFs and manifest
entries. `make paper` builds only `main.tex` and `supplement.tex`; the main
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

This is a private submission candidate, not a public data deposit. Public
accessions, final licences, and a DOI remain author-controlled release gates.
