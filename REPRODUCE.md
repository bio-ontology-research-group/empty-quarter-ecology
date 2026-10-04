# Reproducing the ecology manuscript

The ecology paper uses two levels of verification. The short path checks the
committed claim evidence, rerenders every submitted figure, and rebuilds both
PDFs. The workflow path starts from the checksum-pinned inputs in the data
repository and reruns the configured analysis and knowledge-graph processes.
It is not a complete raw-read-to-paper reconstruction. The September 2026
revision also uses explicit replays whose outputs are the dated evidence
directories under `analysis/v3/*_20260909`; each carries its own manifest and
checksums.

## 1. Obtain the exact repositories

```bash
git clone https://github.com/bio-ontology-research-group/empty-quarter-ecology-reproducibility.git
git clone https://github.com/bio-ontology-research-group/empty-quarter-data-paper.git
cd empty-quarter-data-paper
git checkout "$(awk -F '\t' '$1 == "commit" {print $2}' ../empty-quarter-ecology-reproducibility/DATA_REPOSITORY.lock)"
bash scripts/release/download_bulk_artifacts.sh
bash scripts/release/bootstrap_package_layout.sh .
cd ../empty-quarter-ecology-reproducibility
```

`scripts/release/bootstrap_data_dependency.sh` checks the data commit and all
seven pinned manifest and environment digests before creating relative
compatibility links. It
will not replace an existing path or accept a different data revision.

## 2. Recreate the environment

The data repository carries the authoritative exact Linux/x86-64 environment.
Create it before byte-level figure verification or the complete workflow:

```bash
cd ../empty-quarter-data-paper
make env-linux-exact
export PATH="$PWD/.conda-env/bin:$PATH"
cd ../empty-quarter-ecology-reproducibility
```

The explicit lock fixes every Conda package build, including Matplotlib
`3.9.4` and FreeType `2.14.3`; the small pip overlay is hash-locked and cannot
replace Conda dependencies. A lighter CPython 3.11 environment remains
available for numerical tests that do not render canonical PDFs:

```bash
uv venv --python 3.11 .venv
uv pip sync --python .venv/bin/python \
  ../empty-quarter-data-paper/environment/requirements.lock.txt
```

The editable Conda recipe additionally pins Java, Groovy, R, MAFFT, FastTree,
and the other programs used by the complete workflow. Raptor is built from its
checksum-pinned source archive. Every executed remote workflow records the
versions it actually found and the explicit-lock digest.

## 3. Verify claims, figures, and papers

```bash
make bootstrap DATA_REPO=../empty-quarter-data-paper
make verify PYTHON=.venv/bin/python
make test PYTHON=.venv/bin/python DATA_REPO=../empty-quarter-data-paper
make figures PYTHON=../empty-quarter-data-paper/.conda-env/bin/python
make paper
```

The test suite checks numerical claims against canonical TSV/JSON outputs,
uncertainty intervals and multiplicity decisions, control and PMA boundaries,
rainfall sensitivity, manuscript structure, the citation-custody ledger, cross-paper
author/title consistency, and deterministic result bundles. Figure rendering
uses only committed result files and fails if any regenerated PDF differs from
the reviewed manuscript copy.

These checks do not generate a knowledge graph. They are safe on a local
workstation.

Literature PDFs and other complete source snapshots are not distributed in
this repository. `literature/CITATION_SOURCES.tsv` records their DOI or stable
source, retrieval basis, and the SHA-256 of the copy inspected by the authors.
The source bytes remain in local custody subject to their licences.

## 4. Run the complete workflow remotely

Copy or clone both exact repositories on `ws` or Ontolinator, install the bulk
inputs, create the locked environment, and run:

```bash
cd empty-quarter-ecology-reproducibility
bash workflow/run_on_remote.sh ../empty-quarter-data-paper \
  ./results/remote-validation-$(date -u +%Y%m%dT%H%M%SZ)
```

The institutional profile uses `leechuck-office` (`ws`) or `cbontsr01`
(Ontolinator). External users can run the same wrapper on their own Linux
compute host by setting `EQ_EXECUTION_CONTEXT=external`; the full graph path
requires at least 32 GB Java heap and additional memory and scratch storage.
The wrapper activates the locked environment and verifies the repository lock.
It runs the currently wired core/advanced ecology and graph stages. Nextflow
writes a trace, report, timeline, DAG, source-state record, environment record,
commands, logs, and SHA-256 manifests into a new output directory.

The current full-workflow entry point still copies the frozen pH analysis and
does not regenerate every later taxon, landform, trait and control follow-up.
Those modules have source programs and result tables in `analysis/`. A complete
final manuscript replay requires their integration and a new frozen cross-paper
release after the September source corrections. The claim tests compare many
numbers to frozen outputs; a passing short test run is not evidence of fresh
regeneration from raw reads. Public raw-read accession and upstream prediction/
genome provenance requirements are listed in the companion data repository.

The following committed result directories were produced before the Site 52
coordinate correction and still carry the uncorrected position for that site
(Trip 1 and Trip 3 rows, or the site average over all trips). They are not
cited for coordinate-dependent numbers in the manuscript copy:
`analysis/v3/environment_associations/climate_site_summary.tsv` (the
corrected climate-diversity sensitivity is in
`analysis/v3/biology_context_corrected_20260909/`), the
`analysis/v3/rain_pulse_response*` and `analysis/v3/rain_pulse_sensitivities/`
cohorts (superseded by `analysis/v3/rain_calendar_refit_20260909/`),
`analysis/v3/xrf_community_rescue/xrf_alpha_analysis_table.tsv`,
`analysis/batch-adjacency-2026-09-07/batch_meta.tsv`, and the control-filter
sensitivity under `analysis/rerun-controls-2026-08-30/sensitivity/`, whose
distance-decay comparison was run on the uncorrected coordinates.

Do not use `-resume` after changing source code, input data, manuscript text, or
the data lock. A local Nextflow stub is useful only for checking wiring and is
not scientific or semantic validation.

## 5. Reproducibility boundary

The release can reconstruct the submitted analyses, figures, manuscript, and
derived knowledge graph from the provided tables and archives. It does not yet
reconstruct all canonical amplicon, shotgun, or PMA inputs from public raw
reads. Those accession and upstream-processing records remain explicit release
gates in the data descriptor.
