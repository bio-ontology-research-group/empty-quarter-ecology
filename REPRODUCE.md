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
git clone https://github.com/bio-ontology-research-group/empty-quarter-ecology.git
git clone https://github.com/bio-ontology-research-group/empty-quarter-data-paper.git
cd empty-quarter-data-paper
git checkout "$(awk -F '\t' '$1 == "commit" {print $2}' ../empty-quarter-ecology/DATA_REPOSITORY.lock)"
bash scripts/release/download_bulk_artifacts.sh
bash scripts/release/bootstrap_package_layout.sh .
cd ../empty-quarter-ecology
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
cd ../empty-quarter-ecology
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
cd empty-quarter-ecology
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

Current coordinate-dependent results use Site 52 at 20.82784 N, 53.57835 E.
The October finalization replays the XRF/environment associations, control-filter
sensitivity and batch-adjacency analyses on that geometry. The shared data pin packages the current rainfall inputs under
`metadata/climate/current_analysis_inputs.json`. Its corrected daily Open-Meteo
series is byte-identical to the ecology input at
`analysis/v3/open_meteo_site52_corrected_20261004/daily_weather_canonical_site52_corrected.tsv`.
Its package records the original source hash, corrected request and response,
and verifies that every non-Site-52 row is unchanged. The rainfall refit and
its current sensitivity consumers use this replacement. The five-product Trip 1
comparison uses the corrected grid-cell extraction in
`analysis/v3/rain_event_product_exposures_site52_corrected_20261004/`,
rebuilt from the 21 hash-verified original climate files and also packaged in
the shared data pin. Frozen historical climate inputs retain their acquisition
provenance. The companion data descriptor documents KG v3.0.1
(https://doi.org/10.5281/zenodo.23134168).

The reported rainfall model is `analysis/v3/rain_calendar_refit_20260909/`,
including its complete six-member calendar-year permutation orbit and full
site-bootstrap refits. Other rainfall folders contain diagnostic or superseded
models; they do not replace the reported calendar-randomization analysis.
Likewise, `biology_context/` and `taxon_context/` are superseded by their
`*_corrected_20260909/` counterparts. The latter retain their original run
hashes plus an explicit reconciliation proving that the climate columns
actually consumed are identical after the coordinate-column update.
The pH partition is independently reproducible with
`bash analysis/v3/ph_partition_20261004/reproduce.sh /tmp/ph-reproduction`.
Its packaged inputs include the corrected coordinates; the older lock under
`provenance/` is explicitly historical acquisition context.

Do not use `-resume` after changing source code, input data, manuscript text, or
the data lock. A local Nextflow stub is useful only for checking wiring and is
not scientific or semantic validation.

## 5. Reproducibility boundary

The release can reconstruct the submitted analyses, figures, manuscript, and
derived knowledge graph from the provided tables and archives. It does not yet
reconstruct all canonical amplicon, shotgun, or PMA inputs from public raw
reads. Those accession and upstream-processing records remain explicit release
gates in the data descriptor.

## 6. Review and finalization record

The manuscript in `empty-quarter-amplicon/` is the current Overleaf text, not an
older snapshot. `empty-quarter-amplicon/validation/finalization-20261004/`
contains the complete change list and redlines against Rund's 2 October
Overleaf revision (`a5d5623`). The finalization record identifies input versions,
reruns and validation. `scripts/release/build_rund_redlines.py` generates the
text comparisons from an exported baseline and built current manuscript.

One author-source verification remains: the physical well-to-well analysis
(259 T1 samples, p=0.87) lacks a deposited plate map, code and output. Its
manuscript text is preserved, and `OPEN_VERIFICATION.txt` in the review package
identifies the required inputs. The reproduced library-order adjacency test
is different. A passing test suite is not verification of that missing analysis.
