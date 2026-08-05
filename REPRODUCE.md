# Reproducing the ecology manuscript

The ecology paper uses two levels of verification. The short path checks the
committed claim evidence, rerenders every submitted figure, and rebuilds both
PDFs. The complete path starts from the checksum-pinned inputs in the data
repository and reruns the analysis and knowledge-graph workflow on an approved
remote host.

## 1. Obtain the exact repositories

```bash
git clone git@github.com:bio-ontology-research-group/empty-quarter-data-paper.git
git clone git@github.com:bio-ontology-research-group/empty-quarter-ecology-reproducibility.git
cd empty-quarter-data-paper
git checkout b657df1b82812bfe1d84b37accb7d8cf7dd37878
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

The wrapper refuses every hostname except `leechuck-office` (`ws`) and
`cbontsr01` (Ontolinator). It verifies the repository lock, runs both test
suites, regenerates the core and advanced ecology analyses, reruns the control,
rainfall, function, network, pH, XRF, and knowledge-graph stages, validates the
full taxonomy ABox, regenerates figures, and builds both manuscripts. Nextflow
writes a trace, report, timeline, DAG, source-state record, environment record,
commands, logs, and SHA-256 manifests into a new output directory.

Do not use `-resume` after changing source code, input data, manuscript text, or
the data lock. A local Nextflow stub is useful only for checking wiring and is
not scientific or semantic validation.

## 5. Reproducibility boundary

The release can reconstruct the submitted analyses, figures, manuscript, and
derived knowledge graph from the provided tables and archives. It does not yet
reconstruct all canonical amplicon, shotgun, or PMA inputs from public raw
reads. Those accession and upstream-processing records remain explicit release
gates in the data descriptor.
