SHELL := /usr/bin/env bash
PYTHON ?= python3
DATA_REPO ?= ../empty-quarter-data-paper
SOURCE_DATE_EPOCH ?= 1785888000

.PHONY: bootstrap manifest verify test figures paper clean well-analysis

bootstrap:
	bash scripts/release/bootstrap_data_dependency.sh "$(DATA_REPO)"

manifest:
	$(PYTHON) scripts/release/build_repository_manifest.py . --write

verify:
	$(PYTHON) scripts/release/verify_repository.py .

test: bootstrap
	EQ_DATA_REPO="$$(realpath "$(DATA_REPO)")" \
		$(PYTHON) -m pytest -q tests

figures:
	$(PYTHON) scripts/release/render_figures.py .

paper:
	cd empty-quarter-amplicon && \
		SOURCE_DATE_EPOCH=$(SOURCE_DATE_EPOCH) FORCE_SOURCE_DATE=1 \
		latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	cd empty-quarter-amplicon && \
		SOURCE_DATE_EPOCH=$(SOURCE_DATE_EPOCH) FORCE_SOURCE_DATE=1 \
		latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex

clean:
	cd empty-quarter-amplicon && latexmk -C main.tex && latexmk -C supplement.tex

well-analysis: bootstrap
	$(PYTHON) analysis/well-adjacency-20261005/01_build_bc.py \
		--positions analysis/well-adjacency-20261005/well_positions.csv \
		--metadata analysis/batch-adjacency-2026-09-07/batch_meta.tsv \
		--features "$(DATA_REPO)/metadata/taxonomy/feature-table-trips1-5.tsv" \
		--output analysis/well-adjacency-20261005/outputs
	$(PYTHON) analysis/well-adjacency-20261005/02_neighbor_permtest.py \
		--metadata analysis/batch-adjacency-2026-09-07/batch_meta.tsv \
		--output analysis/well-adjacency-20261005/outputs
