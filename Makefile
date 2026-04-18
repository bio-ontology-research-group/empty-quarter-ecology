# Empty Quarter amplicon — reproducibility pipeline
# Run `make help` for a tour.

PY := .venv/bin/python
QUARTO := QUARTO_PYTHON=$(PY) ~/.local/bin/quarto

NOTEBOOKS := 00_load_and_qc 01_scale_and_phyla 02_assembly \
             03_temporal 04_depth 05_distance_decay 06_function \
             07_network 08_csp_mag \
             09_causal_tier1 10_causal_tier2 11_causal_tier3 \
             12_causal_nonlinear \
             13_vegetation_mediation 14_keystone_knockout \
             15_intervention_scenarios 16_cross_desert \
             17_xrf_compartment \
             99_audit

NB_OUTPUTS := $(addprefix _output/notebooks/,$(addsuffix .html,$(NOTEBOOKS)))

.PHONY: help env cache figures check lock clean

help:  ## print this help
	@awk 'BEGIN{FS=":.*##"} /^[a-z_-]+:.*##/ {printf "%-12s %s\n",$$1,$$2}' $(MAKEFILE_LIST)

env:  ## create .venv and install all deps
	uv venv --python 3.11 .venv
	uv pip install --python $(PY) -e .
	uv pip install --python $(PY) jupyter nbformat ipykernel jupyter-cache

cache: cache/feature_table.parquet  ## run QC notebook and cache analysis-ready tables

cache/feature_table.parquet: notebooks/00_load_and_qc.qmd src/eq/loader.py src/eq/sample_id.py
	$(QUARTO) render notebooks/00_load_and_qc.qmd

figures: cache/feature_table.parquet  ## render all figure notebooks
	@for nb in 01_scale_and_phyla 02_assembly 03_temporal 04_depth \
	           05_distance_decay 06_function 07_network 08_csp_mag \
	           09_causal_tier1 10_causal_tier2 11_causal_tier3 \
	           12_causal_nonlinear 13_vegetation_mediation \
	           14_keystone_knockout 15_intervention_scenarios \
	           16_cross_desert 17_xrf_compartment \
	           99_audit; do \
	   echo "=== rendering $$nb ==="; \
	   $(QUARTO) render notebooks/$$nb.qmd || exit 1; \
	done

lock:  ## generate/refresh uv.lock
	uv lock --python $(PY)

check:  ## sanity checks on caches and figures
	@echo "== cached tables =="
	@ls -lh cache/*.parquet cache/*.json 2>/dev/null || echo "cache/ incomplete"
	@echo "== figures =="
	@bash -c 'ls figures/fig*.pdf 2>/dev/null | wc -l'
	@echo "(expected 18+ figure PDFs — main-text figs 1-9e + 10,14,15 + supp)"
	@echo "== notebook HTML =="
	@ls _output/notebooks/*.html 2>/dev/null | wc -l
	@echo "(expected 17 HTML outputs: 00-15 + 99_audit)"

clean:  ## remove rendered notebooks and caches; keep raw inputs
	rm -rf _output/ _freeze/ .quarto/
	rm -f cache/*.parquet cache/*.json cache/*.tsv
	@echo "cleaned rendered output + caches"
