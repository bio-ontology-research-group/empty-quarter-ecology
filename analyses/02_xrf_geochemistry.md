# 02. XRF geochemistry baseline

**Question.** What is the spatial structure of soil chemistry across the Empty Quarter?

**Method.** Portable XRF on 725 samples covering all 60 sites × 3 compartments. Lithology PCA, per-element Shannon (chemodiversity), site×compartment panels. `scripts/run_xrf_analysis.py`, `scripts/build_xrf_combined.py`, `scripts/render_supplement_xrf_figures.py`.

**Inputs.**
- `data/geochemistry/xrf_lab_table_all_trips.tsv` (raw)
- `cache/xrf_summary_all_trips.tsv` (canonical aggregate)

**Key results.**
- Bimodal chemistry: **sabkha** (high SO₃, Na, Cl, Ca, Mg — evaporite) vs **sandy** (low salt, higher Si).
- Lithology PCA separates compartments along PC1 (Si–Al silicate vs Ca–S evaporite axis).
- Per-element Shannon (chemodiversity) is highest in rhizosphere and tracks Si:Ca ratio (RQ10/RQ16 supplement figures).
- Mn loading separates desert-varnish-rich surface vs deeper layers (RQ11).
- S concentration correlates with SRB/SOB ratios (RQ12).

**Outputs.**
- `cache/xrf_lithology_pca.tsv`
- `cache/xrf_chemodiversity.tsv`
- `cache/xrf_per_compartment.tsv`
- `cache/xrf_site_compartment_panel.tsv`
- `cache/per_element_shannon.tsv`
- `cache/supplement_xrf_stats.txt`
- Supplement figures `supp_RQ10_*`, `supp_RQ11_*`, `supp_RQ12_*`, `supp_RQ16_*`

**Status.** solid. The "sabkha vs sandy" binary is the substrate axis on which the two-strategy hypothesis later sits.

**Cross-refs.** 37 (Nibribacter XRF substitution), 39–42 (two-strategy is sabkha-vs-sandy aligned).
