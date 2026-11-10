# Genome-Resolved Functional and Canonical Taxonomy Integration (`analysis/genome_integration_20261108`)

This directory contains the Stage 4 genome-resolved ecological synthesis and the canonical GTDB release 232 genus-turnover sensitivity analysis integrated into `empty-quarter-amplicon/alt_stats.tex` and `empty-quarter-amplicon/supp_altstat.tex`.

## Directory Layout

- `code/common_science.py`: Shared primary-key validation, symmetric three-component genus-turnover decomposition (`total_diff = between_shared + within_shared + exclusive`), and spatial-block bootstrap helpers.
- `code/stage4_ecology_synthesis.py`: Stage 4 synthesis generating the 41-rule microTrait/KO functional matrices (`task_d_h6_baseline_vs_protection_surface_deep.tsv`, `task_d_all_41_microtrait_and_curated_matrix.tsv`, `task_d_h6_h3_systematic_matrix.tsv`), denominator comparison (`task_d_turnover_allread_vs_mapped_comparison.tsv`), and conditional eastern megabarchan-versus-uruq joint diversity contrast sensitivity (`task_d_eastern_joint_diversity_contrast.tsv`).
- `code/canonical_taxonomy_sensitivity.py`: Exact 1:1 reconciliation between the frozen analytical taxonomy snapshot (`EQ-TAX-LEGACY-ANALYTICAL-v0`, `genome_quality.tsv`) and canonical `gtdb_r232_classification` (`EQ-MAG-REFINED-v1.0.0`, `refined_representatives.tsv`), plus the 10-contiguous-spatial-block (1,000 draws, seed `20261108`) and individual-site genus-turnover decomposition replay for the three focal guilds (`HMM:glycine betaine transporter` [BetL], `HMM:ahpCF` [AhpCF], and `typed:[NiFe]_Group_1l`).
- `results/`:
  - `endpoints.tsv`: Sample-level diversity, trait, and marker endpoints (`n = 546` core-transect shotgun libraries; `n = 545` read-depth-modeled libraries).
  - `PICRUSt_MAG_effect_concordance.tsv`: Paired compartment concordance between PICRUSt2 and MAG-derived KO profiles (`EQ-KO-LEGACY-ANALYTICAL-v0`, `rep_ko_v2.tsv`) across `8,860` shared KOs in the `442` specimen-matched libraries.
  - `spatial_fixed_effects.tsv` and `spatial_fit_diagnostics.tsv`: Primary and sensitivity R-INLA SPDE posterior summaries and numerical audit across `1,800` fits (`207` responses).
  - `turnover_summary.tsv`: Frozen three-component genus-turnover summary under `EQ-TAX-LEGACY-ANALYTICAL-v0`.
  - `task_d_h6_baseline_vs_protection_surface_deep.tsv`: Selected surface-minus-subsurface functional panel (Supplementary Table `tab:stage4_selected_panel`).
  - `task_d_all_41_microtrait_and_curated_matrix.tsv` and `task_d_h6_h3_systematic_matrix.tsv`: Complete 41-trait synthesis matrices across compartments, denominators, and cohorts.
  - `task_d_turnover_allread_vs_mapped_comparison.tsv`: Mapped-community versus all-read denominator turnover comparison.
  - `task_d_eastern_joint_diversity_contrast.tsv`: Conditional joint spatial-block bootstrap sensitivity for eastern megabarchan versus uruq diversity slopes (Supplementary Table `tab:stage4_joint_diversity`; `685/1,000` sample-level and `767/1,000` site--compartment resamples identifiable).
  - `task_d_panel_and_multiplicity_manifest.json` and `summary.json`: Multiplicity family definitions and Stage 4 synthesis summary metadata.
  - `task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv`: Completed canonical GTDB release 232 versus earlier analytical taxonomy genus-turnover decomposition for the three focal guilds (Supplementary Table `tab:canonical_turnover_3guilds`; primary rows filter `is_primary_required_3guild == True` and `denominator_mode == "mapped_community_normalized_weights"`).
  - `taxonomy_reconciliation_summary.json`: Exact 5,438-genome and 237-Thermoplasmatota taxonomy reconciliation summary (`1,898` whole strings and `1,124` genus labels differ; `0` domain differences).
  - `canonical_root_growth_taxonomy_status.json`: Status record distinguishing the completed canonical three-guild genus-turnover decomposition from the exploratory order- and genus-adjusted `gRodon2` root-depletion regression, which remains under `EQ-TAX-LEGACY-ANALYTICAL-v0` (`PENDING: TAXONOMY-SNAPSHOT-SENSITIVITY`).
- `tests/`:
  - `test_science_20261108.py`: Unit and regression tests for `common_science.py` and `stage4_ecology_synthesis.py`.
  - `test_canonical_taxonomy_sensitivity.py`: Unit, adversarial, and staged-output verification tests for `canonical_taxonomy_sensitivity.py` and `task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv`.

## Canonical Replay vs Frozen 41-Trait Panel and Pending `gRodon2` Regression

1. **Completed canonical three-guild turnover replay (`task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv`)**:
   - Reruns the exact symmetric genus-turnover decomposition (`total_diff = between_shared + within_shared + exclusive`) under canonical `gtdb_r232_classification` against the immutable legacy baseline (`EQ-TAX-LEGACY-ANALYTICAL-v0`) using identical CoverM abundances (`sample_by_mag_coverm_relative_abundance.npz`), genome trait calls (`genome_traits.tsv`), and paired visits (`68` visits across `36` sites for Surface minus Deep; `106` visits across `47` sites for Root-adjacent minus Deep), with `1,000` resamples across the `10` contiguous spatial blocks (`transect_km`, seed `20261108`) and individual-site bootstrap intervals as a secondary sensitivity.
   - Genomes without a named genus (`g__`) are kept as separate singletons (`unclassified_g___genomes_kept_as_separate_singletons_never_pooled`), are not treated as independent biological lineages, and have zero within-group change by construction.
2. **Frozen 41-trait panel (`task_d_h6_baseline_vs_protection_surface_deep.tsv`, Supplementary Table `tab:stage4_selected_panel`)**:
   - Retains the original paired trait bootstrap intervals and frozen `EQ-TAX-LEGACY-ANALYTICAL-v0` decomposition components for non-focal traits without re-running non-focal trait decompositions under canonical taxonomy.
3. **Pending `gRodon2` order/genus regression (`canonical_root_growth_taxonomy_status.json`)**:
   - The exploratory order- and genus-adjusted `gRodon2` growth-potential regression for group 1l root depletion remains reported under the earlier taxonomy snapshot (`EQ-TAX-LEGACY-ANALYTICAL-v0`) with a visible `TAXONOMY-SNAPSHOT-SENSITIVITY` marker pending canonical replay.

## Verification Commands

Run the Stage 4 and canonical taxonomy unit and output verification tests from the repository root:

```bash
uv run --python .venv/bin/python -m pytest -q \
  analysis/genome_integration_20261108/tests/test_science_20261108.py \
  analysis/genome_integration_20261108/tests/test_canonical_taxonomy_sensitivity.py \
  tests/test_active_alt_stats_integration.py
```

To replay the canonical three-guild turnover decomposition from upstream inputs:

```bash
uv run --python .venv/bin/python analysis/genome_integration_20261108/code/canonical_taxonomy_sensitivity.py \
  --mode turnover_3guilds_only \
  --orig-arch-root "$EQ_ORIG_ARCH_ROOT" \
  --orig-sci-out "$EQ_ORIG_SCIENCE_OUT" \
  --canon-out "$EQ_CANONICAL_SCIENCE_OUT" \
  --canonical-rep-tsv "$EQ_DATA_REPO/metadata/metagenome/canonical/refined_representatives.tsv"
```
