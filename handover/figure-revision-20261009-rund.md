# Figure revision handoff — Empty Quarter ecology manuscript

Prepared 9 October 2026 against `bio-ontology-research-group/empty-quarter-ecology`,
branch `main`, commit `64b40f54b856a4dea7cb9d2ce56a7afa20e311a9`.

This folder carries the code changes for five revised manuscript figures. Nothing
here has been committed. The changes were previewed outside the pinned
environment, so regenerate every output in the pinned environment before release.

## Contents

| Path | What it is |
|---|---|
| `figure_changes.patch` | All code changes as one patch (four edited files, one new file). Verified to apply cleanly to the commit above. |
| `files/` | The same five files in full, at their repository paths, in case the patch does not apply. |
| `previews/fig1.png` … `fig5.png` | Reference renders of the intended result. Layout references only. |

Files changed by the patch:

- `analysis/v3/make_submission_figures.py` (edited)
- `analysis/v3/render_review_figures.py` (edited: two call signatures)
- `analysis/v3/aitchison_ordination.py` (new)
- `scripts/figures/prepare_bluemarble_crop.py` (edited: crop extent)
- `tests/test_taxon_context.py` (edited: expected crop extent)

## Apply

Run from the repository root.

```sh
# 1. Code
git apply --check figure_handoff/figure_changes.patch
git apply figure_handoff/figure_changes.patch

# 2. Download the Blue Marble source tile (88 MB; do not commit it)
curl -L -o metadata/geodata/world.topo.bathy.200407.3x21600x21600.C1.jpg \
  "https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73751/world.topo.bathy.200407.3x21600x21600.C1.jpg"

# 3. Regenerate the wider satellite crop (38-63 E, 16-25 N, 120 pixels per degree)
python scripts/figures/prepare_bluemarble_crop.py \
  metadata/geodata/world.topo.bathy.200407.3x21600x21600.C1.jpg

# 4. Compute the Aitchison ordination used by Figure 2
python analysis/v3/aitchison_ordination.py \
  --output-dir analysis/v3/aitchison_ordination

# 5. Build the figures
python analysis/v3/make_submission_figures.py \
  --core-dir <core results dir> --output-dir <figure output dir>
```

Step 2: the crop script verifies the tile against SHA-256
`ee8490ab1eb35d620d8d1ad8e69b3234c0b050e4eddb80e7232a2d165e475aa0` and stops if
it differs. The tile is untracked and not git-ignored, so keep it out of commits.

Step 3 overwrites `metadata/geodata/bluemarble_arabia_200407_120ppd.png` and its
`.json` sidecar (3,000 x 1,080 pixels, about 3.5 MB; previously 1,560 x 1,080).

Step 4 expected result: 630 profiles, 60 sites, 200 genera, PC1 14.1% and PC2 9.6%
of variance. It writes `ordination_scores.tsv`, `ordination_loadings.tsv` and
`ordination_summary.json`.

Step 5 requires Python 3.11.14, matplotlib 3.9.4 and FreeType 2.14.3; the script
refuses to run otherwise. New optional arguments, each with a default under
`analysis/v3/`: `--ordination-dir` (`aitchison_ordination`), `--ph-dir`
(`ph_group_linkage_20260909`), `--xrf-dir` (`xrf_community_rescue`),
`--biology-dir` (`biology_context_corrected_20260909`).

## What each figure should look like

All panels carry a bold lower-case letter at the top left and no title.

| Figure | Output file | Function | Panels |
|---|---|---|---|
| 1 | `fig1_landscape.pdf` | `make_landscape_figure` | a: map across the full width, 38.9-61.7 E, whole Rub' al-Khali outline. b: samples per campaign by compartment, no numbers on bars. c: climate along the transect. |
| 2 | `fig_composition_geography.pdf` | `make_composition_geography_figure` (new) | a: Aitchison ordination, one panel per compartment, coloured by transect coordinate, all profiles in grey behind. b: distance decay, full width. |
| 3 | `fig2_soil_position.pdf` | `make_soil_position_figure` | Unchanged apart from letters a-d and no titles. |
| 4 | `fig4_environment_gradients.pdf` | `make_environment_gradient_figure` (new) | a: soil pH, site means per compartment. b: laboratory XRF elemental PC1, site means per compartment. c: Shannon diversity per site by landform. d: climate-diversity associations. e: climate-associated genera. |
| 5 | `fig3_function_controls.pdf` | `make_function_control_figure` | a: predicted pathway contrasts. b: PMA aliquot outcomes. c: contaminant reads removed. The gene-family rank agreement panel was removed. |

## Open items for whoever implements this

1. **Output file names do not match the new figure order** (see the table). Decide
   the final numbering, then rename the outputs and update `main.tex`, the
   manifests and the tests that cite the file names.
2. **Tests that will fail until updated:**
   - `tests/test_ecology_rewrite_claims.py::test_landscape_figure_contains_six_evidence_bearing_panels`
     expects the old six-panel Figure 1 and its titles.
   - `tests/test_cross_desert_context.py::test_soil_position_figure_adds_descriptive_taxon_context`
     expects the removed panel title "Genera contributing most to compartment differences".
   - `tests/test_taxon_context.py::test_landscape_figure_uses_the_committed_satellite_crop`
     compares the crop checksum with `figure_review_manifest.json`, which is stale
     until `render_review_figures.py` is rerun.
3. **Manifests are stale:** `FILE_MANIFEST.tsv` and
   `empty-quarter-amplicon/figures/figure_review_manifest.json` hold old checksums.
4. **Manuscript not touched.** Captions and in-text panel references in
   `empty-quarter-amplicon/main.tex` still describe the old figures. Captions must
   now carry what the removed titles said, including that removed reads in
   Figure 5c are *candidate* contaminants.
5. **Unused inputs.** `main()` still loads and validates the two gene-family
   (KO) tables, which no figure uses now.
6. **No tests yet** for `aitchison_ordination.py` or the two new figure functions.
7. **Design choices to confirm with the author:** PC1 is oriented to increase
   eastwards; Figure 4a-b average over campaigns inside the figure script;
   Figure 4c pools the four single-site landforms as "Other"; Figure 4 uses the
   `ph_group_linkage_20260909` pH table.
