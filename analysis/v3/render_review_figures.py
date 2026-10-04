#!/usr/bin/env python3
"""Regenerate corrected compartment and descriptive PMA figures with provenance."""
import argparse
import json
from pathlib import Path

import pandas as pd

import make_submission_figures as figures


def campaign_matched_alpha(table):
    """Adapt the approved common-profile estimand to the plotting schema."""
    selected = table.loc[
        table["cohort"].eq("common_profiles")
        & table["estimand"].eq("campaign_matched")
    ].copy()
    expected = {"Deep-Surface", "Rhizosphere-Surface", "Rhizosphere-Deep"}
    for metric in ("shannon", "normalized_shannon"):
        rows = selected.loc[selected["metric"].eq(metric)]
        if len(rows) != 3 or set(rows["contrast"]) != expected:
            raise ValueError(f"Incomplete campaign-matched figure cohort: {metric}")
    selected["metric"] = selected["metric"].replace(
        {"normalized_shannon": "evenness_h_over_log_hurlbert"}
    )
    selected["trip"] = "all"
    selected["comparison"] = selected["contrast"]
    selected["mean_difference"] = selected["mean"]
    selected["bootstrap_ci_low"] = selected["ci_low"]
    selected["bootstrap_ci_high"] = selected["ci_high"]
    return selected, selected.copy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--landscape", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source = root / "analysis/v3"
    inputs = {
        "matched_alpha": source / "paired_alpha_sensitivity_20260909/paired_alpha_estimands.tsv",
        "location": source / "compartment_composition/compartment_location_results.tsv",
        "loadings": source / "compartment_composition/paired_displacement_loadings.tsv",
        "pathways": source / "picrust2_ecology/position_profile_tests.tsv",
        "ko_profiles": source / "measured_function_summary_results/per_sample_ko_correlations.tsv",
        "ko_metrics": source / "measured_function_summary_results/summary_metrics.tsv",
        "pma": source / "pma_endpoint_results/pma_pair_endpoints.tsv",
        "controls": source / "control_audit/trip5_removal_fraction_by_profile.tsv",
    }
    if args.landscape:
        inputs.update({
            "alpha": root / "analysis/v2/review/cache/alpha.tsv",
            "coordinates": source / "spatial_turnover_rescue/results/site_coordinates.tsv",
            "distance_pairs": source / "distance_decay_turnover/distance_decay_pairs.tsv",
            "climate_site": source / "environment_associations/climate_site_summary.tsv",
            "climate_alpha": source / "environment_associations/climate_alpha_correlations.tsv",
            "climate_genus": source / "environment_associations/climate_genus_correlations.tsv",
            "boundary": root / "metadata/geodata/empty_quarter_boundary.kml",
            "background": root / "metadata/geodata/bluemarble_arabia_200407_120ppd.png",
            "background_metadata": root / "metadata/geodata/bluemarble_arabia_200407_120ppd.json",
        })
    runtime = figures.require_figure_runtime()
    frames = {key: pd.read_csv(path, sep="\t", index_col=0 if key == "alpha" else None)
              for key, path in inputs.items() if path.suffix == ".tsv"}
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures.setup_style()
    paired, evenness = campaign_matched_alpha(frames["matched_alpha"])
    paths = [output / "fig2_soil_position.pdf", output / "fig3_function_controls.pdf"]
    figures.make_soil_position_figure(
        paired, evenness, frames["location"], frames["loadings"], paths[0]
    )
    figures.make_function_control_figure(
        frames["pathways"], frames["ko_profiles"], frames["ko_metrics"],
        frames["pma"], {}, frames["controls"], paths[1]
    )
    if args.landscape:
        landscape = output / "fig1_landscape.pdf"
        site = frames["climate_site"].drop(columns=["latitude", "longitude", "transect_km"])
        site = site.merge(frames["coordinates"], on="site", validate="one_to_one")
        figures.make_landscape_figure(
            frames["alpha"], frames["coordinates"], inputs["boundary"], inputs["background"],
            frames["distance_pairs"], site, frames["climate_alpha"], frames["climate_genus"], landscape)
        paths.append(landscape)
    manifest = {
        "runtime": runtime,
        "scope": "Figure 2 campaign-matched common-profile diversity contrasts; Figure 3 descriptive PMA groups at two campsites; optional Figure 1 two-column layout and corrected route coordinates with frozen climate exposures",
        "inputs": {key: {"path": str(path.relative_to(root)), "sha256": figures.sha256(path)} for key, path in inputs.items()},
        "generators": {str(path.relative_to(root)): figures.sha256(path) for path in [Path(__file__).resolve(), Path(figures.__file__).resolve()]},
        "outputs": {path.name: {"sha256": figures.sha256(path), "bytes": path.stat().st_size} for path in paths},
    }
    (output / "figure_review_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    legacy = output / "figure_manifest.tsv"
    if legacy.exists():
        table = pd.read_csv(legacy, sep="\t")
        # Replace retired alpha inputs as well as the output hashes. The JSON
        # manifest remains the complete custody record for this renderer.
        table = table.loc[~table["name"].isin([
            "paired_compartment_effects", "paired_evenness_effects", "campaign_matched_alpha"
        ])].copy()
        matched = inputs["matched_alpha"]
        table = pd.concat([table, pd.DataFrame([{
            "role": "input", "name": "campaign_matched_alpha", "file": matched.name,
            "bytes": matched.stat().st_size, "sha256": figures.sha256(matched),
        }])], ignore_index=True)
        for path in paths:
            selected = table["role"].eq("output") & table["file"].eq(path.name)
            table.loc[selected, "sha256"] = figures.sha256(path)
            table.loc[selected, "bytes"] = path.stat().st_size
        table.to_csv(legacy, sep="\t", index=False)


if __name__ == "__main__":
    main()
