#!/usr/bin/env python3
"""Build evidence-bounded ecology-paper figures from canonical TSV outputs.

This script performs no statistical fitting. It reads the exact outputs of the
grouped ecology, spatial, climate, predicted-function, PMA and control analyses;
validates the expected records; and renders only results whose wording is
permitted by the claim verdicts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import ft2font
from matplotlib import patheffects
import numpy as np
import pandas as pd


COLORS = {
    "surface": "#E69F00",
    "deep": "#0072B2",
    "root": "#009E73",
    "positive": "#0072B2",
    "negative": "#D55E00",
    "null": "#777777",
    "observed": "#CC3311",
}
PDF_METADATA = {
    "Title": "Empty Quarter ecology manuscript evidence figure",
    "Creator": "analysis/v3/make_submission_figures.py",
    "CreationDate": None,
    "ModDate": None,
}
EXPECTED_FIGURE_RUNTIME = {
    "python": "3.11.14",
    "matplotlib": "3.9.4",
    "freetype": "2.14.3",
}


def require_figure_runtime() -> dict[str, str]:
    """Fail before rendering when native font metrics are not reproducible."""
    observed = {
        "python": sys.version.split()[0],
        "matplotlib": matplotlib.__version__,
        "freetype": ft2font.__freetype_version__,
    }
    if observed != EXPECTED_FIGURE_RUNTIME:
        raise RuntimeError(
            "Figure runtime differs from environment/conda-linux-64.lock: "
            f"expected {EXPECTED_FIGURE_RUNTIME}, observed {observed}"
        )
    return observed


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_columns(frame: pd.DataFrame, required: set[str], source: Path) -> None:
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: {sorted(missing)}")


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 11.5,
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "savefig.transparent": False,
        }
    )


def asymmetric_errors(frame: pd.DataFrame, estimate: str) -> np.ndarray:
    return np.vstack(
        [
            frame[estimate].to_numpy() - frame["ci_low"].to_numpy(),
            frame["ci_high"].to_numpy() - frame[estimate].to_numpy(),
        ]
    )


def read_kml_polygon(path: Path) -> np.ndarray:
    """Read the first polygon ring from the project orientation boundary."""
    root = ET.parse(path).getroot()
    namespace = {"kml": "http://www.opengis.net/kml/2.2"}
    node = root.find(".//kml:LinearRing/kml:coordinates", namespace)
    if node is None or not node.text:
        raise ValueError(f"No polygon coordinates found in {path}")
    coordinates = []
    for token in node.text.split():
        longitude, latitude, *_ = token.split(",")
        coordinates.append((float(longitude), float(latitude)))
    if len(coordinates) < 4:
        raise ValueError(f"Boundary in {path} has fewer than four vertices")
    return np.asarray(coordinates, dtype=float)


def read_background_image(path: Path) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Read the committed satellite crop and its geographic extent sidecar."""
    sidecar = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    extent_degrees = sidecar["extent_degrees"]
    image = mpimg.imread(path)
    extent = (
        float(extent_degrees["lon_min"]),
        float(extent_degrees["lon_max"]),
        float(extent_degrees["lat_min"]),
        float(extent_degrees["lat_max"]),
    )
    return image, extent


def make_landscape_figure(
    alpha: pd.DataFrame,
    coordinates: pd.DataFrame,
    boundary_path: Path,
    background_path: Path,
    site: pd.DataFrame,
    output: Path,
) -> None:
    """Show the transect map above sampling coverage and site climate."""
    # The map keeps true geographic proportions and spans the full figure
    # width; it extends east and west of the boundary to fill the wide panel.
    fig = plt.figure(figsize=(8.4, 6.5))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.05, 1.0))
    map_ax = fig.add_subplot(grid[0, :])
    coverage_ax = fig.add_subplot(grid[1, 0])
    climate_ax = fig.add_subplot(grid[1, 1])
    for panel_ax, panel_letter in (
        (map_ax, "a"),
        (coverage_ax, "b"),
        (climate_ax, "c"),
    ):
        panel_ax.set_title(panel_letter, loc="left", fontweight="bold")

    coordinates = coordinates.sort_values("transect_km")
    boundary = read_kml_polygon(boundary_path)
    # Show the whole boundary in latitude and twice its longitude range,
    # centred on the boundary.
    map_ylim = (boundary[:, 1].min() - 0.2, boundary[:, 1].max() + 0.2)
    map_half_width = boundary[:, 0].max() - boundary[:, 0].min() + 0.4
    map_centre = (boundary[:, 0].min() + boundary[:, 0].max()) / 2
    map_xlim = (map_centre - map_half_width, map_centre + map_half_width)
    background, extent = read_background_image(background_path)
    if (
        extent[0] > map_xlim[0]
        or extent[1] < map_xlim[1]
        or extent[2] > map_ylim[0]
        or extent[3] < map_ylim[1]
    ):
        raise ValueError(
            f"Background {background_path} covers {extent}, which does not "
            f"contain the map limits {map_xlim + map_ylim}"
        )
    map_ax.imshow(
        background,
        extent=extent,
        origin="upper",
        interpolation="none",
        aspect="auto",
        zorder=0,
    )
    map_ax.plot(
        np.append(boundary[:, 0], boundary[0, 0]),
        np.append(boundary[:, 1], boundary[0, 1]),
        color="#FFFFFF",
        linewidth=0.9,
        linestyle=(0, (3, 2)),
        zorder=1,
    )
    map_ax.text(
        51.1,
        22.7,
        "Rub' al-Khali",
        color="#FFFFFF",
        fontsize=9.0,
        ha="center",
        path_effects=[
            patheffects.withStroke(linewidth=1.6, foreground="#4A3520")
        ],
        zorder=3,
    )
    map_ax.plot(
        coordinates["longitude"],
        coordinates["latitude"],
        color="#FFFFFF",
        linewidth=0.9,
        alpha=0.85,
        zorder=1,
    )
    map_ax.scatter(
        coordinates["longitude"],
        coordinates["latitude"],
        color="#B56A2D",
        s=25,
        edgecolor="#FFFFFF",
        linewidth=0.45,
        zorder=2,
    )
    for site_number, horizontal, vertical in ((1, "left", "top"), (30, "left", "bottom"), (60, "right", "top")):
        row = coordinates.loc[coordinates["site"] == site_number]
        if len(row) != 1:
            raise ValueError(f"Expected one coordinate for site {site_number}")
        map_ax.annotate(
            f"Site {site_number}",
            (float(row.iloc[0]["longitude"]), float(row.iloc[0]["latitude"])),
            xytext=(4 if horizontal == "left" else -4, 4 if vertical == "bottom" else -4),
            textcoords="offset points",
            ha=horizontal,
            va=vertical,
            fontsize=9.0,
            color="#FFFFFF",
            path_effects=[
                patheffects.withStroke(linewidth=1.4, foreground="#4A3520")
            ],
            zorder=3,
        )
    map_ax.set_aspect(1 / np.cos(np.deg2rad(coordinates["latitude"].mean())))
    map_ax.set_xlim(*map_xlim)
    map_ax.set_ylim(*map_ylim)
    map_ax.set(xlabel="Longitude", ylabel="Latitude")

    type_order = ["Surface", "Deep", "Rhizosphere"]
    type_labels = ["Surface", "Shallow subsurface", "Root-adjacent"]
    type_colours = [COLORS["surface"], COLORS["deep"], COLORS["root"]]
    profile_counts = alpha.groupby(["Trip", "Type"]).size().unstack(fill_value=0)
    profile_counts = profile_counts.reindex(
        index=range(1, 6), columns=type_order, fill_value=0
    )
    bottom = np.zeros(len(profile_counts), dtype=float)
    for sample_type, label, colour in zip(type_order, type_labels, type_colours):
        values = profile_counts[sample_type].to_numpy()
        coverage_ax.bar(
            profile_counts.index,
            values,
            bottom=bottom,
            color=colour,
            label=label,
            width=0.72,
        )
        bottom += values
    coverage_ax.set(
        xlabel="Campaign",
        ylabel="Samples",
        xticks=range(1, 6),
        xticklabels=[
            "T1", "T2", "T3", "T4", "T5",
        ],
    )
    coverage_ax.tick_params(axis="x", labelsize=9.0)
    coverage_ax.tick_params(axis="y", labelsize=9.0)
    # Headroom keeps the legend clear of the tallest bar.
    coverage_ax.set_ylim(0, 700)
    coverage_ax.legend(frameon=False, fontsize=9.0, loc="upper right")

    climate_order = [
        "mean_air_temperature_c",
        "mean_monthly_rain_mm",
        "mean_relative_humidity_pct",
    ]
    climate_labels = ["Temperature", "Rain", "Humidity"]
    climate_colours = ["#D55E00", "#0072B2", "#009E73"]
    ordered = site.sort_values("transect_km")
    for variable, label, colour in zip(
        climate_order, climate_labels, climate_colours
    ):
        values = ordered[variable]
        standardized = (values - values.mean()) / values.std(ddof=1)
        climate_ax.plot(
            ordered["transect_km"],
            standardized,
            color=colour,
            linewidth=1.4,
            alpha=0.9,
            label=label,
        )
    climate_ax.axhline(0, color="#999999", linewidth=0.6)
    climate_ax.set(
        xlabel="West-east coordinate (km)",
        ylabel="Climate\n(σ relative to site)",
    )
    climate_ax.legend(frameon=False, fontsize=9.0)

    fig.tight_layout(h_pad=2.0, w_pad=1.8)
    fig.savefig(output, bbox_inches="tight", metadata=PDF_METADATA)
    plt.close(fig)


def make_composition_geography_figure(
    ordination: pd.DataFrame,
    ordination_summary: dict,
    distance_pairs: pd.DataFrame,
    output: Path,
) -> None:
    """Show the Aitchison ordination by compartment above distance decay."""
    fig = plt.figure(figsize=(8.4, 5.6), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.8))
    first_ax = fig.add_subplot(grid[0, 0])
    ordination_axes = [
        first_ax,
        fig.add_subplot(grid[0, 1], sharex=first_ax, sharey=first_ax),
        fig.add_subplot(grid[0, 2], sharex=first_ax, sharey=first_ax),
    ]
    distance_ax = fig.add_subplot(grid[1, :])
    for panel_ax, panel_letter in ((first_ax, "a"), (distance_ax, "b")):
        panel_ax.set_title(panel_letter, loc="left", fontweight="bold")

    variance = ordination_summary["variance_explained"]
    # One joint ordination and one colour scale; each panel highlights one
    # compartment over all profiles in grey. Marker shapes match panel b.
    colour_norm = plt.Normalize(
        ordination["transect_km"].min(), ordination["transect_km"].max()
    )
    compartment_markers = [
        ("Surface", "Surface", "o"),
        ("Deep", "Shallow subsurface", "s"),
        ("Rhizosphere", "Root-adjacent", "^"),
    ]
    unknown = set(ordination["compartment"]) - {
        compartment for compartment, _, _ in compartment_markers
    }
    if unknown:
        raise ValueError(f"Unexpected ordination compartments: {sorted(unknown)}")
    for panel_ax, (compartment, label, marker) in zip(
        ordination_axes, compartment_markers
    ):
        part = ordination[ordination["compartment"] == compartment]
        panel_ax.scatter(
            ordination["pc1"],
            ordination["pc2"],
            color="#DDDDDD",
            s=7,
            linewidth=0,
            zorder=1,
        )
        points = panel_ax.scatter(
            part["pc1"],
            part["pc2"],
            c=part["transect_km"],
            cmap="viridis",
            norm=colour_norm,
            marker=marker,
            s=15,
            edgecolor="white",
            linewidth=0.25,
            zorder=2,
        )
        panel_ax.scatter(
            [], [], marker=marker, s=15, color="#555555", label=label
        )
        panel_ax.legend(
            frameon=False,
            fontsize=9.0,
            loc="upper left",
            handletextpad=0.1,
            borderaxespad=0.1,
        )
    # Headroom keeps the compartment labels clear of the points.
    pc2_low, pc2_high = first_ax.get_ylim()
    first_ax.set_ylim(pc2_low, pc2_high + 0.12 * (pc2_high - pc2_low))
    first_ax.set_ylabel(f"PC2 ({100 * variance['pc2']:.1f}%)")
    ordination_axes[1].set_xlabel(f"PC1 ({100 * variance['pc1']:.1f}%)")
    for panel_ax in ordination_axes[1:]:
        panel_ax.tick_params(labelleft=False)
    colourbar = fig.colorbar(
        points, ax=ordination_axes, fraction=0.03, pad=0.015
    )
    colourbar.set_label("West-east coordinate (km)")

    distance_pairs = distance_pairs.copy()
    upper = float(distance_pairs["geographic_distance_km"].max())
    edges = np.linspace(0.0, upper + np.finfo(float).eps, 11)
    distance_pairs["distance_bin"] = pd.cut(
        distance_pairs["geographic_distance_km"],
        bins=edges,
        labels=False,
        include_lowest=True,
    )
    distance_order = [
        ("Surface", "Surface", COLORS["surface"], "o"),
        ("Deep", "Shallow subsurface", COLORS["deep"], "s"),
        ("Rhizosphere", "Root-adjacent", COLORS["root"], "^"),
    ]
    for compartment, label, color, marker in distance_order:
        part = distance_pairs[distance_pairs["compartment"] == compartment]
        summary = (
            part.groupby("distance_bin", observed=True)
            .agg(
                distance_km=("geographic_distance_km", "mean"),
                dissimilarity=("aitchison_dissimilarity", "mean"),
            )
            .reset_index()
        )
        distance_ax.plot(
            summary["distance_km"],
            summary["dissimilarity"],
            marker=marker,
            color=color,
            label=label,
            linewidth=1.5,
            markersize=4.2,
        )
    distance_ax.set(
        xlabel="Distance between sites (km)",
        ylabel="Aitchison dissimilarity",
    )
    distance_ax.legend(frameon=False, fontsize=9.0)

    fig.savefig(output, bbox_inches="tight", metadata=PDF_METADATA)
    plt.close(fig)


def make_environment_gradient_figure(
    ph_groups: pd.DataFrame,
    xrf_axis: pd.DataFrame,
    coordinates: pd.DataFrame,
    site: pd.DataFrame,
    landforms: pd.DataFrame,
    alpha_correlations: pd.DataFrame,
    genus_correlations: pd.DataFrame,
    output: Path,
) -> None:
    """Show soil and diversity along the transect with climate associations."""
    fig = plt.figure(figsize=(8.4, 8.4))
    grid = fig.add_gridspec(4, 2, height_ratios=(1.0, 1.0, 1.0, 2.0))
    ph_ax = fig.add_subplot(grid[0, :])
    xrf_ax = fig.add_subplot(grid[1, :], sharex=ph_ax)
    shannon_ax = fig.add_subplot(grid[2, :], sharex=ph_ax)
    diversity_ax = fig.add_subplot(grid[3, 0])
    genus_ax = fig.add_subplot(grid[3, 1])
    for panel_ax, panel_letter in zip(
        (ph_ax, xrf_ax, shannon_ax, diversity_ax, genus_ax),
        "abcde",
    ):
        panel_ax.set_title(panel_letter, loc="left", fontweight="bold")

    # Compartment colours and markers match the distance-decay panel.
    compartment_styles = [
        ("Surface", "Surface", COLORS["surface"], "o"),
        ("Deep", "Shallow subsurface", COLORS["deep"], "s"),
        ("Rhizosphere", "Root-adjacent", COLORS["root"], "^"),
    ]
    transect_by_site = coordinates.set_index("site")["transect_km"]
    ph_means = ph_groups.groupby(["site", "compartment"], as_index=False)[
        "ph"
    ].mean()
    xrf_means = (
        xrf_axis.rename(columns={"Site": "site", "Type": "compartment"})
        .groupby(["site", "compartment"], as_index=False)["elemental_pc1"]
        .mean()
    )
    for panel_ax, means, value in (
        (ph_ax, ph_means, "ph"),
        (xrf_ax, xrf_means, "elemental_pc1"),
    ):
        unknown = set(means["compartment"]) - {
            compartment for compartment, _, _, _ in compartment_styles
        }
        if unknown or not set(means["site"]) <= set(transect_by_site.index):
            raise ValueError(f"Unexpected site or compartment in {value} table")
        for compartment, label, colour, marker in compartment_styles:
            part = means[means["compartment"] == compartment]
            panel_ax.scatter(
                part["site"].map(transect_by_site),
                part[value],
                color=colour,
                marker=marker,
                s=16,
                alpha=0.85,
                edgecolor="white",
                linewidth=0.3,
                label=label,
            )
    ph_ax.set_ylabel("Soil pH")
    ph_ax.legend(
        frameon=False,
        fontsize=9.0,
        ncol=3,
        loc="upper left",
        handletextpad=0.1,
        columnspacing=1.0,
    )
    xrf_ax.set_ylabel("Elemental PC1")

    climate_order = [
        "mean_air_temperature_c",
        "mean_monthly_rain_mm",
        "mean_relative_humidity_pct",
    ]
    climate_labels = ["Temperature", "Rain", "Humidity"]
    climate_colours = ["#D55E00", "#0072B2", "#009E73"]

    # Landforms recorded at a single site are pooled so that each colour
    # stands for at least three sites.
    landform_styles = [
        ("sand dune", "Sand dune", "#B08D3C"),
        ("saline pan", "Saline pan", "#56B4E9"),
        ("desert oasis", "Desert oasis", "#CC79A7"),
    ]
    site_landform = site.merge(
        landforms[["site", "landform"]], on="site", validate="one_to_one"
    )
    if len(site_landform) != len(site):
        raise ValueError("A site has no landform record")
    named = {landform for landform, _, _ in landform_styles}
    for landform, label, colour in landform_styles:
        part = site_landform[site_landform["landform"] == landform]
        shannon_ax.scatter(
            part["transect_km"],
            part["shannon"],
            color=colour,
            s=20,
            edgecolor="white",
            linewidth=0.3,
            label=label,
        )
    other = site_landform[~site_landform["landform"].isin(named)]
    if other["landform"].value_counts().max() > 1:
        raise ValueError("A landform with several sites has no colour")
    shannon_ax.scatter(
        other["transect_km"],
        other["shannon"],
        color="#777777",
        s=20,
        edgecolor="white",
        linewidth=0.3,
        label="Other",
    )
    shannon_ax.set(
        xlabel="West-east coordinate (km)", ylabel="Shannon diversity"
    )
    shannon_ax.legend(
        frameon=False,
        fontsize=9.0,
        ncol=4,
        loc="lower left",
        handletextpad=0.1,
        columnspacing=1.0,
    )
    for panel_ax in (ph_ax, xrf_ax):
        panel_ax.tick_params(labelbottom=False)

    response_order = ["shannon", "expected_richness_25k", "normalized_evenness"]
    response_labels = ["Shannon", "Expected\nrichness", "Norm.\nShannon"]
    matrix = (
        alpha_correlations.pivot(
            index="climate_variable", columns="response", values="spearman_rho"
        )
        .loc[climate_order, response_order]
        .to_numpy(dtype=float)
    )
    image_plot = diversity_ax.imshow(
        matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto"
    )
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            diversity_ax.text(
                column,
                row,
                f"{matrix[row, column]:.2f}",
                ha="center",
                va="center",
                color="white" if abs(matrix[row, column]) > 0.55 else "#222222",
                fontsize=9.0,
            )
    diversity_ax.set(
        xticks=np.arange(3),
        xticklabels=response_labels,
        yticks=np.arange(3),
        yticklabels=climate_labels,
    )
    diversity_ax.tick_params(axis="x", labelsize=9.0, rotation=25)
    for label in diversity_ax.get_xticklabels():
        label.set_ha("right")
    colourbar = fig.colorbar(image_plot, ax=diversity_ax, fraction=0.05, pad=0.04)
    colourbar.set_label("Rank association (Spearman $\\rho$)")

    supported = genus_correlations[
        genus_correlations["supported_q_lt_0_05"]
    ].copy()
    selected_genera: set[str] = set()
    for variable in climate_order:
        part = supported[supported["climate_variable"] == variable]
        selected_genera.update(part.nlargest(2, "spearman_rho")["genus"])
        selected_genera.update(part.nsmallest(2, "spearman_rho")["genus"])
    selected = supported[supported["genus"].isin(selected_genera)].pivot(
        index="genus", columns="climate_variable", values="spearman_rho"
    )
    if selected.shape != (10, 3) or selected.isna().any().any():
        raise ValueError(
            "Expected ten fully supported genera from the two-tail selection rule"
        )
    selected["mean_rho"] = selected[climate_order].mean(axis=1)
    selected = selected.sort_values("mean_rho")
    genus_y = np.arange(len(selected))
    offsets = [-0.18, 0.0, 0.18]
    for variable, label, colour, offset in zip(
        climate_order, climate_labels, climate_colours, offsets
    ):
        genus_ax.scatter(
            selected[variable],
            genus_y + offset,
            color=colour,
            s=23,
            label=label,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
    genus_ax.axvline(0, color="#777777", linewidth=0.8)
    genus_ax.set(
        xlabel="Genus–climate rank association\n(Spearman $\\rho$)",
        yticks=genus_y,
        yticklabels=[rf"$\it{{{name}}}$" for name in selected.index],
        xlim=(-0.9, 0.9),
    )
    genus_ax.legend(frameon=False, fontsize=9.0, loc="lower right")

    fig.tight_layout(h_pad=0.5, w_pad=1.8)
    fig.align_ylabels([ph_ax, xrf_ax, shannon_ax])
    fig.savefig(output, bbox_inches="tight", metadata=PDF_METADATA)
    plt.close(fig)


def make_soil_position_figure(
    paired: pd.DataFrame,
    evenness: pd.DataFrame,
    location: pd.DataFrame,
    loadings: pd.DataFrame,
    output: Path,
) -> None:
    """Plot paired soil-position effects and their leading taxon contrasts."""
    fig = plt.figure(figsize=(10.6, 6.1))
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 1.15))
    axes = [fig.add_subplot(grid[0, column]) for column in range(3)]
    loading_ax = fig.add_subplot(grid[1, :])
    for panel_ax, panel_letter in zip([*axes, loading_ax], "abcd"):
        panel_ax.set_title(panel_letter, loc="left", fontweight="bold")

    comparison_order = [
        "Rhizosphere-Deep",
        "Rhizosphere-Surface",
        "Deep-Surface",
    ]
    labels = [
        "Root-adjacent − shallow",
        "Root-adjacent − surface",
        "Shallow − surface",
    ]
    y = np.arange(len(comparison_order))

    primary_location = location[
        (location["analysis"] == "primary")
        & (location["contrast"].isin(comparison_order))
    ].copy()
    primary_location["contrast"] = pd.Categorical(
        primary_location["contrast"], comparison_order, ordered=True
    )
    primary_location = primary_location.sort_values("contrast")
    if primary_location["contrast"].astype(str).tolist() != comparison_order:
        raise ValueError("Missing a primary paired composition contrast")
    location_low = primary_location["standardized_displacement_ci_low"]
    location_high = primary_location["standardized_displacement_ci_high"]
    axes[0].errorbar(
        primary_location["standardized_displacement"],
        y,
        xerr=np.vstack(
            [
                primary_location["standardized_displacement"] - location_low,
                location_high - primary_location["standardized_displacement"],
            ]
        ),
        fmt="o",
        color=COLORS["root"],
        capsize=3,
    )
    axes[0].axvline(0, color="#777777", linewidth=0.8)
    axes[0].set(
        yticks=y,
        yticklabels=labels,
        xlabel="Consistency of composition shift\n(0 = cancelling, 1 = aligned)",
    )

    all_shannon = paired[
        (paired["trip"].astype(str) == "all")
        & (paired["metric"] == "shannon")
    ].copy()
    all_shannon["comparison"] = pd.Categorical(
        all_shannon["comparison"], comparison_order, ordered=True
    )
    all_shannon = all_shannon.sort_values("comparison")
    if all_shannon["comparison"].astype(str).tolist() != comparison_order:
        raise ValueError("Missing an all-campaign Shannon comparison")
    axes[1].errorbar(
        all_shannon["mean_difference"],
        y,
        xerr=asymmetric_errors(all_shannon, "mean_difference"),
        fmt="o",
        color="#333333",
        capsize=3,
    )
    axes[1].axvline(0, color="#777777", linewidth=0.8)
    axes[1].set(
        yticks=y,
        xlabel="Paired Shannon difference",
    )
    axes[1].tick_params(axis="y", labelleft=False)

    normalized = evenness[
        evenness["metric"] == "evenness_h_over_log_hurlbert"
    ].copy()
    normalized["contrast"] = pd.Categorical(
        normalized["contrast"], comparison_order, ordered=True
    )
    normalized = normalized.sort_values("contrast")
    if normalized["contrast"].astype(str).tolist() != comparison_order:
        raise ValueError("Missing a normalized-evenness contrast")
    errors = np.vstack(
        [
            normalized["mean_difference"].to_numpy()
            - normalized["bootstrap_ci_low"].to_numpy(),
            normalized["bootstrap_ci_high"].to_numpy()
            - normalized["mean_difference"].to_numpy(),
        ]
    )
    axes[2].errorbar(
        normalized["mean_difference"],
        y,
        xerr=errors,
        fmt="o",
        color=COLORS["root"],
        capsize=3,
    )
    axes[2].axvline(0, color="#777777", linewidth=0.8)
    axes[2].set(
        yticks=y,
        xlabel="Paired normalized\nShannon difference",
    )
    axes[2].tick_params(axis="y", labelleft=False)
    axes[2].set_xticks([-0.05, 0.0, 0.025], labels=["−0.05", "0", "0.025"])

    for axis in axes:
        axis.invert_yaxis()
        axis.tick_params(labelsize=11.5)
        axis.xaxis.label.set_size(11.5)

    # This panel is descriptive: it names the genera contributing most to the
    # three paired CLR displacement vectors without turning their loadings into
    # univariate differential-abundance tests.
    loadings = loadings[loadings["contrast"].isin(comparison_order)].copy()
    selected_genera: list[str] = []
    for comparison in comparison_order:
        ranked = (
            loadings[loadings["contrast"] == comparison]
            .assign(abs_loading=lambda frame: frame["mean_clr_difference"].abs())
            .nlargest(4, "abs_loading")
        )
        for genus in ranked["genus"]:
            if genus not in selected_genera:
                selected_genera.append(genus)
    if len(selected_genera) < 6:
        raise ValueError("Too few genera selected for the position-loading panel")
    loading_matrix = (
        loadings[loadings["genus"].isin(selected_genera)]
        .pivot(index="contrast", columns="genus", values="mean_clr_difference")
        .reindex(index=comparison_order, columns=selected_genera)
    )
    if loading_matrix.isna().any().any():
        raise ValueError("Position-loading panel has missing contrast values")
    maximum = max(2.5, float(np.abs(loading_matrix.to_numpy()).max()))
    image_plot = loading_ax.imshow(
        loading_matrix.to_numpy(),
        cmap="RdBu_r",
        vmin=-maximum,
        vmax=maximum,
        aspect="auto",
    )
    for row in range(loading_matrix.shape[0]):
        for column in range(loading_matrix.shape[1]):
            value = float(loading_matrix.iloc[row, column])
            loading_ax.text(
                column,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=12.0,
                color="white" if abs(value) > maximum * 0.58 else "#222222",
            )
    loading_ax.set(
        xticks=np.arange(len(selected_genera)),
        xticklabels=[rf"$\it{{{name}}}$" for name in selected_genera],
        yticks=np.arange(len(comparison_order)),
        yticklabels=labels,
    )
    loading_ax.tick_params(axis="x", rotation=28, labelsize=12.0)
    loading_ax.tick_params(axis="y", labelsize=12.0)
    colourbar = fig.colorbar(
        image_plot, ax=loading_ax, fraction=0.025, pad=0.02
    )
    colourbar.set_label("Relative-abundance difference\n(log ratio)")

    fig.tight_layout(h_pad=2.0, w_pad=1.6)
    fig.savefig(output, bbox_inches="tight", metadata=PDF_METADATA)
    plt.close(fig)


def make_function_control_figure(
    position: pd.DataFrame,
    pma_pairs: pd.DataFrame,
    pma_summary: dict[str, Any],
    removal: pd.DataFrame,
    output: Path,
) -> None:
    """Show predicted pathway structure, PMA outcomes and control removal."""
    fig, (pathway_ax, pma_ax, removal_ax) = plt.subplots(
        1, 3, figsize=(10.6, 3.4)
    )
    for panel_ax, panel_letter in zip((pathway_ax, pma_ax, removal_ax), "abc"):
        panel_ax.set_title(panel_letter, loc="left", fontweight="bold")

    primary = position[
        position["analysis"].eq("primary")
        & ~position["contrast"].eq("omnibus_three_positions")
    ].copy()
    contrast_order = [
        "Rhizosphere-Deep",
        "Rhizosphere-Surface",
        "Deep-Surface",
    ]
    contrast_labels = [
        "Root-adjacent − shallow",
        "Root-adjacent − surface",
        "Shallow − surface",
    ]
    primary["contrast"] = pd.Categorical(
        primary["contrast"], contrast_order, ordered=True
    )
    primary = primary.sort_values("contrast")
    if primary["contrast"].astype(str).tolist() != contrast_order:
        raise ValueError("Missing a primary PICRUSt2 position contrast")
    y = np.arange(3)
    pathway_low = primary["standardized_ci_low"]
    pathway_high = primary["standardized_ci_high"]
    pathway_ax.errorbar(
        primary["standardized_displacement"],
        y,
        xerr=np.vstack(
            [
                primary["standardized_displacement"] - pathway_low,
                pathway_high - primary["standardized_displacement"],
            ]
        ),
        fmt="o",
        color=COLORS["root"],
        capsize=3,
    )
    pathway_ax.axvline(0, color="#777777", linewidth=0.8)
    pathway_ax.set(
        yticks=y,
        yticklabels=contrast_labels,
        xlabel="Consistency of pathway shift\n(0 = cancelling, 1 = aligned)",
    )
    pathway_ax.invert_yaxis()

    pma_groups = pma_pairs.assign(group=pma_pairs["pair_id"].str.extract(r"^(C[12][RS])", expand=False))
    expected_groups = {"C1R", "C2R", "C2S"}
    if set(pma_groups["group"]) != expected_groups or len(pma_groups) != 9:
        raise ValueError("Expected three PMA aliquot groups at two campsites")
    for group, color in zip(("C1R", "C2R", "C2S"), ("#0072B2", "#009E73", "#D55E00")):
        values = pma_groups[pma_groups["group"].eq(group)]
        if len(values) != 3:
            raise ValueError(f"Expected three aliquot comparisons for {group}")
        for row in values.itertuples():
            pma_ax.plot(
                [0, 1],
                [row.untreated_expected_rarefied_richness, row.treated_expected_rarefied_richness],
                color=color, linewidth=0.9, alpha=0.5, marker="o", markersize=3.4,
            )
        pma_ax.plot(
            [0, 1],
            [values["untreated_expected_rarefied_richness"].mean(), values["treated_expected_rarefied_richness"].mean()],
            color=color, linewidth=2.1, marker="D", markersize=4,
            label=f"{group} mean", zorder=4,
        )
    pma_ax.set(
        xticks=[0, 1],
        xticklabels=["Untreated", "PMA treated"],
        ylabel="Expected richness",
    )
    pma_ax.legend(frameon=False, fontsize=8.0)

    biological = removal[
        removal["role"].eq("compatible_biological_profile")
    ].copy()
    fractions = np.sort(
        biological["candidate_contaminant_read_fraction"].to_numpy(dtype=float)
        * 100
    )
    removal_ax.scatter(
        np.arange(1, len(fractions) + 1),
        np.maximum(fractions, 0.0001),
        s=10,
        color="#009E73",
        alpha=0.75,
        edgecolor="none",
    )
    removal_ax.axhline(
        np.median(fractions), color="#CC3311", linewidth=1.2, linestyle="--"
    )
    fraction_q1, fraction_q3 = np.quantile(fractions, [0.25, 0.75])
    removal_ax.axhspan(
        fraction_q1,
        fraction_q3,
        color="#CC3311",
        alpha=0.10,
        linewidth=0,
    )
    removal_ax.set_yscale("log")
    removal_ax.set(
        xlabel="Trip 5 profiles,\nordered by fraction",
        ylabel="Contaminant reads removed (%)",
    )
    removal_ax.text(
        0.03,
        0.95,
        f"median {np.median(fractions):.2f}%\n"
        f"middle 50% {fraction_q1:.2f}–{fraction_q3:.2f}%\n"
        f"maximum {np.max(fractions):.1f}%",
        transform=removal_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
    )

    fig.tight_layout(w_pad=2.2)
    fig.savefig(output, bbox_inches="tight", metadata=PDF_METADATA)
    plt.close(fig)


def main() -> None:
    figure_runtime = require_figure_runtime()
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-dir", type=Path, required=True)
    # Retained as optional compatibility arguments for older workflow calls.
    parser.add_argument("--network-dir", type=Path, default=None)
    parser.add_argument("--functional-dir", type=Path, default=None)
    parser.add_argument("--environment-dir", type=Path, default=None)
    parser.add_argument("--picrust-dir", type=Path, default=None)
    parser.add_argument("--rain-dir", type=Path, default=None)
    parser.add_argument("--ordination-dir", type=Path, default=None)
    parser.add_argument("--control-dir", type=Path, default=None)
    parser.add_argument("--ph-dir", type=Path, default=None)
    parser.add_argument("--xrf-dir", type=Path, default=None)
    parser.add_argument("--biology-dir", type=Path, default=None)
    parser.add_argument("--pma-dir", type=Path, default=None)
    parser.add_argument("--measured-function-dir", type=Path, default=None)
    parser.add_argument("--boundary-file", type=Path, default=None)
    parser.add_argument("--background-image", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    core = args.core_dir.resolve()
    local_v3 = Path(__file__).resolve().parent
    environment_dir = (
        args.environment_dir.resolve()
        if args.environment_dir is not None
        else local_v3 / "environment_associations"
    )
    picrust_dir = (
        args.picrust_dir.resolve()
        if args.picrust_dir is not None
        else local_v3 / "picrust2_ecology"
    )
    control_dir = (
        args.control_dir.resolve()
        if args.control_dir is not None
        else local_v3 / "control_audit"
    )
    ph_dir = (
        args.ph_dir.resolve()
        if args.ph_dir is not None
        else local_v3 / "ph_group_linkage_20260909"
    )
    xrf_dir = (
        args.xrf_dir.resolve()
        if args.xrf_dir is not None
        else local_v3 / "xrf_community_rescue"
    )
    biology_dir = (
        args.biology_dir.resolve()
        if args.biology_dir is not None
        else local_v3 / "biology_context_corrected_20260909"
    )
    pma_dir = (
        args.pma_dir.resolve()
        if args.pma_dir is not None
        else local_v3 / "pma_endpoint_results"
    )
    measured_function_dir = (
        args.measured_function_dir.resolve()
        if args.measured_function_dir is not None
        else local_v3 / "measured_function_summary_results"
    )
    project_root = Path(__file__).resolve().parents[2]
    if args.boundary_file is not None:
        boundary_file = args.boundary_file.resolve()
    else:
        boundary_candidates = (
            project_root / "data/metadata/misc/boundary.kml",
            project_root / "metadata/geodata/empty_quarter_boundary.kml",
        )
        boundary_file = next(
            (candidate for candidate in boundary_candidates if candidate.is_file()),
            boundary_candidates[0],
        )
    background_file = (
        args.background_image.resolve()
        if args.background_image is not None
        else project_root / "metadata/geodata/bluemarble_arabia_200407_120ppd.png"
    )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    rain_response_dir = (
        args.rain_dir.resolve()
        if args.rain_dir is not None
        else local_v3 / "rain_pulse_response"
    )
    ordination_dir = (
        args.ordination_dir.resolve()
        if args.ordination_dir is not None
        else local_v3 / "aitchison_ordination"
    )
    distance_decay_dir = core / "distance_decay_turnover"
    if not (distance_decay_dir / "distance_decay_pairs.tsv").is_file():
        distance_decay_dir = core.parent / "distance_decay_turnover"
    if not (distance_decay_dir / "distance_decay_pairs.tsv").is_file():
        distance_decay_dir = local_v3 / "distance_decay_turnover"

    input_paths = {
        "alpha_table": core / "cache/alpha.tsv",
        "empty_quarter_orientation_boundary": boundary_file,
        "landscape_background_image": background_file,
        "landscape_background_sidecar": background_file.with_suffix(".json"),
        "spatial_site_coordinates": (
            core
            / "spatial_turnover_rescue/results/site_coordinates.tsv"
        ),
        "pma_summary": core / "pma_endpoints/pma_summary.json",
        "paired_compartment_effects": (
            core / "claim_rescue/paired_compartment_effects.tsv"
        ),
        "paired_evenness_effects": (
            core / "evenness_decomposition/paired_contrasts.tsv"
        ),
        "paired_composition_location": (
            core
            / "compartment_composition/compartment_location_results.tsv"
        ),
        "paired_composition_loadings": (
            core
            / "compartment_composition/paired_displacement_loadings.tsv"
        ),
        "rain_response_figure": rain_response_dir / "rain_pulse_response.pdf",
        "rain_response_decision": rain_response_dir / "analysis_decision.json",
        "spatial_claim_verdict": (
            core / "spatial_turnover_rescue/results/claim_verdict.json"
        ),
        "distance_decay_pairs": (
            distance_decay_dir / "distance_decay_pairs.tsv"
        ),
        "ordination_scores": ordination_dir / "ordination_scores.tsv",
        "ordination_summary": ordination_dir / "ordination_summary.json",
        "climate_site_summary": (
            environment_dir / "climate_site_summary.tsv"
        ),
        "climate_alpha_correlations": (
            environment_dir / "climate_alpha_correlations.tsv"
        ),
        "climate_genus_correlations": (
            environment_dir / "climate_genus_correlations.tsv"
        ),
        "environment_decision": (
            environment_dir / "analysis_decision.json"
        ),
        "picrust_position_tests": (
            picrust_dir / "position_profile_tests.tsv"
        ),
        "picrust_decision": picrust_dir / "analysis_decision.json",
        "ko_validation": (
            measured_function_dir / "per_sample_ko_correlations.tsv"
        ),
        "ko_metrics": measured_function_dir / "summary_metrics.tsv",
        "pma_pairs": pma_dir / "pma_pair_endpoints.tsv",
        "control_removal": (
            control_dir / "trip5_removal_fraction_by_profile.tsv"
        ),
        "control_sensitivity_summary": (
            control_dir / "sensitivity_inputs/summary.json"
        ),
        "ph_group_table": ph_dir / "ph_group_analysis_table.tsv",
        "laboratory_xrf_axis": xrf_dir / "laboratory_xrf_axis.tsv",
        "site_landforms": biology_dir / "site_landforms.tsv",
    }
    missing = [str(path) for path in input_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing canonical figure inputs:\n" + "\n".join(missing))

    alpha = pd.read_csv(input_paths["alpha_table"], sep="\t", index_col=0)
    coordinates = pd.read_csv(
        input_paths["spatial_site_coordinates"], sep="\t"
    )
    paired = pd.read_csv(input_paths["paired_compartment_effects"], sep="\t")
    evenness = pd.read_csv(input_paths["paired_evenness_effects"], sep="\t")
    location = pd.read_csv(
        input_paths["paired_composition_location"], sep="\t"
    )
    loadings = pd.read_csv(
        input_paths["paired_composition_loadings"], sep="\t"
    )
    distance_pairs = pd.read_csv(
        input_paths["distance_decay_pairs"], sep="\t"
    )
    ph_groups = pd.read_csv(input_paths["ph_group_table"], sep="\t")
    xrf_axis = pd.read_csv(input_paths["laboratory_xrf_axis"], sep="\t")
    landforms = pd.read_csv(input_paths["site_landforms"], sep="\t")
    require_columns(
        ph_groups,
        {"trip", "site", "compartment", "ph"},
        input_paths["ph_group_table"],
    )
    require_columns(
        xrf_axis,
        {"Trip", "Site", "Type", "elemental_pc1"},
        input_paths["laboratory_xrf_axis"],
    )
    require_columns(
        landforms, {"site", "landform"}, input_paths["site_landforms"]
    )
    ordination = pd.read_csv(input_paths["ordination_scores"], sep="\t")
    ordination_summary = json.loads(
        input_paths["ordination_summary"].read_text(encoding="utf-8")
    )
    require_columns(
        ordination,
        {"campaign", "site", "compartment", "transect_km", "pc1", "pc2"},
        input_paths["ordination_scores"],
    )
    if (
        len(ordination) != ordination_summary["n_groups"]
        or ordination_summary["n_genera"] != 200
        or ordination[["pc1", "pc2", "transect_km"]].isna().any().any()
    ):
        raise ValueError("Ordination scores do not match their summary")
    climate_site = pd.read_csv(
        input_paths["climate_site_summary"], sep="\t"
    )
    climate_alpha = pd.read_csv(
        input_paths["climate_alpha_correlations"], sep="\t"
    )
    climate_genus = pd.read_csv(
        input_paths["climate_genus_correlations"], sep="\t"
    )
    picrust_position = pd.read_csv(
        input_paths["picrust_position_tests"], sep="\t"
    )
    ko_validation = pd.read_csv(input_paths["ko_validation"], sep="\t")
    ko_metrics = pd.read_csv(input_paths["ko_metrics"], sep="\t")
    pma_pairs = pd.read_csv(input_paths["pma_pairs"], sep="\t")
    control_removal = pd.read_csv(
        input_paths["control_removal"], sep="\t"
    )
    spatial_verdict = json.loads(input_paths["spatial_claim_verdict"].read_text())
    pma_summary = json.loads(input_paths["pma_summary"].read_text())
    rain_response_verdict = json.loads(
        input_paths["rain_response_decision"].read_text()
    )
    environment_verdict = json.loads(
        input_paths["environment_decision"].read_text()
    )
    picrust_verdict = json.loads(input_paths["picrust_decision"].read_text())
    control_sensitivity = json.loads(
        input_paths["control_sensitivity_summary"].read_text()
    )

    require_columns(
        alpha,
        {"Trip", "Site", "Type"},
        input_paths["alpha_table"],
    )
    require_columns(
        coordinates,
        {"site", "latitude", "longitude", "transect_km"},
        input_paths["spatial_site_coordinates"],
    )
    require_columns(
        paired,
        {
            "trip",
            "metric",
            "comparison",
            "mean_difference",
            "ci_low",
            "ci_high",
        },
        input_paths["paired_compartment_effects"],
    )
    require_columns(
        evenness,
        {
            "metric",
            "contrast",
            "mean_difference",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
        },
        input_paths["paired_evenness_effects"],
    )
    require_columns(
        location,
        {
            "analysis",
            "contrast",
            "displacement",
            "displacement_ci_low",
            "displacement_ci_high",
            "standardized_displacement",
            "standardized_displacement_ci_low",
            "standardized_displacement_ci_high",
            "permutation_p",
        },
        input_paths["paired_composition_location"],
    )
    require_columns(
        loadings,
        {"contrast", "genus", "mean_clr_difference"},
        input_paths["paired_composition_loadings"],
    )
    require_columns(
        distance_pairs,
        {
            "site_a",
            "site_b",
            "geographic_distance_km",
            "compartment",
            "compartment_label",
            "aitchison_dissimilarity",
        },
        input_paths["distance_decay_pairs"],
    )
    require_columns(
        climate_site,
        {
            "site",
            "transect_km",
            "mean_air_temperature_c",
            "mean_monthly_rain_mm",
            "mean_relative_humidity_pct",
        },
        input_paths["climate_site_summary"],
    )
    require_columns(
        climate_alpha,
        {
            "climate_variable",
            "response",
            "spearman_rho",
            "q_global_9",
        },
        input_paths["climate_alpha_correlations"],
    )
    require_columns(
        climate_genus,
        {
            "climate_variable",
            "genus",
            "spearman_rho",
            "q_global_600",
            "supported_q_lt_0_05",
        },
        input_paths["climate_genus_correlations"],
    )
    require_columns(
        picrust_position,
        {
            "analysis",
            "contrast",
            "displacement",
            "ci_low",
            "ci_high",
            "standardized_displacement",
            "standardized_ci_low",
            "standardized_ci_high",
            "permutation_p",
        },
        input_paths["picrust_position_tests"],
    )
    require_columns(
        ko_validation,
        {"sample", "n_shared_kos", "spearman_rho"},
        input_paths["ko_validation"],
    )
    require_columns(
        pma_pairs,
        {
            "pair_id",
            "treated_expected_rarefied_richness",
            "untreated_expected_rarefied_richness",
        },
        input_paths["pma_pairs"],
    )
    require_columns(
        control_removal,
        {"profile_id", "role", "candidate_contaminant_read_fraction"},
        input_paths["control_removal"],
    )

    # Fail closed if a verdict changes: the plot captions and manuscript text
    # would then require scientific review rather than automatic reuse.
    expected_verdicts = {
        "spatial": (
            spatial_verdict["status"],
            "broad_geographic_structure_supported",
        ),
        "pma": (pma_summary["status"], "paired_endpoints_only"),
        "rain_response": (
            rain_response_verdict["analysis_status"],
            "temporally_localized_association_borderline",
        ),
        "environment": (
            environment_verdict["status"],
            "observational_climate_associations_supported",
        ),
        "picrust2": (
            picrust_verdict["status"],
            "predicted_functional_structure_supported",
        ),
        "control_sensitivity": (
            control_sensitivity["profiles_below_rarefaction_depth_after_filter"],
            0,
        ),
    }
    changed = {
        name: observed
        for name, (observed, expected) in expected_verdicts.items()
        if observed != expected
    }
    if changed:
        raise ValueError(f"Claim verdict changed; review figures first: {changed}")

    setup_style()
    output_paths = {
        "landscape": output / "fig1_landscape.pdf",
        "composition_geography": output / "fig2_composition_geography.pdf",
        "soil_position": output / "fig3_soil_position.pdf",
        "campaign_rainfall_supplement": (
            output / "figS_campaign_rainfall.pdf"
        ),
        "function_controls": output / "fig5_function_controls.pdf",
        "environment_gradients": output / "fig4_environment_gradients.pdf",
    }
    make_landscape_figure(
        alpha,
        coordinates,
        boundary_file,
        background_file,
        climate_site,
        output_paths["landscape"],
    )
    make_composition_geography_figure(
        ordination,
        ordination_summary,
        distance_pairs,
        output_paths["composition_geography"],
    )
    make_soil_position_figure(
        paired,
        evenness,
        location,
        loadings,
        output_paths["soil_position"],
    )
    shutil.copyfile(
        input_paths["rain_response_figure"],
        output_paths["campaign_rainfall_supplement"],
    )
    make_function_control_figure(
        picrust_position,
        pma_pairs,
        pma_summary,
        control_removal,
        output_paths["function_controls"],
    )
    make_environment_gradient_figure(
        ph_groups,
        xrf_axis,
        coordinates,
        climate_site,
        landforms,
        climate_alpha,
        climate_genus,
        output_paths["environment_gradients"],
    )

    rows = []
    for name, path in sorted(input_paths.items()):
        rows.append(
            {
                "role": "input",
                "name": name,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    for name, path in sorted(output_paths.items()):
        rows.append(
            {
                "role": "output",
                "name": name,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(
        output / "figure_manifest.tsv", sep="\t", index=False
    )
    (output / "figure_runtime.json").write_text(
        json.dumps(
            {
                "schema_version": "figure-runtime-v1",
                **figure_runtime,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
