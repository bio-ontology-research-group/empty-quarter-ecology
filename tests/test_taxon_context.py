"""Descriptive taxon-context module: outputs, checksums and manuscript numbers."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis/v3/taxon_context"
PAPER = ROOT / "empty-quarter-amplicon"
MAIN = (PAPER / "main.tex").read_text(encoding="utf-8")
SUPPLEMENT = (PAPER / "supplement.tex").read_text(encoding="utf-8")
FRAGMENT = PAPER / "generated/taxon_context_tables.tex"


def flat(text: str) -> str:
    return " ".join(text.split()).replace("$", "").replace("{,}", ",")


def test_checksums_are_current() -> None:
    for line in (RESULTS / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert hashlib.sha256((RESULTS / name).read_bytes()).hexdigest() == expected, name


def test_manifest_scope_and_cohort() -> None:
    manifest = json.loads((RESULTS / "run_manifest.json").read_text())
    assert manifest["status"] == "descriptive_context_complete"
    assert manifest["cohort"]["core_site_profiles"] == 1227
    assert manifest["cohort"]["genus_sum_max_abs_difference_vs_cache"] == 0
    assert manifest["cohort"]["primary_genera"] == 200
    assert manifest["parameters"]["rarefaction_depth"] == 8000
    assert manifest["parameters"]["detection_prevalence"] == 0.10
    assert "no feature-table merging" in manifest["scope"]
    for third in ("west", "central", "east"):
        assert len(manifest["cohort"]["transect_thirds"][third]) == 20


def test_phylum_and_genus_numbers_match_main_text() -> None:
    phyla = pd.read_csv(RESULTS / "phylum_composition.tsv", sep="\t").set_index("phylum")
    main = flat(MAIN)
    for phylum, printed in (
        ("Pseudomonadota", "26.9"),
        ("Bacillota", "22.5"),
        ("Actinomycetota", "18.7"),
        ("Chloroflexota", "8.6"),
        ("Bacteroidota", "7.1"),
        ("Planctomycetota", "3.2"),
        ("Gemmatimonadota", "2.7"),
        ("Acidobacteriota", "2.3"),
        ("Cyanobacteriota", "0.4"),
    ):
        assert f"{100 * phyla.loc[phylum, 'mean_relative_abundance']:.1f}" == printed
        assert f"{phylum} ({printed}\\,\\%)" in main or f"{phylum} held {printed}\\,\\%" in main
    combined = phyla.loc[["Acidobacteriota", "Verrucomicrobiota"], "mean_relative_abundance"].sum()
    assert f"{100 * combined:.1f}" == "3.3"
    assert "together held 3.3\\,\\%" in main

    genera = pd.read_csv(RESULTS / "genus_composition.tsv", sep="\t").set_index("genus")
    for genus, printed in (
        ("Domibacillus", "6.2"),
        ("Massilia", "3.5"),
        ("Bacillus", "3.0"),
        ("Microvirga", "2.4"),
        ("Flavisolibacter", "2.0"),
    ):
        assert f"{100 * genera.loc[genus, 'mean_relative_abundance']:.1f}" == printed
        assert f"\\textit{{{genus}}} ({printed}\\,\\%" in main
    leaders = pd.read_csv(RESULTS / "leading_genera_by_stratum.tsv", sep="\t")
    for stratum in ("surface", "shallow_subsurface", "root_adjacent"):
        top = leaders[(leaders["stratum"] == stratum) & (leaders["rank_in_stratum"] == 1)]
        assert top["genus"].iloc[0] == "Domibacillus"
    east = leaders[(leaders["stratum"] == "east_third") & (leaders["rank_in_stratum"] == 1)]
    assert east["genus"].iloc[0] == "Bacillus"

    classes = pd.read_csv(RESULTS / "class_composition.tsv", sep="\t").set_index("class")
    for name, printed in (
        ("Bacilli", "21.0"),
        ("Gammaproteobacteria", "14.7"),
        ("Actinobacteria", "13.3"),
        ("Alphaproteobacteria", "12.2"),
    ):
        assert f"{100 * classes.loc[name, 'mean_relative_abundance']:.1f}" == printed
        assert printed in main


def test_replacement_numbers_match_text() -> None:
    table = pd.read_csv(RESULTS / "transect_replacement.tsv", sep="\t")
    supported = table[table["supported_q_lt_0_05"]]
    assert len(supported) == 124
    assert (supported["direction"] == "decreases_eastward").sum() == 77
    assert (supported["direction"] == "increases_eastward").sum() == 47
    by_phylum = supported.groupby(["direction", "phylum"]).size()
    assert by_phylum[("decreases_eastward", "Actinomycetota")] == 27
    assert by_phylum[("decreases_eastward", "Pseudomonadota")] == 20
    assert by_phylum[("decreases_eastward", "Bacteroidota")] == 9
    assert by_phylum[("decreases_eastward", "Acidobacteriota")] == 4
    assert by_phylum[("increases_eastward", "Bacillota")] == 17
    assert by_phylum[("increases_eastward", "Pseudomonadota")] == 16
    rows = table.set_index("genus")
    assert f"{rows.loc['Cellulomonas', 'spearman_rho_route_position']:.2f}" == "-0.79"
    assert f"{rows.loc['Bacillus', 'spearman_rho_route_position']:.2f}" == "0.77"
    assert f"{rows.loc['Halalkalibacter', 'spearman_rho_route_position']:.2f}" == "0.85"
    assert f"{100 * rows.loc['Bacillus', 'mean_relative_abundance_west_third']:.1f}" == "1.6"
    assert f"{100 * rows.loc['Bacillus', 'mean_relative_abundance_east_third']:.1f}" == "6.2"
    main = flat(MAIN)
    assert "124 changed monotonically along the route" in main
    assert "77 declined from west to east and 47 increased" in main
    assert "\\textit{Cellulomonas}, Spearman \\rho=-0.79" in main
    assert "rose from 1.6 to 6.2\\,\\% of reads" in main
    assert "\\textit{Halalkalibacter} (\\rho=0.85)" in main

    gradients = pd.read_csv(RESULTS / "site_gradients.tsv", sep="\t").set_index("variable")
    for variable, printed in (
        ("mean_air_temperature_c", "0.98"),
        ("mean_relative_humidity_pct", "0.92"),
        ("mean_ph", "0.82"),
    ):
        assert f"{gradients.loc[variable, 'spearman_rho_route_position']:.2f}" == printed
    assert "0.98 for temperature, 0.92 for humidity and 0.82 for pH" in main


def test_overlap_numbers_match_text() -> None:
    overlap = pd.read_csv(RESULTS / "genus_set_overlap.tsv", sep="\t")
    west_east = overlap.iloc[0]
    assert west_east["set_a"].endswith("west third") and west_east["set_b"].endswith("east third")
    assert int(west_east["n_shared"]) == 267
    assert int(west_east["n_genera_a"] + west_east["n_genera_b"] - west_east["n_shared"]) == 418
    assert f"{west_east['jaccard']:.2f}" == "0.64"
    assert int(west_east["top50_a_detected_in_b"]) == 50
    assert int(west_east["top50_b_detected_in_a"]) == 44
    pit = overlap.iloc[3]
    assert "Atacama pit all depths" in pit["set_b"]
    assert int(pit["n_shared"]) == 43
    assert int(pit["n_genera_a"] + pit["n_genera_b"] - pit["n_shared"]) == 397
    assert f"{pit['jaccard']:.2f}" == "0.11"
    assert int(pit["top50_a_detected_in_b"]) == 17
    assert int(pit["top50_b_detected_in_a"]) == 31
    assert round(west_east["jaccard"] / pit["jaccard"]) == 6
    main = flat(MAIN)
    assert "shared 267 of 418 genera (Jaccard 0.64; 0.68 and 0.83 for adjacent thirds)" in main
    assert "shared 43 of 397 genera (Jaccard 0.11)" in main
    assert "only 17 of the 50 leading Empty Quarter genera" in main
    supplement = flat(SUPPLEMENT)
    assert "Jaccard 0.639" in supplement and "Jaccard 0.108" in supplement
    assert "31 of the 50 leading pit genera" in supplement


def test_pathway_numbers_match_text() -> None:
    classes = pd.read_csv(RESULTS / "pathway_class_share.tsv", sep="\t").set_index("pathway_class")
    for name, printed in (
        ("biosynthesis", "66.8"),
        ("energy_central_metabolism", "17.2"),
        ("degradation_utilization", "10.7"),
    ):
        assert f"{100 * classes.loc[name, 'share_of_predicted_pathway_abundance']:.1f}" == printed
    dominance = pd.read_csv(RESULTS / "pathway_dominance.tsv", sep="\t").set_index("pathway")
    assert dominance.loc["PWY-3781", "rank"] == 1
    assert f"{100 * dominance.loc['PWY-3781', 'mean_relative_abundance']:.1f}" == "1.7"
    assert dominance.loc["CALVIN-PWY", "rank"] == 44
    assert f"{100 * dominance.loc['CALVIN-PWY', 'mean_relative_abundance']:.1f}" == "0.6"
    assert dominance.loc["P23-PWY", "rank"] == 122
    assert "PWY-101" not in dominance.index
    supported = pd.read_csv(RESULTS / "supported_pathways_by_class.tsv", sep="\t")
    route = supported[supported["family"] == "route_correlation_supported"]
    assert int(route["n_supported"].sum()) == 92
    biosynthesis = route[route["pathway_class"] == "biosynthesis"].iloc[0]
    assert int(biosynthesis["n_supported"]) == 63 and int(biosynthesis["n_positive"]) == 54
    compartment = supported[supported["family"] == "compartment_contrast_supported"]
    assert int(compartment["n_supported"].sum()) == 270
    main = flat(MAIN)
    assert "Biosynthesis pathways held 66.8\\,\\% of predicted pathway abundance" in main
    assert "Calvin--Benson--Bassham cycle ranked 44th (0.6\\,\\%)" in main
    assert "reductive TCA cycle 122nd (0.4\\,\\%)" in main
    assert "54 of 63" in main


def test_generated_tables_are_current_and_included() -> None:
    assert "\\input{generated/taxon_context_tables.tex}" in SUPPLEMENT
    rendered = subprocess.run(
        [
            sys.executable,
            str(ROOT / "analysis/v3/render_taxon_context_tex.py"),
            "--output",
            "/dev/stdout",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    # The script prints the output path after writing; /dev/stdout receives
    # the fragment followed by that path line.
    fragment = FRAGMENT.read_text(encoding="utf-8")
    assert fragment.strip() in rendered
    for label in (
        "tab:taxa-phyla",
        "tab:taxa-genera",
        "tab:taxa-replacement",
        "tab:route-gradients",
        "tab:genus-overlap",
        "tab:pathway-classes",
        "tab:pathway-dominance",
    ):
        assert f"\\label{{{label}}}" in fragment
        assert f"\\ref{{{label}}}" in SUPPLEMENT


def test_landscape_figure_uses_the_committed_satellite_crop() -> None:
    sidecar = json.loads(
        (ROOT / "metadata/geodata/bluemarble_arabia_200407_120ppd.json").read_text()
    )
    crop = ROOT / "metadata/geodata/bluemarble_arabia_200407_120ppd.png"
    assert hashlib.sha256(crop.read_bytes()).hexdigest() == sidecar["output_sha256"]
    assert sidecar["extent_degrees"] == {"lon_min": 44.0, "lon_max": 57.0, "lat_min": 16.0, "lat_max": 25.0}
    assert sidecar["pixels_per_degree"] == 120
    manifest = pd.read_csv(PAPER / "figures/figure_manifest.tsv", sep="\t")
    row = manifest[manifest["name"] == "landscape_background_image"].iloc[0]
    assert row["sha256"] == sidecar["output_sha256"]
    assert "NASA Blue Marble" in MAIN
