"""Descriptive taxon-context module: outputs, checksums and manuscript numbers."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from manuscript_paths import PAPER

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis/v3/taxon_context_corrected_20260909"
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
    composition = (PAPER / "generated/taxon_composition_tables.tex").read_text()
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
        if phylum in ("Pseudomonadota", "Bacillota", "Actinomycetota"):
            assert f"{phylum} ({printed}\\,\\%)" in main
        elif phylum == "Cyanobacteriota":
            # The abundance is retained as context in the SI; the revised
            # main text no longer uses it to classify the community's metabolism.
            assert f"Cyanobacteriota contributed {printed}\\,\\%" in flat(SUPPLEMENT)
        else:
            row = next(line for line in composition.splitlines() if line.startswith(phylum + " &"))
            assert row.split(" & ")[2] == printed
    combined = phyla.loc[["Acidobacteriota", "Verrucomicrobiota"], "mean_relative_abundance"].sum()
    assert f"{100 * combined:.1f}" == "3.3"
    # The combined 3.3% remains an artifact check; the shortened main text
    # reports three leading phyla/genera and relocates classes to the SI.

    genera = pd.read_csv(RESULTS / "genus_composition.tsv", sep="\t").set_index("genus")
    for genus, printed in (
        ("Domibacillus", "6.2"),
        ("Massilia", "3.5"),
        ("Bacillus", "3.0"),
        ("Microvirga", "2.4"),
        ("Flavisolibacter", "2.0"),
    ):
        assert f"{100 * genera.loc[genus, 'mean_relative_abundance']:.1f}" == printed
        if genus in ("Domibacillus", "Massilia", "Bacillus"):
            assert f"\\textit{{{genus}}} ({printed}\\,\\%" in main
        else:
            row = next(line for line in composition.splitlines() if line.startswith(f"\\textit{{{genus}}} &"))
            assert float(row.split(" & ")[2]) == round(100 * genera.loc[genus, "mean_relative_abundance"], 2)
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
        assert f"{name} ({printed}\\,\\%)" in flat(SUPPLEMENT)


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
    assert "124 changed significantly after correcting for multiple testing" in main
    assert "77 that declined along the transect" in main
    assert "47 that increased" in main
    assert "increased from 1.6 to 6.2\\,\\%" in main
    assert "\\textit{Halalkalibacter} (\\rho=0.85)" in main

    gradients = pd.read_csv(RESULTS / "site_gradients.tsv", sep="\t").set_index("variable")
    for variable, printed in (
        ("mean_air_temperature_c", "0.98"),
        ("mean_relative_humidity_pct", "0.92"),
        ("mean_ph", "0.82"),
    ):
        assert f"{gradients.loc[variable, 'spearman_rho_route_position']:.2f}" == printed
    assert "0.98 temperature, 0.92 humidity" in main
    # The pH gradient uses independently rejoined assay groups and is displayed
    # in the SI table after shortening the main-text environmental summary.
    assert "archived-soil pH (equal-weight assay-group means within site) & 60 & 0.82 & 7.85 & 8.23" in (PAPER / "generated/taxon_composition_tables.tex").read_text()


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
    supplement = flat(SUPPLEMENT)
    assert "shared 267 of 418 detected genera" in supplement
    assert "shared 43 of 397 genera" in supplement
    assert "17 of the 50 leading Empty Quarter genera" in supplement
    assert "Jaccard 0.639" in supplement and "Jaccard 0.108" in supplement
    assert "31 of the 50 leading pit genera" in supplement
    overlap_rows = (PAPER / "generated/taxon_overlap_tables.tex").read_text()
    west_east_row = next(line for line in overlap_rows.splitlines() if "west third" in line and "east third" in line)
    assert west_east_row.endswith(" & 50 " + r"\\")
    # The table prints A→B only; the revised prose must also preserve B→A=44.
    overlap_section = supplement.split("Genus-set overlap with the Atacama pit", 1)[1].split("Sampling designs", 1)[0]
    assert "All 50 leading western genera were detected in the east" in overlap_section
    assert "44 of the 50 leading eastern genera were detected in the west" in overlap_section


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
    assert "Biosynthesis pathways accounted for 66.8\\,\\% of predicted pathway abundance" in main
    assert "Calvin--Benson--Bassham cycle ranked 44th (0.6\\,\\%)" in main
    assert "reductive TCA cycle 122nd (0.4\\,\\%)" in main
    assert "63 were biosynthesis pathways (54 increasing eastward)" in flat(SUPPLEMENT)


def test_generated_tables_are_current_and_included() -> None:
    names = ("taxon_composition_tables.tex", "taxon_overlap_tables.tex", "taxon_pathway_tables.tex")
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp)
        subprocess.run([sys.executable, str(ROOT / "analysis/v3/render_taxon_context_tex.py"),
                        "--results", str(RESULTS), "--output", str(target / "all.tex"),
                        "--split-dir", str(target)], check=True, capture_output=True, text=True)
        assert (target / "all.tex").read_text() == FRAGMENT.read_text()
        for name in names:
            assert f"\\input{{generated/{name}}}" in SUPPLEMENT
            assert (target / name).read_bytes() == (PAPER / "generated" / name).read_bytes()
    fragment = "\n".join((PAPER / "generated" / name).read_text() for name in names)
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
    manifest = json.loads((PAPER / "figures/figure_review_manifest.json").read_text())
    row = manifest["inputs"]["background"]
    assert row["sha256"] == sidecar["output_sha256"]
    assert "NASA Blue Marble" in MAIN
