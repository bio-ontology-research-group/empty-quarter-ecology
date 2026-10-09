"""Descriptive taxon-context module: outputs, checksums and manuscript numbers."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from manuscript_paths import PAPER
import manuscript_text as mt
from manuscript_text import fmt, pct

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
    expected_phyla = {
        "Pseudomonadota": "26.9", "Bacillota": "22.5", "Actinomycetota": "18.7",
        "Chloroflexota": "8.6", "Bacteroidota": "7.1", "Planctomycetota": "3.2",
        "Gemmatimonadota": "2.7", "Acidobacteriota": "2.3", "Cyanobacteriota": "0.4",
    }
    for phylum, printed in expected_phyla.items():
        assert pct(phyla.loc[phylum, "mean_relative_abundance"]) == printed
    main = mt.text("main")
    for phylum in ("Pseudomonadota", "Bacillota", "Actinomycetota"):
        assert f"{phylum} was the most abundant phylum ({expected_phyla[phylum]}%)" in main or (
            f"{phylum} ({expected_phyla[phylum]}%)" in main
        ), phylum
    combined = phyla.loc[["Acidobacteriota", "Verrucomicrobiota"], "mean_relative_abundance"].sum()
    assert pct(combined) == "3.3"
    assert f"together accounted for {pct(combined)}%" in main
    supplement = mt.text("supplement")
    bacillota_actino = phyla.loc[["Bacillota", "Actinomycetota"], "mean_relative_abundance"].sum()
    assert f"together accounted for {pct(bacillota_actino)}% of reads" in supplement
    cyano = phyla.loc["Cyanobacteriota"]
    assert f"Cyanobacteriota contributed {expected_phyla['Cyanobacteriota']}%" in supplement
    assert (
        f"surface soil ({fmt(100 * cyano['mean_surface'], 2)}%) than in shallow-subsurface "
        f"({fmt(100 * cyano['mean_shallow_subsurface'], 2)}%) or root-adjacent soil "
        f"({fmt(100 * cyano['mean_root_adjacent'], 2)}%)"
    ) in supplement
    over_two = int((phyla["mean_relative_abundance"] > 0.02).sum())
    assert f"{['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine'][over_two]} phyla exceeded 2%" in supplement

    # Supplementary phylum table: every printed row matches the artifact.
    columns = ("mean_relative_abundance", "prevalence", "mean_surface", "mean_shallow_subsurface",
               "mean_root_adjacent", "mean_west_third", "mean_east_third")
    for row in mt.table_rows("supplement", "tab:taxa_phyla")[1:]:
        source = phyla.loc[row[0]]
        assert row[1] == f"{int(source['n_asvs']):,}", row
        expected = [pct(source[c]) if c != "prevalence" else fmt(100 * source[c], 0) for c in columns]
        assert row[2:] == expected, row

    genera = pd.read_csv(RESULTS / "genus_composition.tsv", sep="\t").set_index("genus")
    for genus, printed in (
        ("Domibacillus", "6.2"), ("Massilia", "3.5"), ("Bacillus", "3.0"),
        ("Microvirga", "2.4"), ("Flavisolibacter", "2.0"),
    ):
        assert pct(genera.loc[genus, "mean_relative_abundance"]) == printed
        assert f"{genus} ({printed}%" in main, genus
    unclassified = pct(genera.loc["unclassified_genus", "mean_relative_abundance"])
    assert unclassified == "35.6"
    assert f"ASVs without a genus assignment made up {unclassified}%" in main
    assert f"Reads without a genus assignment made up {unclassified}% of reads" in supplement
    columns = ("mean_relative_abundance", "prevalence", "mean_surface", "mean_shallow_subsurface",
               "mean_root_adjacent", "mean_west_third", "mean_central_third", "mean_east_third")
    for row in mt.table_rows("supplement", "tab:taxa_genera")[1:]:
        source = genera.loc[row[0]]
        assert row[1] == source["phylum"], row
        expected = [fmt(100 * source[c], 0) if c == "prevalence" else fmt(100 * source[c], 2) for c in columns]
        assert row[2:] == expected, row

    leaders = pd.read_csv(RESULTS / "leading_genera_by_stratum.tsv", sep="\t")
    for stratum in ("surface", "shallow_subsurface", "root_adjacent", "west_third", "central_third"):
        top = leaders[(leaders["stratum"] == stratum) & (leaders["rank_in_stratum"] == 1)]
        assert top["genus"].iloc[0] == "Domibacillus", stratum
    east = leaders[(leaders["stratum"] == "east_third") & (leaders["rank_in_stratum"] == 1)]
    assert east["genus"].iloc[0] == "Bacillus"
    assert "Domibacillus (6.2%) was most abundant in every compartment" in main
    assert "Domibacillus was the most abundant genus in every compartment and in the western and central thirds" in supplement
    assert "Bacillus was the most abundant in the eastern third" in supplement

    # The class-level summary remains an artifact check (no class table is printed).
    classes = pd.read_csv(RESULTS / "class_composition.tsv", sep="\t").set_index("class")
    for name, printed in (
        ("Bacilli", "21.0"),
        ("Gammaproteobacteria", "14.7"),
        ("Actinobacteria", "13.3"),
        ("Alphaproteobacteria", "12.2"),
    ):
        assert pct(classes.loc[name, "mean_relative_abundance"]) == printed


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
    main = mt.text("main")
    assert "124 genera were significantly associated with the transect (q<0.05), where 77 decreased and 47 increased" in main
    actino = fmt(100 * 27 / 77, 0)
    bacillota = fmt(100 * 17 / 47, 0)
    assert f"mainly Actinomycetota ({actino}%)" in main and f"mainly Bacillota ({bacillota}%)" in main
    for genus in ("Cellulomonas", "Marmoricola", "Halalkalibacter", "Aquibacillus", "Pseudalkalibacillus"):
        assert f"{genus} (ρ={fmt(rows.loc[genus, 'spearman_rho_route_position'], 2)})" in main, genus
    assert f"western and eastern thirds (ρ={fmt(rows.loc['Bacillus', 'spearman_rho_route_position'], 2)})" in main
    west = pct(rows.loc["Bacillus", "mean_relative_abundance_west_third"])
    east = pct(rows.loc["Bacillus", "mean_relative_abundance_east_third"])
    assert (west, east) == ("1.6", "6.2")
    assert f"Bacillus increased from {west} to {east}% mean relative abundance" in main
    supplement = mt.text("supplement")
    assert "decreasing genera were mainly Actinomycetota (27) and Pseudomonadota (20)" in supplement
    assert "increasing genera were mainly Bacillota (17) and Pseudomonadota (16)" in supplement
    # Supplementary replacement table rows.
    printed = mt.table_rows("supplement", "tab:taxa_replacement")[1:]
    assert len(printed) == 20
    for row in printed:
        source = rows.loc[row[0]]
        assert row[1] == source["phylum"]
        assert row[2] == fmt(source["spearman_rho_route_position"], 2), row
        assert row[3] == _sci(source["q_bh_200"]), row
        assert row[4] == (
            f"{fmt(source['east_minus_west_mean_clr'], 2)} [{fmt(source['east_minus_west_ci_low'], 2)}, "
            f"{fmt(source['east_minus_west_ci_high'], 2)}]"
        ), row
        assert row[5:] == [_share(source["mean_relative_abundance_west_third"]),
                           _share(source["mean_relative_abundance_east_third"])], row

    gradients = pd.read_csv(RESULTS / "site_gradients.tsv", sep="\t").set_index("variable")
    temperature, humidity, ph = (
        fmt(gradients.loc[v, "spearman_rho_route_position"], 2)
        for v in ("mean_air_temperature_c", "mean_relative_humidity_pct", "mean_ph")
    )
    assert f"{temperature} for temperature, {humidity} for humidity and {ph} for pH" in main
    labels = {
        "49-month mean air temperature": "mean_air_temperature_c",
        "49-month mean monthly rainfall": "mean_monthly_rain_mm",
        "49-month mean relative humidity": "mean_relative_humidity_pct",
        "archived-soil pH (equal-weight assay-group means within site)": "mean_ph",
    }
    for row in mt.table_rows("supplement", "tab:route_gradients")[1:]:
        source = gradients.loc[labels[row[0]]]
        assert row[1:] == [
            str(int(source["n_sites"])),
            fmt(source["spearman_rho_route_position"], 2),
            fmt(source["mean_west_third"], 2),
            fmt(source["mean_east_third"], 2),
        ], row


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
    # The genus-set overlap and the between-desert comparison are not part of
    # the current manuscript; no overlap or cross-desert claim is printed.
    combined = mt.combined()
    assert "Jaccard" not in combined
    for value in ("267 of 418", "43 of 397"):
        assert value not in combined


def test_pathway_numbers_match_text() -> None:
    classes = pd.read_csv(RESULTS / "pathway_class_share.tsv", sep="\t").set_index("pathway_class")
    printed = {name: pct(classes.loc[name, "share_of_predicted_pathway_abundance"])
               for name in ("biosynthesis", "energy_central_metabolism", "degradation_utilization")}
    assert printed == {"biosynthesis": "66.8", "energy_central_metabolism": "17.2", "degradation_utilization": "10.7"}
    dominance = pd.read_csv(RESULTS / "pathway_dominance.tsv", sep="\t").set_index("pathway")
    assert dominance.loc["PWY-3781", "rank"] == 1
    assert pct(dominance.loc["PWY-3781", "mean_relative_abundance"]) == "1.7"
    assert dominance.loc["CALVIN-PWY", "rank"] == 44
    assert pct(dominance.loc["CALVIN-PWY", "mean_relative_abundance"]) == "0.6"
    assert dominance.loc["P23-PWY", "rank"] == 122
    assert "PWY-101" not in dominance.index
    supported = pd.read_csv(RESULTS / "supported_pathways_by_class.tsv", sep="\t")
    route = supported[supported["family"] == "route_correlation_supported"]
    assert int(route["n_supported"].sum()) == 92
    biosynthesis = route[route["pathway_class"] == "biosynthesis"].iloc[0]
    assert int(biosynthesis["n_supported"]) == 63 and int(biosynthesis["n_positive"]) == 54
    compartment = supported[supported["family"] == "compartment_contrast_supported"]
    assert int(compartment["n_supported"].sum()) == 270
    main = mt.text("main")
    assert (
        f"Of the {len(dominance)} predicted MetaCyc pathways, biosynthesis pathways accounted for "
        f"{printed['biosynthesis']}% of predicted pathway abundance, energy and central metabolism for "
        f"{printed['energy_central_metabolism']}%, and degradation and utilization for {printed['degradation_utilization']}%"
    ) in main
    assert f"Aerobic respiration (cytochrome c) was the most abundant pathway ({pct(dominance.loc['PWY-3781', 'mean_relative_abundance'])}%)" in main
    rtca = pct(dominance.loc["P23-PWY", "mean_relative_abundance"])
    assert "Calvin–Benson–Bassham cycle ranked 44th (0.6%)" in main
    assert f"reductive TCA cycle 122nd ({rtca}%)" in main
    assert "92 of the 200 pathways were correlated with transect position" in main
    assert "Of 600 pathway–compartment contrasts, 270 were significant" in main
    increased = int(route["n_positive"].sum())
    decreased = int(route["n_negative"].sum())
    assert f"92 were significant (Benjamini–Hochberg q<0.05): {increased} increased and {decreased} decreased eastward" in mt.text("supplement")


def test_generated_tables_are_current_and_included() -> None:
    names = ("taxon_composition_tables.tex", "taxon_overlap_tables.tex", "taxon_pathway_tables.tex")
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp)
        subprocess.run([sys.executable, str(ROOT / "analysis/v3/render_taxon_context_tex.py"),
                        "--results", str(RESULTS), "--output", str(target / "all.tex"),
                        "--split-dir", str(target)], check=True, capture_output=True, text=True)
        # The renderer still reproduces its committed fragments byte for byte.
        assert (target / "all.tex").read_text() == FRAGMENT.read_text()
        for name in names:
            assert (target / name).read_bytes() == (PAPER / "generated" / name).read_bytes()
    # The supplement prints its own composition tables (tables/community_composition.tex);
    # their rows are compared with the same artifacts in the tests above.
    assert r"\input{tables/community_composition.tex}" in mt.strip_comments(SUPPLEMENT)
    referenced = mt.refs("supplement") | {r.removeprefix("si-") for r in mt.refs("main")}
    for label in ("tab:taxa_phyla", "tab:taxa_genera", "tab:taxa_replacement"):
        assert label in mt.labels("supplement")
        assert label in referenced, label


def test_landscape_figure_uses_the_committed_satellite_crop() -> None:
    sidecar = json.loads(
        (ROOT / "metadata/geodata/bluemarble_arabia_200407_120ppd.json").read_text()
    )
    crop = ROOT / "metadata/geodata/bluemarble_arabia_200407_120ppd.png"
    assert hashlib.sha256(crop.read_bytes()).hexdigest() == sidecar["output_sha256"]
    assert sidecar["extent_degrees"] == {"lon_min": 38.0, "lon_max": 63.0, "lat_min": 16.0, "lat_max": 25.0}
    assert sidecar["pixels_per_degree"] == 120
    manifest = json.loads((PAPER / "figures/figure_review_manifest.json").read_text())
    row = manifest["inputs"]["background"]
    assert row["sha256"] == sidecar["output_sha256"]
    assert "NASA Blue Marble" in MAIN


def _sci(value: float) -> str:
    """Two significant figures as m×10^k, the supplement's q format."""
    exponent = int(f"{value:.1e}".split("e")[1])
    return f"{value / 10 ** exponent:.1f}×10^{exponent}"


def _share(value: float) -> str:
    """Percent with two decimals; shares below 0.005% print as <0.01."""
    return "<0.01" if 100 * value < 0.005 else fmt(100 * value, 2)
