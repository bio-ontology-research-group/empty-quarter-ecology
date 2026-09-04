"""Biology-context and trait-gene modules: checksums, outputs and manuscript numbers."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIO = ROOT / "analysis/v3/biology_context"
TRAIT = ROOT / "analysis/v3/trait_genes"
PAPER = ROOT / "empty-quarter-amplicon"
MAIN = (PAPER / "main.tex").read_text(encoding="utf-8")
SUPPLEMENT = (PAPER / "supplement.tex").read_text(encoding="utf-8")


def flat(text: str) -> str:
    return " ".join(text.split()).replace("$", "").replace("{,}", ",")


def checksums_current(directory: Path) -> None:
    for line in (directory / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        assert hashlib.sha256((directory / name).read_bytes()).hexdigest() == expected, name


def test_checksums_are_current() -> None:
    checksums_current(BIO)
    checksums_current(TRAIT)


def test_landform_numbers_match_text() -> None:
    manifest = json.loads((BIO / "run_manifest.json").read_text())
    sites = manifest["landform"]["sites_per_landform"]
    assert sites["sand dune"] == 44 and sites["saline pan"] == 9 and sites["desert oasis"] == 3
    assert manifest["landform"]["saline_pan_sites_by_third"] == {"east": 5, "west": 3, "central": 1}
    omni = {row["adjustment"]: row for row in manifest["landform"]["omnibus"]}
    assert f"{omni['none']['pseudo_f']:.2f}" == "4.56" and f"{omni['none']['permutation_p']:.4f}" == "0.0013"
    assert f"{omni['route_linear_quadratic']['pseudo_f']:.2f}" == "1.74"
    assert f"{omni['route_linear_quadratic']['permutation_p']:.3f}" == "0.036"
    assert manifest["landform"]["supported_genera_unadjusted"] == 16
    assert manifest["landform"]["supported_genera_route_adjusted"] == 0
    route = pd.read_csv(BIO / "route_model_dune_sensitivity.tsv", sep="\t").set_index("cohort")
    assert f"{100 * route.loc['all_60_sites', 'route_r2']:.1f}" == "40.2"
    assert f"{100 * route.loc['sand_dune_sites', 'route_r2']:.1f}" == "34.2"
    assert int(route.loc["sand_dune_sites", "supported_route_genera_dune_sites"]) == 71
    assert int(route.loc["sand_dune_sites", "supported_in_both_same_sign"]) == 70
    alpha = pd.read_csv(BIO / "alpha_by_landform.tsv", sep="\t").set_index("landform")
    assert f"{alpha.loc['desert oasis', 'mean_shannon']:.2f}" == "4.25"
    assert f"{alpha.loc['sand dune', 'mean_shannon']:.2f}" == "5.77"
    assert alpha["dune_vs_pan_mannwhitney_p_shannon"].iloc[0] >= 0.28
    assert alpha["dune_vs_pan_mannwhitney_p_expected_richness_25k"].iloc[0] >= 0.28
    clim = pd.read_csv(BIO / "climate_diversity_dune_sensitivity.tsv", sep="\t")
    temp_hum = clim[clim["climate_variable"].isin(["mean_air_temperature_c", "mean_relative_humidity_pct"])]
    assert (temp_hum["q_bh_9_dune_44"] < 0.05).all()
    main = flat(MAIN)
    assert "44 sand-dune sites, 9 saline pans (3 in the west, 1 central, 5 in the east), 3 oases in the east" in main
    assert "pseudo-F=4.56, p=0.0013; 16 genera at q<0.05" in main
    assert "pseudo-F=1.74, p=0.036; no genus at q<0.05" in main
    assert "mean Shannon 4.25 against 5.77 at dune sites" in main
    assert "34.2\\,\\% of composition (40.2\\,\\% on all sites), 70 of the 71 genera" in main


def test_compartment_genus_family_matches_text() -> None:
    fam = pd.read_csv(BIO / "compartment_genus_family.tsv", sep="\t")
    assert len(fam) == 600
    sup = fam[fam["supported_q_lt_0_05"]]
    assert len(sup) == 295
    counts = sup.groupby(["contrast", "higher_in"]).size()
    assert counts[("Deep-Surface", "shallow_subsurface")] == 33 and counts[("Deep-Surface", "surface")] == 39
    assert counts[("Rhizosphere-Surface", "root_adjacent")] == 54 and counts[("Rhizosphere-Surface", "surface")] == 55
    assert counts[("Rhizosphere-Deep", "root_adjacent")] == 53 and counts[("Rhizosphere-Deep", "shallow_subsurface")] == 61
    manifest = json.loads((BIO / "run_manifest.json").read_text())
    assert manifest["compartment_family"]["max_abs_deviation_from_committed_loadings"] < 1e-6

    def higher(genus: str, contrast: str, compartment: str) -> bool:
        row = sup[(sup["genus"] == genus) & (sup["contrast"] == contrast)]
        return len(row) == 1 and row["higher_in"].iloc[0] == compartment

    for g in ("Deinococcus", "Kineococcus", "Kocuria", "Planococcus", "Rufibacter"):
        assert higher(g, "Deep-Surface", "surface"), g
    for g in ("Deinococcus", "Rubrobacter", "Geodermatophilus", "Blastococcus"):
        assert higher(g, "Rhizosphere-Surface", "surface"), g
    for g in ("Pelagibacterium", "Neorhizobium", "Devosia", "Metabacillus", "Pseudomonas", "TM7a", "Domibacillus"):
        assert higher(g, "Rhizosphere-Surface", "root_adjacent") and higher(g, "Rhizosphere-Deep", "root_adjacent"), g
    for g in ("Brevundimonas", "Pantoea"):
        assert higher(g, "Rhizosphere-Deep", "root_adjacent"), g
    for g in ("Nitrospira", "MND1"):
        assert higher(g, "Deep-Surface", "shallow_subsurface") and higher(g, "Rhizosphere-Deep", "shallow_subsurface"), g
    for g in ("Gaiella", "Symbiobacterium", "Caldilinea"):
        assert higher(g, "Rhizosphere-Deep", "shallow_subsurface"), g
    assert higher("Bradyrhizobium", "Rhizosphere-Deep", "shallow_subsurface")
    assert higher("Bradyrhizobium", "Rhizosphere-Surface", "surface")
    main = flat(MAIN)
    assert "(600 tests) supported 295 differences after correction" in main
    supplement = flat(SUPPLEMENT)
    assert "supported 295 differences: 72 in the shallow-subsurface--surface contrast (33 higher below the surface, 39 higher at the surface), 109 in the root-adjacent--surface contrast (54 and 55) and 114" in supplement


def test_gradient_and_core_numbers_match_text() -> None:
    manifest = json.loads((BIO / "run_manifest.json").read_text())
    assert manifest["xrf_axis"]["positive_loading_elements"] == ["Ca", "Mg", "Na", "S", "Cl", "Fe", "Ti"]
    assert manifest["xrf_axis"]["negative_loading_elements"] == ["Si"]
    assert f"{manifest['xrf_axis']['spearman_rho_axis_vs_route']:.2f}" == "0.73"
    assert manifest["xrf_axis"]["supported_genera"] == 118
    assert manifest["ph"]["supported_genera"] == 117
    assert manifest["ph"]["supported_also_route_supported"] == 104
    ph = pd.read_csv(BIO / "ph_genus_correlations.tsv", sep="\t").set_index("genus")
    up = ph.loc[["Polygonibacillus", "Sediminibacillus", "Halalkalibacter", "Gracilibacillus", "Halomonas", "Aquibacillus"], "spearman_rho_site_ph"]
    assert up.min() >= 0.705 and up.max() <= 0.815 and ph.loc[up.index, "supported_q_lt_0_05"].all()
    down = ph.loc[["Pirellula", "Roseisolibacter", "Gemmatimonas", "Steroidobacter", "Gaiella", "Sphingomonas", "Streptomyces"], "spearman_rho_site_ph"]
    assert down.max() <= -0.595 and down.min() >= -0.685
    core = manifest["core"]["core_genera_per_compartment"]
    assert core == {"surface": 82, "shallow_subsurface": 69, "root_adjacent": 85}
    assert len(manifest["core"]["core_shared_by_all_three"]) == 59
    table = pd.read_csv(BIO / "core_genera_by_compartment.tsv", sep="\t").set_index("compartment")
    assert [int(table.loc[c, "median_rarefied_genera_per_site"]) for c in ("root_adjacent", "surface", "shallow_subsurface")] == [299, 274, 284]
    main = flat(MAIN)
    assert "root-adjacent soil had 85 genera present at 90\\,\\% or more of the sites, surface soil 82 and shallow-subsurface soil 69" in main
    assert "median number of genera per site was 299, 274 and 284" in main
    assert "117 of the 200 genera correlated with site-mean pH (q<0.05), 104 of them also with route position" in main
    assert "118 genera tracked it" in main
    assert "Spearman \\rho=0.73 with route position" in main
    taxon = pd.read_csv(ROOT / "analysis/v3/taxon_context/genus_composition.tsv", sep="\t").set_index("genus")
    assert f"{100 * taxon.loc['unclassified_genus', 'mean_relative_abundance']:.1f}" == "35.6"
    phyla = pd.read_csv(ROOT / "analysis/v3/taxon_context/phylum_composition.tsv", sep="\t").set_index("phylum")
    archaea = phyla.loc[[p for p in ("Halobacteriota", "Thermoplasmatota", "Thermoproteota") if p in phyla.index], "mean_relative_abundance"].sum()
    assert f"{100 * archaea:.2f}" == "0.02"
    cyano = phyla.loc["Cyanobacteriota"]
    assert [f"{100 * cyano[c]:.2f}" for c in ("mean_surface", "mean_shallow_subsurface", "mean_root_adjacent")] == ["0.40", "0.41", "0.51"]
    assert "Reads without a genus assignment made up 35.6\\,\\% of the total, archaea 0.02\\,\\%" in main


def test_trait_gene_numbers_match_text() -> None:
    manifest = json.loads((TRAIT / "run_manifest.json").read_text())
    assert manifest["genomes"]["analysed"] == 975
    assert manifest["libraries"]["parsed"] == 150 and manifest["libraries"]["matched_to_picrust"] == 119
    genome = pd.read_csv(TRAIT / "genome_trait_summary.tsv", sep="\t").set_index("trait")
    expected = {
        "coxL_CO_dehydrogenase": "38.5", "NiFe_hydrogenase_large": "22.6", "rbcL_RuBisCO": "24.6",
        "psbA_photosystem_II": "0.2", "nifH_nitrogenase": "1.0", "crtB_phytoene_synthase": "49.3",
    }
    for trait, printed in expected.items():
        assert f"{100 * genome.loc[trait, 'fraction_of_genomes']:.1f}" == printed, trait
    assert 70 <= 100 * genome.loc["treY_treZ_trehalose", "fraction_of_genomes"] <= 73
    assert 70 <= 100 * genome.loc["otsA_otsB_trehalose", "fraction_of_genomes"] <= 73
    assert round(100 * genome.loc["ectABC_ectoine", "fraction_of_genomes"]) == 36
    assert round(100 * genome.loc["uvrA_excision_repair", "fraction_of_genomes"]) == 95
    assert round(100 * genome.loc["recA_recombination_repair", "fraction_of_genomes"]) == 92
    assert round(100 * genome.loc["katE_katG_catalase", "fraction_of_genomes"]) == 77
    assert round(100 * genome.loc["sodA_superoxide_dismutase", "fraction_of_genomes"]) == 83
    hyd = genome.loc["NiFe_hydrogenase_large"]
    assert [round(100 * hyd[c]) for c in ("mean_share_shallow_subsurface", "mean_share_surface", "mean_share_root_adjacent")] == [32, 17, 15]
    assert round(100 * genome.loc["crtB_phytoene_synthase", "mean_share_surface"]) == 61
    assert round(100 * genome.loc["phrB_photolyase", "mean_share_surface"]) == 34
    assert genome.loc["NiFe_hydrogenase_large", "leading_carrier_phyla"].startswith("Actinomycetota")
    assert genome.loc["coxL_CO_dehydrogenase", "leading_carrier_phyla"].startswith("Actinomycetota")
    contrasts = pd.read_csv(TRAIT / "picrust_trait_compartment_contrasts.tsv", sep="\t").set_index(["trait", "contrast"])
    assert f"{contrasts.loc[('NiFe_hydrogenase_large', 'root_adjacent-shallow_subsurface'), 'fold_change']:.2f}" == "0.67"
    assert f"{contrasts.loc[('coxL_CO_dehydrogenase', 'root_adjacent-shallow_subsurface'), 'fold_change']:.2f}" == "0.74"
    spo = contrasts.loc[[("spo0A_sporulation", "root_adjacent-surface"), ("spo0A_sporulation", "root_adjacent-shallow_subsurface")], "fold_change"]
    assert 1.4 <= spo.min() and spo.max() < 1.5
    for trait in ("katE_katG_catalase", "sodA_superoxide_dismutase", "dps_DNA_protection"):
        for contrast in ("root_adjacent-surface", "root_adjacent-shallow_subsurface"):
            row = contrasts.loc[(trait, contrast)]
            assert row["fold_change"] > 1 and row["supported_q_lt_0_05"], (trait, contrast)
    summary = pd.read_csv(TRAIT / "picrust_trait_summary.tsv", sep="\t").set_index("trait")
    assert f"{summary.loc['otsA_otsB_trehalose', 'spearman_rho_route']:.2f}" == "-0.71"
    assert f"{summary.loc['crtB_phytoene_synthase', 'spearman_rho_route']:.2f}" == "-0.59"
    assert f"{summary.loc['ectABC_ectoine', 'spearman_rho_route']:.2f}" == "0.36"
    assert f"{summary.loc['spo0A_sporulation', 'spearman_rho_route']:.2f}" == "0.56"
    assert f"{summary.loc['sodA_superoxide_dismutase', 'spearman_rho_route']:.2f}" == "0.70"
    agree = pd.read_csv(TRAIT / "source_agreement.tsv", sep="\t").set_index("trait")["spearman_rho_genome_vs_picrust"]
    named = agree.loc[["coxL_CO_dehydrogenase", "NiFe_hydrogenase_large", "ectABC_ectoine", "otsA_otsB_trehalose", "spo0A_sporulation", "katE_katG_catalase", "sodA_superoxide_dismutase", "dps_DNA_protection"]]
    assert f"{named.min():.2f}" == "0.51" and f"{named.max():.2f}" == "0.76"
    main = flat(MAIN)
    assert "38.5\\,\\% carried the aerobic CO dehydrogenase \\textit{coxL} and 22.6\\,\\% a group 1 [NiFe]-hydrogenase" in main
    assert "hydrogenase carriers made up 32\\,\\% of the recruited community in shallow-subsurface soil against 17\\,\\% at the surface and 15\\,\\% near roots" in main
    assert "Spearman \\rho=0.51--0.76 across 119 matched libraries" in main
    assert "(0.67 and 0.74 of the shallow-subsurface value)" in main


def test_generated_fragments_are_current_and_included() -> None:
    assert "\\input{generated/biology_context_tables.tex}" in SUPPLEMENT
    assert "\\input{generated/trait_gene_tables.tex}" in SUPPLEMENT
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "bio.tex"
        out_traits = Path(tmp) / "trait.tex"
        subprocess.run(
            [sys.executable, str(ROOT / "analysis/v3/render_biology_context_tex.py"), "--output", str(out), "--output-traits", str(out_traits)],
            check=True, capture_output=True, text=True,
        )
        assert out.read_text() == (PAPER / "generated/biology_context_tables.tex").read_text()
        assert out_traits.read_text() == (PAPER / "generated/trait_gene_tables.tex").read_text()
    for label in ("tab:landform-alpha", "tab:landform-sensitivity", "tab:compartment-genera", "tab:xrf-genera", "tab:ph-genera", "tab:core-genera", "tab:compartment-pathways", "tab:trait-genomes", "tab:trait-picrust"):
        assert f"\\ref{{{label}}}" in SUPPLEMENT, label
