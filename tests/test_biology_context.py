"""Biology-context and trait-gene modules: checksums, outputs and manuscript numbers."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import re

import numpy as np
import pandas as pd
from manuscript_paths import PAPER
import manuscript_text as mt
from manuscript_text import fmt, pct

ROOT = Path(__file__).resolve().parents[1]
BIO = ROOT / "analysis/v3/biology_context_corrected_20260909"
TRAIT = ROOT / "analysis/v3/trait_genes_20260909"
PH = ROOT / "analysis/v3/ph_context_reanalysis_20260909/ph_genus_correlations.tsv"
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
    assert f"{omni['route_linear_quadratic']['pseudo_f']:.2f}" == "1.78"
    assert f"{omni['route_linear_quadratic']['permutation_p']:.3f}" == "0.032"
    assert manifest["landform"]["supported_genera_unadjusted"] == 16
    assert manifest["landform"]["supported_genera_route_adjusted"] == 0
    route = pd.read_csv(BIO / "route_model_dune_sensitivity.tsv", sep="\t").set_index("cohort")
    assert f"{100 * route.loc['all_60_sites', 'route_r2']:.1f}" == "40.1"
    assert f"{100 * route.loc['sand_dune_sites', 'route_r2']:.1f}" == "34.2"
    assert int(route.loc["sand_dune_sites", "supported_route_genera_dune_sites"]) == 71
    assert int(route.loc["sand_dune_sites", "supported_in_both_same_sign"]) == 70
    alpha = pd.read_csv(BIO / "alpha_by_landform.tsv", sep="\t").set_index("landform")
    assert f"{alpha.loc['desert oasis', 'mean_shannon']:.2f}" == "4.25"
    assert f"{alpha.loc['sand dune', 'mean_shannon']:.2f}" == "5.77"
    tests = [alpha[f"dune_vs_pan_mannwhitney_p_{m}"].iloc[0] for m in ("shannon", "expected_richness_25k", "normalized_evenness")]
    assert min(tests) >= 0.28
    clim = pd.read_csv(BIO / "climate_diversity_dune_sensitivity.tsv", sep="\t")
    temp_hum = clim[clim["climate_variable"].isin(["mean_air_temperature_c", "mean_relative_humidity_pct"])]
    assert (temp_hum["q_bh_9_dune_44"] < 0.05).all()

    main = mt.text("main")
    assert "44 sand dunes, 9 saline pans (3 western, 1 central, 5 eastern), 3 eastern oases" in main
    assert "pseudo-F=4.56, p=0.0013; 16 genera at q<0.05" in main
    assert "pseudo-F=1.78, p=0.032), and no genera remained differentially abundant" in main
    assert f"did not differ in alpha diversity (p≥{fmt(min(tests), 2)})" in main
    assert "mean 4.25, compared with 5.77 at sand dune sites" in main
    assert "34.2% of compositional variation (40.1% across all sites), and 70 of the 71 genera" in main
    supplement = mt.text("supplement")
    assert "gave p=" + ", ".join(fmt(t, 2) for t in tests[:2]) + f" and {fmt(tests[2], 2)}, respectively" in supplement
    # Landform table rows (Shannon two decimals, richness integer, normalized three).
    names = {"sand dune": "sand dune", "saline pan": "saline pan", "desert oasis": "desert oasis",
             "aeolian lake": "aeolian lake", "dune slack": "dune slack", "gravel": "gravel", "oil spill": "oilspill"}
    rows = mt.table_rows("supplement", "tab:landform_alpha")[1:]
    assert len(rows) == 7
    for row in rows:
        source = alpha.loc[names[row[0]]]
        assert row[1] == str(int(source["n_sites"]))
        assert row[3:] == [fmt(source["mean_shannon"], 2), fmt(source["mean_expected_richness_25k"], 0),
                           fmt(source["mean_normalized_evenness"], 3)], row
    sensitivity = mt.table_rows("supplement", "tab:landform_sensitivity")
    assert ["None", "4.56", "0.0013", "16"] in sensitivity
    assert ["Transect position", "1.78", "0.032", "0"] in sensitivity
    assert ["All (n=60)", "0.401", "", "124"] in sensitivity and ["Dune (n=44)", "0.342", "", "71"] in sensitivity
    # Climate-diversity correlations, all sites and dune sites only.
    variables = {"Temperature": "mean_air_temperature_c", "Humidity": "mean_relative_humidity_pct",
                 "Rainfall": "mean_monthly_rain_mm"}
    measures = {"Shannon": "shannon", "Rarefied richness": "expected_richness_25k",
                "Normalized Shannon": "normalized_evenness"}
    current = None
    checked = 0
    for row in mt.table_rows("supplement", "tab:landform_climate"):
        if row[0] in variables:
            current = variables[row[0]]
        if row[1] not in measures:
            continue
        source = clim[(clim.climate_variable == current) & (clim.diversity_measure == measures[row[1]])].iloc[0]
        assert row[2:] == [fmt(source["spearman_rho_all_60"], 2), _q2(source["q_bh_9_all_60"]),
                           fmt(source["spearman_rho_dune_44"], 2), _q2(source["q_bh_9_dune_44"])], row
        checked += 1
    assert checked == 9
    dune_q = clim.set_index(["climate_variable", "diversity_measure"])["q_bh_9_dune_44"]
    assert (
        f"rainfall–Shannon (q={fmt(dune_q[('mean_monthly_rain_mm', 'shannon')], 3)}) and rainfall–normalized Shannon "
        f"(q={fmt(dune_q[('mean_monthly_rain_mm', 'normalized_evenness')], 3)}) did not"
    ) in supplement


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
    main = mt.text("main")
    assert "(600 tests) supported 295 differences after correction" in main
    for genus in ("Deinococcus", "Rubrobacter", "Kineococcus", "Kocuria", "Neorhizobium", "Devosia",
                  "Pseudomonas", "Brevundimonas", "Bradyrhizobium", "Nitrospira", "MND1", "Gaiella"):
        assert genus in main
    supplement = mt.text("supplement")
    assert (
        "Of the 600 genus-level tests, 295 were supported: 72 for shallow-subsurface versus surface "
        "(33 higher below the surface, 39 higher at the surface), 109 for root-adjacent versus surface "
        "(54 and 55) and 114 for root-adjacent versus shallow-subsurface (53 and 61)"
    ) in supplement
    # Supplementary table: |ΔCLR| and the direction of every listed genus.
    table = mt.section("supplement", r"\labeltab:compartment_genera", r"\endtable")
    printed_q = re.search(r"All genera shown have q≤(\d+\.\d+)", mt.text("supplement")).group(1)
    blocks = {
        ("Deep-Surface", "shallow_subsurface"): "Higher in shallow-subsurface than surface",
        ("Deep-Surface", "surface"): "Higher in surface than shallow-subsurface",
        ("Rhizosphere-Surface", "root_adjacent"): "Higher in root-adjacent than surface",
        ("Rhizosphere-Surface", "surface"): "Higher in surface than root-adjacent",
        ("Rhizosphere-Deep", "root_adjacent"): "Higher in root-adjacent than shallow-subsurface",
        ("Rhizosphere-Deep", "shallow_subsurface"): "Higher in shallow-subsurface than root-adjacent",
    }
    rows = [r for r in mt.table_rows("supplement", "tab:compartment_genera")[1:]]
    sections = re.split(r"\\multicolumn3l", table)[1:]
    assert len(sections) == 6
    for index, ((contrast, compartment), heading) in enumerate(blocks.items()):
        assert sections[index].startswith(heading), (index, heading)
    for row_index, row in enumerate(rows):
        pair = row_index // 8
        for side, (genus, phylum, value) in enumerate((row[0:3], row[3:6])):
            contrast, compartment = list(blocks)[2 * pair + side]
            source = fam[(fam.genus == genus) & (fam.contrast == contrast)].iloc[0]
            assert source["higher_in"] == compartment and source["supported_q_lt_0_05"], (genus, contrast)
            assert value == fmt(abs(source["mean_clr_difference"]), 2), (genus, contrast)
            assert source["q_bh_600"] <= float(printed_q) + 5e-4


def test_gradient_and_core_numbers_match_text() -> None:
    manifest = json.loads((BIO / "run_manifest.json").read_text())
    assert manifest["xrf_axis"]["positive_loading_elements"] == ["Ca", "Mg", "Na", "S", "Cl", "Fe", "Ti"]
    assert manifest["xrf_axis"]["negative_loading_elements"] == ["Si"]
    assert f"{manifest['xrf_axis']['spearman_rho_axis_vs_route']:.2f}" == "0.73"
    assert manifest["xrf_axis"]["supported_genera"] == 118
    assert manifest["ph"]["supported_genera"] == 117
    assert manifest["ph"]["supported_also_route_supported"] == 104
    ph = pd.read_csv(PH, sep="\t").set_index("genus")
    up = ph.loc[["Polygonibacillus", "Sediminibacillus", "Halalkalibacter", "Gracilibacillus", "Halomonas", "Aquibacillus"], "spearman_rho_site_ph"]
    assert up.min() >= 0.705 and up.max() <= 0.815 and ph.loc[up.index, "supported_q_lt_0_05"].all()
    down = ph.loc[["Pirellula", "Roseisolibacter", "Gemmatimonas", "Steroidobacter", "Gaiella", "Sphingomonas", "Streptomyces"], "spearman_rho_site_ph"]
    assert down.max() <= -0.595 and down.min() >= -0.685
    core = manifest["core"]["core_genera_per_compartment"]
    assert core == {"surface": 82, "shallow_subsurface": 69, "root_adjacent": 85}
    assert len(manifest["core"]["core_shared_by_all_three"]) == 59
    table = pd.read_csv(BIO / "core_genera_by_compartment.tsv", sep="\t").set_index("compartment")
    assert [int(table.loc[c, "median_rarefied_genera_per_site"]) for c in ("root_adjacent", "surface", "shallow_subsurface")] == [299, 274, 284]

    def median_text(value: float) -> set[str]:
        # A median of an even number of sites can end in .5; either neighbour
        # integer is an acceptable whole-genus rounding.
        return {str(int(np.floor(value))), str(int(np.ceil(value)))}

    main = mt.text("main")
    supplement = mt.text("supplement")
    decay = {row[0]: row for row in mt.table_rows("main", "tab:decay")[1:]}
    labels = {"surface": "Surface", "shallow_subsurface": "Shallow-subsurface", "root_adjacent": "Root-adjacent"}
    for compartment, label in labels.items():
        assert decay[label][5] == str(core[compartment])
        assert decay[label][6] in median_text(table.loc[compartment, "median_rarefied_genera_per_site"])
    assert "85 genera at ≥ 90% of sites, compared with 82 in surface and 69 in shallow-subsurface soil" in supplement
    assert "59 core genera were shared by all three compartments" in supplement
    assert "The three compartment cores shared 59 genera" in main
    names = {"surface": "surface", "shallow_subsurface": "shallow-subsurface", "root_adjacent": "root-adjacent"}
    for compartment, label in names.items():
        row = mt.table_row("supplement", "tab:core_genera", label)
        source = table.loc[compartment]
        assert row[1:4] == [str(int(source["n_sites_pooled"])), str(int(source["n_genera_detected_any_site"])),
                            str(int(source["n_core_genera_ge_90pct_sites"]))]
        assert row[4] in median_text(source["median_rarefied_genera_per_site"])
        assert row[5] == fmt(source["mean_site_occupancy_of_detected_genera"], 3)

    # pH genera: counts, the named genera and their printed rho ranges.
    assert "117 were correlated with mean site pH (q<0.05), 104 of which were also associated with the transect" in main
    assert "117 correlated with site-mean pH (q<0.05), 104 of which also correlated with transect position" in supplement
    named = re.search(
        r"including ([^()]+?), increased with pH \(ρ=([-\d.]+)–([-\d.]+)\), whereas ([^()]+?) decreased "
        r"\(ρ=([-\d.]+) to ([-\d.]+)",
        main,
    )
    assert named, "pH genus sentence not found"
    increasing, low_up, high_up, decreasing, high_down, low_down = named.groups()
    assert (low_up, high_up) == _named_range(_split_names(increasing), ph, "spearman_rho_site_ph"), (
        "printed range must span the named increasing genera"
    )
    lo, hi = _named_range(_split_names(decreasing), ph, "spearman_rho_site_ph")
    assert (high_down, low_down) == (hi, lo), "printed range must span the named decreasing genera"
    for row in mt.table_rows("supplement", "tab:ph_genera")[1:]:
        for genus, phylum, rho in (row[0:3], row[3:6]):
            assert ph.loc[genus, "phylum"] == phylum
            assert fmt(ph.loc[genus, "spearman_rho_site_ph"], 2) == rho, genus
            assert ph.loc[genus, "supported_q_lt_0_05"]

    # XRF axis: counts and the genus table.
    xrf = pd.read_csv(BIO / "xrf_axis_genus_correlations.tsv", sep="\t").set_index("genus")
    assert int(xrf["supported_q_lt_0_05"].sum()) == manifest["xrf_axis"]["supported_genera"]
    assert f"{manifest['xrf_axis']['supported_genera']} genera correlated with the score" in main
    assert f"Scores increased from west to east along the transect (ρ={fmt(manifest['xrf_axis']['spearman_rho_axis_vs_route'], 2)};" in main
    for row in mt.table_rows("supplement", "tab:xrf-genera")[1:]:
        genus, phylum, rho, q = row
        assert xrf.loc[genus, "phylum"] == phylum
        assert fmt(xrf.loc[genus, "spearman_rho_elemental_axis"], 2) == rho, genus
        assert _q2(xrf.loc[genus, "q_bh_200"]) == q, genus

    taxon = pd.read_csv(ROOT / "analysis/v3/taxon_context/genus_composition.tsv", sep="\t").set_index("genus")
    assert f"{100 * taxon.loc['unclassified_genus', 'mean_relative_abundance']:.1f}" == "35.6"
    phyla = pd.read_csv(ROOT / "analysis/v3/taxon_context/phylum_composition.tsv", sep="\t").set_index("phylum")
    archaea = phyla.loc[[p for p in ("Halobacteriota", "Thermoplasmatota", "Thermoproteota") if p in phyla.index], "mean_relative_abundance"].sum()
    assert f"{100 * archaea:.2f}" == "0.02"
    cyano = phyla.loc["Cyanobacteriota"]
    assert [f"{100 * cyano[c]:.2f}" for c in ("mean_surface", "mean_shallow_subsurface", "mean_root_adjacent")] == ["0.40", "0.41", "0.51"]
    assert "ASVs without a genus assignment made up 35.6% of reads per sample on average" in main
    assert f"Archaea contributed {100 * archaea:.2f}%" in supplement


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

    main = mt.text("main")
    results = mt.section("main", "We examined mechanisms relevant to survival", r"\begintable")
    assert (
        f"{pct(genome.loc['coxL_CO_dehydrogenase', 'fraction_of_genomes'])}% of MAGs carried the aerobic CO dehydrogenase gene coxL "
        f"and {pct(genome.loc['NiFe_hydrogenase_large', 'fraction_of_genomes'])}% a [NiFe]-hydrogenase"
    ) in results
    assert "membrane-bound group 1 or soluble group 3d" in results
    assert f"nifH was only in {pct(genome.loc['nifH_nitrogenase', 'fraction_of_genomes'])}%" in results
    assert f"(ρ={named.min():.2f}-{named.max():.2f} across 119 samples)" in results
    assert "combines hyaB, hhyL and hoxH" in mt.text("supplement")
    # Directions along the transect stated in the Results agree with the q<0.05 rhos.
    for trait in ("otsA_otsB_trehalose", "crtB_phytoene_synthase"):
        assert summary.loc[trait, "spearman_rho_route"] < 0 and summary.loc[trait, "q_bh_route"] < 0.05
    for trait in ("ectABC_ectoine", "spo0A_sporulation", "sodA_superoxide_dismutase"):
        assert summary.loc[trait, "spearman_rho_route"] > 0 and summary.loc[trait, "q_bh_route"] < 0.05
    shares = ("mean_share_surface", "mean_share_shallow_subsurface", "mean_share_root_adjacent")
    assert genome.loc["NiFe_hydrogenase_large", list(shares)].idxmax() == "mean_share_shallow_subsurface"
    for trait in ("crtB_phytoene_synthase", "phrB_photolyase"):
        assert genome.loc[trait, list(shares)].idxmax() == "mean_share_surface", trait
    assert genome.loc["spo0A_sporulation", list(shares)].idxmax() == "mean_share_root_adjacent"

    def star(value: str, q: float) -> str:
        return value + ("*" if q < 0.05 else "")

    # Main Table 3 and Supplementary Tables (genomes, PICRUSt2): every row.
    for row in mt.table_rows("main", "tab:markers")[1:]:
        trait = _trait(row[0])
        g = genome.loc[trait]
        assert row[1:5] == [pct(g["fraction_of_genomes"])] + [pct(g[c]) for c in shares], row
        assert row[5] == star(fmt(summary.loc[trait, "spearman_rho_route"], 2), summary.loc[trait, "q_bh_route"]), row
        assert row[6] == fmt(agree.loc[trait], 2), row
    genome_rows = mt.table_rows("supplement", "tab:traits_genomes")[1:]
    assert len(genome_rows) == 20
    abbreviations = {"Act.": "Actinomycetota", "Pse.": "Pseudomonadota", "Bac.": "Bacillota"}
    for row in genome_rows:
        trait = _trait(row[0])
        g = genome.loc[trait]
        assert row[1:5] == [pct(g["fraction_of_genomes"])] + [pct(g[c]) for c in shares], row
        phylum, count = re.fullmatch(r"(\S+) \((\d+)\)", row[5]).groups()
        assert g["leading_carrier_phyla"].startswith(f"{abbreviations[phylum]} {count}"), row
    picrust_rows = mt.table_rows("supplement", "tab:traits_picrust")[1:]
    assert len(picrust_rows) == 20
    for row in picrust_rows:
        trait = _trait(row[0])
        folds = [contrasts.loc[(trait, c)] for c in ("shallow_subsurface-surface", "root_adjacent-surface",
                                                     "root_adjacent-shallow_subsurface")]
        assert row[1:4] == [star(fmt(f["fold_change"], 2), f["q_bh"]) for f in folds], row
        assert row[4] == star(fmt(summary.loc[trait, "spearman_rho_route"], 2), summary.loc[trait, "q_bh_route"]), row
        assert row[5] == fmt(agree.loc[trait], 2), row
    # MAG carrier shares are descriptive; a compartment "difference" claim
    # for them must not be presented as a tested result.
    assert "Carrier abundance differed among compartments" not in main, (
        "MAG carrier shares by compartment are descriptive (no test)"
    )


def test_generated_fragments_are_current_and_included() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "bio.tex"
        out_traits = Path(tmp) / "trait.tex"
        subprocess.run(
            [sys.executable, str(ROOT / "analysis/v3/render_biology_context_tex.py"),
             "--biology-dir", str(BIO), "--trait-dir", str(TRAIT),
             "--ph-correlations", str(PH), "--output", str(out), "--output-traits", str(out_traits)],
            check=True, capture_output=True, text=True,
        )
        # The renderer still reproduces its committed fragments byte for byte.
        assert out.read_text() == (PAPER / "generated/biology_context_tables.tex").read_text()
        assert out_traits.read_text() == (PAPER / "generated/trait_gene_tables.tex").read_text()
    # The supplement prints its own tables; their rows are compared with the
    # artifacts above. Every one is an active input and is cited.
    active = mt.strip_comments(SUPPLEMENT)
    for name in ("landforms", "pH", "xrf", "beta_diversity", "function_tables"):
        assert f"\\input{{tables/{name}.tex}}" in active, name
    referenced = mt.refs("supplement") | {r.removeprefix("si-") for r in mt.refs("main")}
    for label in ("tab:landform_alpha", "tab:landform_sensitivity", "tab:landform_climate", "tab:compartment_genera",
                  "tab:xrf-genera", "tab:ph_genera", "tab:core_genera", "tab:pathways_compartment",
                  "tab:traits_genomes", "tab:traits_picrust"):
        assert label in mt.labels("supplement"), label
        assert label in referenced, label


def _q2(value: float) -> str:
    """Two significant figures; values below 0.01 as m×10^k (supplement style)."""
    if value >= 0.01:
        digits = 2 - int(np.floor(np.log10(value))) - 1
        return fmt(value, digits)
    exponent = int(f"{value:.1e}".split("e")[1])
    return f"{value / 10 ** exponent:.1f}×10^{exponent}"


_TRAIT_KEYS = {
    "coxL": "coxL_CO_dehydrogenase", "NiFe": "NiFe_hydrogenase_large", "rbcL": "rbcL_RuBisCO",
    "psbA": "psbA_photosystem_II", "nifH": "nifH_nitrogenase", "amoA": "amoA_pmoA_ammonia_monooxygenase",
    "nirK": "nirK_nirS_nitrite_reductase", "nosZ": "nosZ_N2O_reductase", "ureC": "ureC_urease",
    "otsA": "otsA_otsB_trehalose", "treY": "treY_treZ_trehalose", "ectA": "ectABC_ectoine",
    "crtB": "crtB_phytoene_synthase", "spo0A": "spo0A_sporulation", "phrB": "phrB_photolyase",
    "uvrA": "uvrA_excision_repair", "recA": "recA_recombination_repair", "katE": "katE_katG_catalase",
    "sodA": "sodA_superoxide_dismutase", "dps": "dps_DNA_protection",
}


def _trait(label: str) -> str:
    for key, trait in _TRAIT_KEYS.items():
        if key in label:
            return trait
    raise KeyError(label)


def _named_range(genera: list[str], table: pd.DataFrame, column: str) -> tuple[str, str]:
    values = table.loc[genera, column]
    return fmt(values.min(), 2), fmt(values.max(), 2)


def _split_names(text: str) -> list[str]:
    return [name.strip() for name in re.split(r",\s*(?:and\s+)?|\s+and\s+", text) if name.strip()]
