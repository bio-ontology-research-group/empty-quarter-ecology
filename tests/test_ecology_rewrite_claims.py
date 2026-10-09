"""Regression checks for the active ecology rewrite and retained archives.

These tests tie the consolidated manuscript prose to the canonical
machine-readable verdicts, including the corrected-coordinate, group-linked
pH and calendar-refitted rainfall outputs. Archived diagnostics remain tested
without requiring their superseded values in current prose.
"""

import hashlib
import json
import re
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from manuscript_paths import PAPER, locked_data_bytes
import manuscript_text as mt
from manuscript_text import fmt, has_number, pct


ROOT = Path(__file__).resolve().parents[1]
MAIN_PATH = PAPER / "main.tex"
SUPPLEMENT_PATH = PAPER / "supplement.tex"
PH_SHARED_PATH = PAPER / "ph_shared_v1.tex"
PH_VALUES_PATH = PAPER / "generated/ph_shared_v1_values.tex"


def _read_manuscript_source(path: Path) -> str:
    if path.exists():
        return path.read_text(encoding="utf-8")
    pytest.fail(f"required active manuscript source is missing: {path}")


@pytest.fixture(scope="module")
def main_tex() -> str:
    return _read_manuscript_source(MAIN_PATH)


@pytest.fixture(scope="module")
def supplement_tex() -> str:
    return _read_manuscript_source(SUPPLEMENT_PATH)


@pytest.fixture(scope="module")
def ph_shared_tex() -> str:
    return _read_manuscript_source(PH_SHARED_PATH)


@pytest.fixture(scope="module")
def ph_values_tex() -> str:
    return _read_manuscript_source(PH_VALUES_PATH)


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def _flat(text: str) -> str:
    return " ".join(text.split())


def _without_value_math(text: str) -> str:
    """Normalize only numeric typesetting for prose-level source checks."""
    return text.replace("$", "").replace("{,}", ",")


def _all_row(relative_path: str) -> pd.Series:
    table = pd.read_csv(ROOT / relative_path, sep="\t")
    return table[
        (table["annotation_basis"] == "ko_copy_count")
        & (table["group"] == "All")
    ].iloc[0]


def test_author_notes_are_resolved_and_no_paragraph_headings_remain(
    main_tex, supplement_tex
):
    combined = main_tex + supplement_tex
    # Unavailable scientific metadata is documented in prose, not TODO boxes.
    todos = re.findall(r"\\todo\{([^}]*)\}", combined)
    assert todos == [], todos
    assert combined.count(r"\todo{") == 0
    assert r"\paragraph{" not in combined
    assert r"\subparagraph{" not in combined
    for author_note in (
        "=>",
        "-> repair",
        "profile'' undefined",
        "Need to confirm that",
        "Why Trip 5?",
    ):
        assert author_note not in combined


def test_author_order_keeps_first_author_then_alphabetical_with_senior_authors_last(
    main_tex,
):
    authors = [
        name.removesuffix(r"\thanks")
        for name in re.findall(r"^\\author(?:\[[0-9,]+\])?\{([^{}]+)", main_tex, re.M)
    ]
    assert authors == [
        "Rund Tawfiq",
        "Marwa Abdelhakim",
        "Sulaiman M. Alajel",
        "Mohammed Alarawi",
        "Hind Aldakhil",
        "Abderahmane Derouiche",
        "Daniela I. Drautz-Moses",
        "Michel Dumontier",
        "Raik Grünberg",
        # Maxat Kulmanov left the author list at his request (30 Aug 2026).
        "Alejandra Lopez-Velazquez",
        "Susana Martinez Arbas",
        "Kexin Niu",
        "Krishnakumar Sivakumar",
        "Tiannyu Wang",
        "Xiang Zhao",
        "Jood Zubair",
        "Magnus Rueping",
        "Robert Hoehndorf",
    ]


def _surname_sort_key(surname: str) -> str:
    """Alphabetical key ignoring LaTeX accents and hyphen/space variants."""
    return (
        surname.replace(r"\\\"", "").replace("ü", "u").replace("-", " ").lower()
    )


def test_data_paper_uses_the_same_author_order_rule():
    source = _read_manuscript_source(ROOT / "data-paper/sn-article.tex")
    authors = re.findall(
        r"^\\author\*?\[[0-9]+\]\{\\fnm\{([^{}]+)\} \\sur\{([^{}]+)\}\}",
        source,
        re.M,
    )
    # The data paper is read from the pinned clone (DATA_REPOSITORY.lock),
    # which lags the live author list (it still names Maxat Kulmanov, who
    # left the ecology author list on 30 Aug 2026).  Compare the ORDER RULE,
    # not an identical list: first author, alphabetical block by surname,
    # then the two senior authors last.
    assert len(authors) >= 17
    assert authors[0] == ("Rund", "Tawfiq")
    assert authors[-2:] == [("Magnus", "Rueping"), ("Robert", "Hoehndorf")]
    middle = [surname for _, surname in authors[1:-2]]
    assert middle == sorted(middle, key=_surname_sort_key), middle
    assert ("Marwa", "Abdelhakim") in authors
    assert ("Raik", r"Gr\"unberg") in authors
    assert ("Xiang", "Zhao") in authors
    assert "Bioscience Core Lab" in source


def test_author_affiliations_match_the_confirmed_institutional_hierarchy(main_tex):
    ecology_by_name = {
        name: affiliation
        for affiliation, name in re.findall(
            r"^\\author\[([0-9,]+)\]\{([^{}]+)", main_tex, re.M
        )
    }
    for name in (
        "Rund Tawfiq",
        "Marwa Abdelhakim",
        "Mohammed Alarawi",
        "Abderahmane Derouiche",
        "Alejandra Lopez-Velazquez",
        "Kexin Niu",
        "Krishnakumar Sivakumar",
        "Jood Zubair",
    ):
        assert "1" in ecology_by_name[name].split(",")
    assert ecology_by_name["Rund Tawfiq"] == "1,2"
    assert ecology_by_name["Susana Martinez Arbas"] == "9"
    assert "Sano Centre for Computational Personalised Medicine" in main_tex
    assert "Maxat Kulmanov" not in ecology_by_name
    # Hind Aldakhil moved to NLFDP/MEWA, Riyadh (her request, 2 Sep 2026).
    assert ecology_by_name["Hind Aldakhil"] == "8"
    assert ecology_by_name["Michel Dumontier"] == "4"
    # Raik Gruenberg's affiliation 7 is now the BioMed Division at KAUST.
    assert ecology_by_name["Raik Grünberg"] == "7"
    assert ecology_by_name["Tiannyu Wang"] == "6"
    assert ecology_by_name["Magnus Rueping"] == "6"
    assert "Bio-Ontology Research Group (BORG)" in main_tex
    assert "Mathematical Sciences and Engineering (CEMSE) Division" in main_tex
    assert "Physical Science and Engineering (PSE) Division" in main_tex
    assert re.search(
        r"^\\affil\[7\]\{Biomedical Sciences Division \(BioMed\)", main_tex, re.M
    )
    assert re.search(
        r"^\\affil\[8\]\{National Livestock \\& Fisheries Development Program \(NLFDP\)",
        main_tex,
        re.M,
    )
    assert "Ministry of Environment, Water and Agriculture (MEWA), Riyadh" in main_tex
    assert "Institute of Data Science, Department of Advanced Computing" in main_tex

    data_source = locked_data_bytes("paper/sn-article.tex").decode("utf-8")
    data_author_pairs = re.findall(
        r"^\\author\*?\[([0-9]+)\]\{\\fnm\{([^{}]+)\} \\sur\{([^{}]+)\}\}",
        data_source,
        re.M,
    )
    data_by_name = {
        f"{given} {surname}": affiliation
        for affiliation, given, surname in data_author_pairs
    }
    # Affiliation numbering differs between papers; read the pinned companion
    # manuscript independently of unrelated edits in the local data checkout.
    for name in (
        "Rund Tawfiq",
        "Marwa Abdelhakim",
        "Mohammed Alarawi",
        "Abderahmane Derouiche",
        "Kexin Niu",
        "Krishnakumar Sivakumar",
    ):
        assert data_by_name[name] == "1"
    for optional_borg in ("Maxat Kulmanov", "Alejandra Lopez-Velazquez",
                          "Alejandra Lopez Velazquez", "Jood Zubair",
                          "Jood Kamal Zubair"):
        if optional_borg in data_by_name:
            assert data_by_name[optional_borg] == "1"
    assert data_by_name["Michel Dumontier"] == "5"
    assert data_by_name["Raik Gr\\\"unberg"] == "6"
    assert data_by_name["Tiannyu Wang"] == "8"
    assert data_by_name["Magnus Rueping"] == "8"
    assert "Bio-Ontology Research Group (BORG)" in data_source
    assert "Mathematical Sciences and Engineering (CEMSE) Division" in data_source
    assert "Physical Science and Engineering (PSE) Division" in data_source
    assert "Biomedical Sciences Division (BioMed)" in data_source
    assert "Institute of Data Science, Department of Advanced Computing" in data_source


def test_active_manuscript_prose_follows_robert_forbidden_word_list():
    # Rund's prose is preserved; stylistic preferences (thus, very,
    # interestingly, utilize, unique) are not enforced. Words that overstate
    # the evidence remain forbidden in claims.
    forbidden = re.compile(
        r"\b(?:clearly|obviously|prove[sn]?|proof|demonstrat(?:e|ed|es|ing)|"
        r"comprehensive|rigorous|definitively|conclusively)\b",
        re.I,
    )
    prose = mt.combined().replace("robust covariance", "sandwich covariance")
    match = forbidden.search(prose)
    assert match is None, prose[max(0, match.start() - 80): match.end() + 80]


def test_main_source_contains_its_prose_and_has_no_tex_fragment_includes(ph_values_tex):
    main_source = mt.strip_comments(_read_manuscript_source(MAIN_PATH))
    preamble, body = main_source.split(r"\begin{document}", 1)
    assert re.search(r"\\(?:input|include|subfile)\s*\{", body) is None
    assert re.search(r"\\(?:input|include|subfile)\s*\{", preamble) is None
    assert r"\newcommand" not in body
    assert re.search(r"\\PH[A-Za-z]+", main_source) is None
    # pH accounting stated in prose matches the frozen shared dataset.
    values = dict(re.findall(r"\\newcommand\{\\(PH[A-Za-z]+)\}\{([^}]*)\}", ph_values_tex))
    summary = _json("analysis/v3/ph_group_linkage_20260909/summary.json")
    counts = summary["counts"]
    dispositions = summary["ingest_dispositions"]
    assert int(values["PHAdmitted"]) == counts["accepted_ph_specimens"] == 712
    assert int(values["PHMatchedSpecimens"]) == counts["group_linked_ph_measurements"] == 709
    assert int(values["PHCompositionGroups"]) == counts["composition_groups"] == 560
    assert dispositions["QUARANTINED_DATE"] == 36 and dispositions["QUARANTINED_QC"] == 19
    measured = 712 + 36 + 19
    assert measured == 767
    assert "We measured pH from 767 stored soil samples" in mt.text("main")
    soil_ph = mt.section("supplement", "Soil pH", "Soil XRF")
    assert "767 measurements, of which 712 met quality-control criteria" in soil_ph
    assert "36 with ambiguous dates and 19 that failed one or more session criteria" in soil_ph
    coverage = mt.section("supplement", r"\labeltab:ph-coverage", r"\endtable")
    rows = re.findall(r"(\d) & ([\d,]+) & ([\d,]+) & ([\d,]+) & ([\d,]+) & ([\d,]+) \\\\", coverage)
    assert len(rows) == 5
    columns = list(zip(*[[int(v.replace(",", "")) for v in row[1:]] for row in rows]))
    target, measured_by, passed, missing, excluded = (sum(c) for c in columns)
    assert (measured_by, passed, excluded) == (767, 712, 55)
    assert target == measured_by + missing
    assert "709 pH measurements were linked to 560 campaign × site × compartment groups" in mt.text("supplement")
    assert "702 profiles" not in mt.combined()


def test_shared_ph_helper_contains_only_generated_scalar_constants(
    ph_shared_tex, ph_values_tex
):
    shared_active = [
        line.strip()
        for line in ph_shared_tex.splitlines()
        if line.strip() and not line.lstrip().startswith("%")
    ]
    assert shared_active == [r"\input{generated/ph_shared_v1_values.tex}"]
    assert r"\newcommand" not in ph_shared_tex
    assert r"\PHEcology" not in ph_shared_tex + ph_values_tex

    value_lines = [
        line.strip()
        for line in ph_values_tex.splitlines()
        if line.strip() and not line.lstrip().startswith("%")
    ]
    assert value_lines
    for line in value_lines:
        match = re.fullmatch(r"\\newcommand\{\\PH[A-Za-z]+\}\{(.+)\}", line)
        assert match is not None, f"non-scalar generated pH command: {line}"
        body = match.group(1)
        assert len(body.split()) == 1, f"prose-bearing generated pH command: {line}"
        assert len(line) < 120


def test_control_method_explains_training_scope():
    audit = _json("analysis/v3/control_audit/summary.json")
    assert len(audit["training_extraction_blanks"]) == 17
    assert audit["mapped_biological_profiles_in_canonical_table"] == 217
    assert audit["mapped_biological_profiles_in_workbook"] == 220
    assert audit["primary_candidate_contaminant_features"] == 351
    assert audit["primary_minimum_blank_prevalence_count"] == 2
    assert audit["primary_prevalence_score_threshold"] == 0.1
    assert audit["positive_controls_in_training"] == 0
    t4 = json.loads((PAPER / "validation/t4-control-diversity/summary.json").read_text())
    assert t4["n_profiles"] == 95 and t4["called_asvs"] == 7
    kit = mt.section("supplement", "DNA extraction kit controls", "Microbial community standard controls")
    assert "present in at least 2 blanks" in kit and "Fisher test (p<0.1)" in kit
    assert "17 extraction blanks linked to 217 samples" in kit
    assert "351 candidate contaminant ASVs" in kit
    assert "6 extraction blanks linked to 95 samples" in kit and "7 candidate contaminant ASVs" in kit
    assert f"ρ={fmt(t4['shannon_spearman'], 5)}" in kit
    table = mt.section("supplement", r"\labeltab:controls", r"\endtable")
    assert "17 mapped T5 blanks linked to 217 biological profiles" in table
    assert "Six blanks linked to 95 biological profiles" in table
    # Positive standards estimate genus recovery; they do not train the screen.
    assert "used for genus recovery" in table
    main = mt.text("main")
    assert "217 samples from Trip 5" in main
    assert "Control samples did not cover all campaigns and sample processing stages" in main


def test_geography_first_results_order_and_regional_novelty_opener():
    introduction = mt.section("main", r"\section*Introduction", r"\section*Results")
    assert "world's largest continuous sand desert" in introduction
    assert "has not yet been surveyed across its landscape or over repeated visits" in introduction
    # Unsupported novelty or population claims stay out of the opener.
    assert "The only direct microbial study" not in introduction
    assert "38% of the world's population" not in introduction
    assert not re.search(r"\b(?:the )?first (?:survey|study|description)\b", introduction, re.I)
    headings = re.findall(r"\\subsection\*\{([^}]+)\}", mt.source("main"))
    assert headings
    assert all("baseline" not in heading.lower() for heading in headings)
    assert not any(heading.startswith("No ") for heading in headings)
    # Results lead with composition and geography before environment and function.
    results = mt.section("main", r"\section*Results", r"\section*Discussion")
    order = [results.index(marker) for marker in (
        "Microbial communities change across the landscape",
        "Compartment differences in composition and diversity within sites",
        "Climate covaries with the geographic pattern",
        "Predicted functional potential follows geography and compartment",
    )]
    assert order == sorted(order)


def test_abstract_leads_with_science_and_keeps_resource_subordinate():
    flat = mt.abstract()
    assert "dissimilarity increased with geographic distance in all sampled compartments" in flat
    assert "Replacement of taxa between sites accounted for most of the dissimilarity" in flat
    assert "root-adjacent soil had the lowest alpha diversity and harbored communities distinct from bulk soil" in flat
    assert "baseline for the Empty Quarter soil microbiome" in flat
    assert "knowledge graph" not in flat.lower()
    assert flat.index("geographic distance") < flat.index("linked data resource")
    turnover = pd.read_csv(
        ROOT / "analysis/v3/distance_decay_turnover/turnover_nestedness_components.tsv", sep="\t"
    )
    assert (turnover["turnover_share_of_sorensen"] > 0.5).all()


def test_pH_attenuation_is_bounded_and_negative_diagnostics_are_supplementary():
    summary = _json("analysis/v3/ph_group_linkage_20260909/summary.json")
    counts = summary["counts"]
    assert counts["group_linked_ph_measurements"] == 709
    assert counts["site_campaign_position_groups"] == 563
    assert counts["composition_groups"] == 560
    composition = summary["composition_primary"]
    geography = _json("analysis/v3/ph_group_linkage_20260909/summary.json")["geographic_same_cohort"]
    # The two archived residual-response R2 values have different
    # denominators; their difference is not an explained fraction.
    assert geography["legacy_partial_r2_fields_are_residual_response_r2"] is True
    partition_dir = ROOT / "analysis/v3/ph_partition_20261004/results"
    partitions = pd.read_csv(partition_dir / "partitions.tsv", sep="\t").set_index("analysis")
    site = partitions.loc["site_averaged_common_adjustment"]
    grouped = partitions.loc["primary_equal_site_weight"]
    fractions = ("unique_ph", "unique_geography", "shared", "unexplained")
    assert abs(sum(site[f] for f in fractions) - 1) < 1e-9
    main = mt.section("main", "Soil pH tracks the west-to-east change", "Soil chemical composition")
    supplement = mt.section("supplement", "Soil pH\\labelsup:beta_ph", "Laboratory XRF elemental axis")
    printed_site = [pct(site[f], 2) for f in fractions]
    for text in (main, supplement):
        for value in printed_site:
            assert has_number(text, value), value
    for value in (pct(grouped[f], 2) for f in fractions):
        assert has_number(supplement, value), value
    assert "common variance denominator" in main and "common response denominator" in supplement
    assert "about half" not in main and "accounted for about half" not in mt.combined()
    # The conditional within-site association is small and not significant.
    assert f"partial RDA, {pct(composition['partial_r2'], 2)}% of variation, p={fmt(composition['p'], 3)}" in main
    assert f"pH explained {pct(composition['partial_r2'], 3)}% of the remaining compositional variation (p={fmt(composition['p'], 3)})" in supplement
    # Sensitivity ranges for the shared component and leave-one-site-out.
    between = pd.read_csv(partition_dir / "between_site_bootstrap_intervals.tsv", sep="\t").set_index(["method", "fraction"])
    whole = between.loc[("site_cluster", "shared")]
    block = between.loc[("noncircular_block_15", "shared")]
    assert f"ranged from {pct(whole['percentile_2_5'], 2)} to {pct(whole['percentile_97_5'], 2)}%" in supplement
    assert f"from {pct(block['percentile_2_5'], 2)} to {pct(block['percentile_97_5'], 2)}% with 15-site blocks" in supplement
    loo = pd.read_csv(partition_dir / "leave_site_out.tsv", sep="\t")["between_site_unique_ph"]
    assert f"{pct(loo.min(), 2)}–{pct(loo.max(), 2)}%" in supplement
    assert "analysis/v3/ph_partition_20261004/" in supplement
    assert "fig:environment-spatial" not in mt.source("main")


def test_paired_composition_claims_match_the_canonical_verdict():
    verdict = _json("analysis/v3/compartment_composition/claim_verdict.json")
    omnibus = verdict["omnibus"]
    contrasts = verdict["contrasts"]

    assert omnibus["n_sites"] == 60
    assert omnibus["n_blocks"] == 170
    assert round(omnibus["pseudo_f"], 2) == 17.28
    assert omnibus["permutation_p"] == 0.001
    main = mt.section("main", "Compartment differences in composition and diversity within sites", "We next identified the genera")
    supplement = mt.section("supplement", "Compartment differences\\labelsup:beta_compartment", "Genus-level contrasts")
    for text in (main, supplement):
        assert f"pseudo-F={omnibus['pseudo_f']:.2f}, p={fmt(omnibus['permutation_p'], 3)}" in text
        assert "170" in text
    labels = {
        "Deep-Surface": "shallow-subsurface versus surface",
        "Rhizosphere-Surface": "root-adjacent versus surface",
        "Rhizosphere-Deep": "root-adjacent versus shallow-subsurface",
    }
    for contrast, label in labels.items():
        result = contrasts[contrast]
        assert result["permutation_p"] == 0.001
        assert result["q_within_primary_family"] == 0.001
        assert result["direction_stable"] and result["leave_one_campaign_supported"]
        assert f"{result['standardized_displacement']:.3f} ({label})" in supplement
        assert has_number(main, fmt(result["standardized_displacement"], 2))
        assert f"{result['n_blocks']} blocks for {label}" in supplement or (
            f"{result['n_blocks']} for {label}" in supplement
        )
    assert "every pair of compartments differed in composition (q=0.001" in main
    assert "fig3_soil_position.pdf" in mt.source("main")
    methods = mt.section("main", "Compartment differences.", "Turnover and nestedness")
    assert "paired sign flip permutation tests (999) on the mean within-site difference vector" in methods
    assert "assume sign symmetry of the site-level difference vectors under the null" in main
    # Compartments are operational: root distance and host plants were not recorded.
    assert "distance from root were not recorded" in mt.text("main")


def test_evenness_decomposition_is_numerically_and_semantically_bounded():
    verdict = _json("analysis/v3/evenness_decomposition/claim_verdict.json")
    results = verdict["primary_results"]
    root_surface = results["Rhizosphere-Surface"]["evenness_sensitivity"]
    root_shallow = results["Rhizosphere-Deep"]["evenness_sensitivity"]

    assert f"{root_surface['mean_difference']:.5f}" == "-0.03148"
    assert f"{root_shallow['mean_difference']:.5f}" == "-0.04398"
    assert verdict["input"]["blocks_with_evenness_sensitivity"] == 617
    for value in ("-0.03148", "-0.04398"):
        assert value not in mt.text("main")
    # The index is defined explicitly as a mixed-depth ratio.
    assert "H/\\log(E[S_25k])" in mt.text("main")
    indices = mt.section("supplement", "Diversity indices and rarefaction", "Marginal and paired compartment contrasts")
    assert "H/\\log(E[S_25k])" in indices
    assert "mixed-depth ratio" in indices and "can exceed 1" in indices
    manifest = _json("analysis/v3/paired_alpha_sensitivity_20260909/manifest.json")
    low, high = manifest["normalized_shannon_range"]
    assert f"ranged from {fmt(low, 3)} to {fmt(high, 3)}" in indices
    assert f"ρ = {fmt(manifest['normalized_shannon_depth_spearman'], 3)}" in indices
    assert "617 of the 633" not in mt.text("main")


def test_geographic_transport_detail_is_relocated_and_interpretation_bounded(
    main_tex, supplement_tex
):
    verdict = _json("analysis/v3/geographic_prediction/claim_verdict.json")
    assert verdict["primary_arm_supported"] is False
    assert verdict["sensitivity_arm_supported"] is True
    assert round(verdict["group_level_equal_weight_skill"], 4) == -0.0017
    assert round(verdict["group_level_pooled_skill"], 4) == -0.0005
    assert round(verdict["site_level_block_skill"], 4) == 0.2528
    assert verdict["strictest_group_level_null_p_value"] == 0.406
    assert verdict["n_group_level_folds_with_positive_skill"] == 10
    assert verdict["n_group_level_folds"] == 18

    # Only the corrected-coordinate transport diagnostic is reported.
    alias = verdict["collection_order_alias"]
    assert alias["campaigns_with_abs_rho_at_least_0_99"] == 5
    for stale in ("archived pre-correction prediction", "earlier coordinate version"):
        assert stale not in _flat(supplement_tex)
    assert "collection order" in _flat(main_tex)
    for legacy in ("$R^{2}=-0.0017$", "$R^{2}=0.2528$", "-0.0015", "0.2526"):
        assert legacy not in main_tex


def test_moran_claim_is_bounded_to_the_tested_neighbourhood_scale():
    verdict = _json(
        "analysis/v3/spatial_resolution_sensitivity/claim_verdict.json"
    )
    table = pd.read_csv(
        ROOT / "analysis/v3/spatial_resolution_sensitivity/moran_k_sensitivity.tsv",
        sep="\t",
    )
    assert verdict["moran_k_status"] == (
        "residual_autocorrelation_depends_on_neighbour_count"
    )
    assert table["neighbours_k"].tolist() == [3, 4, 5, 6, 8, 10]
    assert round(table.iloc[0]["residual_moran_i"], 4) == 0.1101
    assert round(table.iloc[-1]["residual_moran_i"], 4) == -0.0098
    assert table.iloc[-1]["permutation_p"] == 0.243
    assert verdict["neighbour_counts_with_detected_autocorrelation"] == [3, 4, 5, 6]
    assert verdict["neighbour_counts_without_detected_autocorrelation"] == [8, 10]
    primary = pd.read_csv(
        ROOT / "analysis/v3/spatial_turnover_rescue/results/spatial_model_results.tsv", sep="\t"
    ).query("analysis == 'primary' and taxon_count == 200 and trend_degree == 2").iloc[0]
    moran = f"Moran's I={fmt(primary['residual_moran_i'], 3)}, p={fmt(primary['residual_moran_p'], 3)}"
    # The Moran statistic is reported with the neighbourhood it was tested on.
    assert moran in mt.text("main")
    assert moran + ", five-nearest-neighbour graph" in mt.text("supplement")
    assert "5-nearest neighbor graph" in mt.text("main")
    combined = mt.combined()
    for stale in ("0.1104", "-0.0097", "p=0.241 at", "through k=8", "earlier coordinate version"):
        assert stale not in combined

    # Inference under specified spatial covariances (73 models).
    current = _json("analysis/v3/spatial_covariance_sensitivity_20260909/summary.json")
    assert current["model_count"] == 73
    assert current["original_descriptive_r2"] == pytest.approx(0.40073421629965933)
    assert current["family_p_min"] == 0.0001 and current["family_p_max"] == 1.0
    errors = mt.section("supplement", "Residual spatial autocorrelation", r"\\begintable")
    # Neighbour-count sensitivity of the residual Moran statistic.
    detected = table[table.neighbours_k.isin([3, 4, 5, 6])]
    assert (
        f"With 3–6 neighbours, I={fmt(detected.residual_moran_i.min(), 3)}–"
        f"{fmt(detected.residual_moran_i.max(), 3)} (p≤{fmt(detected.permutation_p.max(), 3)})"
    ) in errors
    k8, k10 = (table.set_index("neighbours_k").loc[k, "permutation_p"] for k in (8, 10))
    assert f"with 8 and 10 neighbours, p={fmt(k8, 3)} and {fmt(k10, 3)}" in errors
    assert "this gave 73 models" in errors
    tests = pd.read_csv(
        ROOT / "analysis/v3/spatial_covariance_sensitivity_20260909/spatial_covariance_tests.tsv", sep="\t"
    ).dropna(subset=["range_km"])
    table_s = mt.section("supplement", r"\labeltab:spatial_covariance", r"\endtable")
    kernels = {"exponential": "Exponential", "matern_3_2": "Mat\\'ern 3/2", "matern_5_2": "Mat\\'ern 5/2"}
    for kernel, label in kernels.items():
        block = table_s.split(label + " &", 1)[1].split(r"\addlinespace", 1)[0]
        for range_km, rows in tests[tests.kernel == kernel].groupby("range_km"):
            line = block.split(f"{range_km:.1f} &", 1)[1].split(r"\\", 1)[0]
            printed = [cell.strip() for cell in line.split("&")]
            expected = [f"{p:.4f}" for p in rows.sort_values("group_noise_fraction")["rotation_p"]]
            assert printed == expected, (kernel, range_km)
    shortest = tests.range_km.min()
    supported = tests[(tests.range_km == shortest) | (tests.group_noise_fraction == 0.9)]
    assert (supported.rotation_p == 0.0001).all()
    assert f"supported (p=0.0001) at the {shortest:.1f}-km range and whenever 90% of variance was independent noise" in errors
    half = tests[(tests.group_noise_fraction == 0.5) & (tests.range_km > 200)]
    assert (
        f"With 50% noise, p ranged from {half.rotation_p.min():.4f} to {fmt(half.rotation_p.max(), 3)} "
        f"at ranges of {half.range_km.min():.0f} km or more"
    ) in errors
    low = tests[(tests.group_noise_fraction <= 0.1) & (tests.range_km > 50)]
    assert (
        f"with ≤10% noise and ranges of {round(low.range_km.min(), -1):.0f} km or more, "
        f"p ranged from {low.rotation_p.min():.4f} to {fmt(low.rotation_p.max(), 3)}"
    ) in errors
    # Main text: support is lost at ranges of about 50 km, not only over
    # hundreds of kilometres.
    statement = mt.section("main", "When we accounted for this spatial autocorrelation", "We therefore report")
    assert "hundreds of kilometres" not in statement
    assert f"ranges of {round(low.range_km.min(), -1):.0f} km or more" in statement


def test_landscape_figure_contains_three_survey_panels(main_tex):
    figure_script = (
        ROOT / "analysis/v3/make_submission_figures.py"
    ).read_text(encoding="utf-8")
    study_figure = figure_script.split(
        "def make_landscape_figure", 1
    )[1].split("def make_soil_position_figure", 1)[0]

    # Figure 1 after the October 2026 revision: map, samples per campaign and
    # climate along the transect; distance decay moved to Figure 2 and the
    # climate associations to Figure 4.
    for content in ("map_ax", "coverage_ax", "climate_ax"):
        assert content in study_figure
    assert "make_composition_geography_figure" in figure_script
    assert "make_environment_gradient_figure" in figure_script
    assert "Explicit analysis cohorts" not in study_figure
    assert "fig:overview-revised" not in main_tex
    assert "(c) The samples available" not in main_tex


def test_xrf_non_detection_is_not_written_as_evidence_of_absence():
    claim = _json("analysis/v3/xrf_community_rescue/xrf_claim_summary.json")
    assert round(claim["alpha_primary"]["p"], 2) == 0.99
    combined = mt.combined()
    # A null Shannon association is not reported as evidence of absence.
    for phrase in (
        "did not track within-site Shannon diversity",
        "Its adjusted association with Shannon diversity was null",
        "It had no conditional association with Shannon diversity",
        "no adjusted Shannon association;",
    ):
        assert phrase not in combined
    # The reported association is the bounded within-site composition model.
    assert "Within sites, the elemental axis explained a small but significant share" in mt.text("main")
    assert "less than 1% of compositional variation within sites" in mt.text("main")


def test_depth_adjusted_claims_match_the_canonical_verdict():
    verdict = _json("analysis/v3/depth_extraction/claim_verdict.json")
    interaction = verdict["campaign_by_position_interaction"]
    supported = verdict["contrasts"]["Rhizosphere-Deep"]
    direction_only = verdict["contrasts"]["Deep-Surface"]
    sensitivity_dependent = verdict["contrasts"]["Rhizosphere-Surface"]

    assert f"{interaction['unadjusted_wald_p']:.5f}" == "0.00831"
    assert f"{interaction['depth_adjusted_wald_p']:.5f}" == "0.17476"
    assert f"{supported['depth_adjusted_estimate']:.3f}" == "-0.273"
    assert f"{supported['depth_adjusted_ci'][0]:.3f}" == "-0.424"
    assert f"{supported['depth_adjusted_ci'][1]:.3f}" == "-0.122"
    assert f"{direction_only['depth_adjusted_estimate']:.3f}" == "0.103"
    assert sensitivity_dependent["status"] == "sensitivity_dependent"
    assert all(c["direction_stable_across_models"] for c in verdict["contrasts"].values())
    assert sum(sensitivity_dependent["interval_excludes_zero_by_model"].values()) == 3
    assert all(supported["interval_excludes_zero_by_model"].values())
    assert verdict["join_audit"]["profiles_with_recorded_kit"] == 858
    confounds = mt.section("supplement", "Potential confounds", "Climate associations")
    assert "root-adjacent versus shallow subsurface contrast remained significant in every model" in confounds
    assert "root-adjacent versus surface contrast in only three of the six" in confounds
    assert "extraction kit was recorded for 858 of the samples" in confounds
    assert "partly confounded with campaign" in confounds
    # Normalized Shannon was adjusted for campaign and read depth only (no kit term).
    evenness = pd.read_csv(ROOT / "analysis/v3/evenness_decomposition/evenness_depth_sensitivity.tsv", sep="\t")
    assert not evenness["formula"].str.contains("kit").any()
    assert evenness["formula"].str.contains("log_sequencing_depth").any()
    sentence = re.search(r"[^.]*normalized Shannon contrasts[^.]*\.", confounds, re.I)
    assert sentence and "depth" in sentence.group(0), (
        "state that normalized Shannon contrasts were adjusted for campaign and read depth (not kit)"
    )
    # Main text: direction is stable; no claim beyond direction for all endpoints.
    main = mt.text("main")
    assert "did not change the direction of any contrast" in main
    assert "expected richness retained a depth-adjusted interaction" not in mt.combined()


def test_consolidation_removes_untraceable_between_site_xrf_numbers():
    flat = mt.text("main")
    assert "PC1 and Shannon diversity had ρ=-0.68" not in flat
    assert "partial ρ=-0.30" not in flat
    assert "block size increased from 3 to 20 sites" not in flat
    verdict = _json("analysis/v3/xrf_community_clr/claim_verdict.json")["primary"]
    assert f"{pct(verdict['partial_r2'], 2)}%" in flat
    assert f"{pct(verdict['partial_r2'], 3)}%" in mt.text("supplement")


def test_campaign_omission_is_bounded_as_an_influence_analysis():
    cohort = pd.read_csv(
        ROOT / "analysis/v3/compartment_composition/cohort_accounting.tsv",
        sep="\t",
    )
    retained = cohort[cohort["retained"] == True]  # noqa: E712
    counts = retained.groupby("campaign").size()
    shares = counts / counts.sum() * 100
    assert counts.tolist() == [159, 16, 169, 176, 110]
    assert [round(value, 1) for value in shares] == [
        25.2,
        2.5,
        26.8,
        27.9,
        17.5,
    ]
    preparation = mt.section("supplement", "Data preparation, resampling and multiple testing", "Transect model")
    assert "Campaigns contributed " + ", ".join(fmt(v, 1) for v in shares[:-1]) + f" and {fmt(shares.iloc[-1], 1)}%" in preparation
    spatial = pd.read_csv(
        ROOT / "analysis/v3/spatial_turnover_rescue/results/spatial_model_results.tsv", sep="\t"
    )
    loco = spatial[spatial.analysis == "leave_one_campaign_out"]
    assert len(loco) == 5 and (loco.permutation_p <= 0.001).all()
    assert "omitted any one campaign" in mt.text("main")
    assert "after excluding each campaign in turn; all fits gave p=0.001" in mt.text("supplement")


def test_assay_aware_control_filter_is_bounded_and_headlines_are_stable():
    audit = _json("analysis/v3/control_audit/summary.json")
    sensitivity = _json(
        "analysis/v3/control_audit/sensitivity_inputs/summary.json"
    )
    headlines = _json(
        "analysis/v3/control_sensitivity/headline_result_sensitivity.json"
    )
    spillover = _json(
        "analysis/v3/control_audit/positive_control_spillover_summary.json"
    )

    assert audit["canonical_features"] == 351472
    assert audit["primary_candidate_contaminant_features"] == 351
    assert audit["mapped_biological_profiles_in_canonical_table"] == 217
    assert len(audit["training_extraction_blanks"]) == 17
    assert audit["positive_controls_in_training"] == 0
    assert audit["positive_control_profiles"] == 7
    assert sensitivity["mapped_profiles"] == 217
    assert sensitivity["removed_read_fraction"]["pooled"] == pytest.approx(
        0.02186632117803215
    )
    assert sensitivity["removed_read_fraction"]["median"] == pytest.approx(
        0.004027353747150651
    )
    assert sensitivity["removed_read_fraction"]["maximum"] == pytest.approx(
        0.5659508027771022
    )
    assert sensitivity["shannon"]["spearman_before_after"] == pytest.approx(
        0.994959530620969
    )
    assert sensitivity["profiles_below_rarefaction_depth_after_filter"] == 0
    assert headlines["headline_metrics_compared"] == 25
    assert headlines["all_headline_verdicts_stable"] is True
    assert headlines["verdict_changes"] == []
    assert spillover["interpretation_limit"].startswith(
        "Exact-ASV overlap does not distinguish"
    )
    removed = sensitivity["removed_read_fraction"]
    main = mt.text("main")
    kit = mt.section("supplement", "DNA extraction kit controls", "Microbial community standard controls")
    assert f"{pct(removed['pooled'], 2)}% of pooled reads in 217 samples from Trip 5" in main
    assert f"{pct(removed['pooled'], 1)}% of pooled reads and a median of {pct(removed['median'], 2)}% per profile" in kit
    assert f"maximum {pct(removed['maximum'], 2)}%" in kit
    assert f"ρ = {fmt(sensitivity['shannon']['spearman_before_after'], 3)}" in main
    assert "no profile fell below the rarefaction depth" in mt.text("supplement")
    headline = pd.read_csv(ROOT / "analysis/v3/control_sensitivity/headline_result_sensitivity.tsv", sep="\t")
    shannon_q = headline[(headline.claim == "paired Shannon distribution") & (headline.metric == "Rhizosphere-Surface_q")].iloc[0]
    assert f"from {fmt(shannon_q['canonical_value'], 3)} to {fmt(shannon_q['control_adjusted_value'], 3)}" in mt.text("supplement")
    combined = mt.combined()
    for stale in ("14,822", "2.8684"):
        assert stale not in combined
    assert has_number(main, "351,472")
    assert "retained all ASVs in downstream analyses" in main


def test_functional_top_k_sensitivity_matches_all_three_canonical_tables():
    rows = [
        _all_row(
            "analysis/v3/functional_redundancy_sensitivity/"
            "results-k1000-p999/functional_redundancy_null.tsv"
        ),
        _all_row(
            "analysis/v3/functional_redundancy_results/"
            "functional_redundancy_null.tsv"
        ),
        _all_row(
            "analysis/v3/functional_redundancy_sensitivity/"
            "results-k5000-p999/functional_redundancy_null.tsv"
        ),
    ]
    observed = [f"{row['observed_functional_median_bray']:.4f}" for row in rows]
    nulls = [f"{row['null_median']:.4f}" for row in rows]
    assert observed == ["0.1005", "0.1339", "0.1721"]
    assert nulls == ["0.0671", "0.0873", "0.1174"]
    assert all(row["upper_tail_p"] == 0.001 for row in rows)
    # Observed functional dissimilarity exceeds the null in every top-k set, so
    # no functional-redundancy claim is supported; none is made.
    assert all(row["observed_functional_median_bray"] > row["null_median"] for row in rows)
    assert "functional redundancy" not in mt.combined().lower()


def test_distance_decay_is_surfaced_and_uses_whole_site_permutations():
    current = _json("analysis/v3/distance_decay_turnover/claim_verdict.json")
    assert current["site_pairs"] == 1770
    assert current["matched_sites"] == 60
    assert current["standardised_depth"] == 12865
    assert current["permutations"] == 9999
    assert current["omnibus_p"] == 0.0051
    slopes = pd.read_csv(ROOT / "analysis/v3/distance_decay_turnover/distance_decay_slopes.tsv", sep="\t")
    supplement = mt.text("supplement")
    decay = mt.section("supplement", "Distance decay\\labelsup:beta_decay", "Compartment differences")
    main_table = mt.section("main", r"\labeltab:decay", r"\endtable")
    main = mt.text("main")
    assert "(1,770 site pairs)" in decay
    for row in slopes[slopes.family == "aitchison"].itertuples():
        assert f"{row.slope_per_100km:.3f}" in decay
        interval = (
            f"{fmt(row.slope_per_100km, 2)} [{fmt(row.jackknife_ci_low_per_100km, 2)}, "
            f"{fmt(row.jackknife_ci_high_per_100km, 2)}]"
        )
        assert interval in main_table, interval
        assert row.two_sided_p == 0.0001
    assert "p=0.0051" in decay and "p=0.0051" in main
    contrast = slopes.query("family == 'contrast' and response == 'Deep-Rhizosphere'").iloc[0]
    assert contrast.max_t_adjusted_p == 0.0039
    assert "p_adj=0.0039" in decay and "p_adj=0.004" in main
    other = slopes.query("family == 'contrast' and response == 'Surface-Rhizosphere'").iloc[0]
    assert "p_adj=0.173" in main and other.max_t_adjusted_p == pytest.approx(0.1728)
    for row in slopes[slopes.family.isin(["simpson_turnover", "nestedness"])].itertuples():
        assert f"{row.slope_per_100km:.4f}" in supplement
        assert f"({row.slope_per_100km:.4f})" in main_table
    shares = pd.read_csv(
        ROOT / "analysis/v3/distance_decay_turnover/turnover_nestedness_components.tsv", sep="\t"
    )["turnover_share_of_sorensen"]
    span = f"{pct(shares.min(), 0)}–{pct(shares.max(), 0)}%"
    assert span in main and span + " of Sørensen dissimilarity" in supplement
    methods = mt.section("main", "Distance decay.", "Compartment differences.")
    assert "permuting site labels (9,999 permutations) rather than individual distances" in methods
    assert "applied each permutation to all three compartment matrices at once" in methods


def test_new_methodological_citations_have_byte_verifiable_source_custody():
    custody = pd.read_csv(
        ROOT / "literature/CITATION_SOURCES.tsv", sep="\t", dtype=str
    ).set_index("cite_key")
    expected = {
        "baselga2010partition": (
            "https://doi.org/10.1111/j.1466-8238.2009.00490.x",
            "e7c5d191ee0fc84d1de1bea2b3f4ac395a6be3853b4770ba8fba548efd7a9c0f",
        ),
        "klindworth2013primers": (
            "https://doi.org/10.1093/nar/gks808",
            "d04632f48140b963b2c17184fc188f4b659b5003f7802c963d6f1d222b7bbbc0",
        ),
        "guillot2013mantel": (
            "https://doi.org/10.1111/2041-210X.12018",
            "66844824822228bbb5141db10aa0dcd7ecbd5e75e1d29c9a9a60cff2202163a3",
        ),
        "gloor2017microbiome": (
            "https://doi.org/10.3389/fmicb.2017.02224",
            "32cf8e632d4648fabd1bf02124891f403299df5f205cc2a5ff95c55394e2a0f1",
        ),
        "kurtz2015spieceasi": (
            "https://doi.org/10.1371/journal.pcbi.1004226",
            "bda074fd8f1ffa410040a52673e8f4a1c8866ea16cac0156a95716b0bdd29890",
        ),
    }
    for cite_key, (identifier, digest) in expected.items():
        assert custody.loc[cite_key, "doi_or_identifier"] == identifier
        assert custody.loc[cite_key, "local_sha256"] == digest
        assert custody.loc[cite_key, "distribution"] == "local custody only"

    tracked = [
        path.decode("utf-8")
        for path in subprocess.check_output(
            ["git", "-C", str(ROOT), "ls-files", "-z"]
        ).split(b"\0")
        if path
    ]
    source_texts = [
        path for path in tracked if path.lower().endswith((".nxml", ".html"))
    ]
    literature_pdfs = [
        path
        for path in tracked
        if path.lower().endswith(".pdf")
        and (
            path.startswith(("literature/", "review-literature/"))
            or "/repository/resources/pdfs/" in path
        )
    ]
    assert source_texts == []
    assert literature_pdfs == []
