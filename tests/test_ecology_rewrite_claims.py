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
from manuscript_paths import PAPER


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
        for name in re.findall(r"^\\author(?:\[[0-9]+\])?\{([^{}]+)", main_tex, re.M)
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
            r"^\\author\[([0-9]+)\]\{([^{}]+)", main_tex, re.M
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
        assert ecology_by_name[name] == "1"
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

    data_source = _read_manuscript_source(ROOT / "data-paper/sn-article.tex")
    data_author_pairs = re.findall(
        r"^\\author\*?\[([0-9]+)\]\{\\fnm\{([^{}]+)\} \\sur\{([^{}]+)\}\}",
        data_source,
        re.M,
    )
    data_by_name = {
        f"{given} {surname}": affiliation
        for affiliation, given, surname in data_author_pairs
    }
    # The pinned data-paper clone lags the live list (it may still carry
    # Maxat Kulmanov, Hind Aldakhil at BORG, and older name spellings); only
    # the authors present in both are compared.
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
    assert data_by_name["Michel Dumontier"] == "7"
    assert data_by_name["Raik Gr\\\"unberg"] == "6"
    assert data_by_name["Tiannyu Wang"] == "5"
    assert data_by_name["Magnus Rueping"] == "5"
    assert "Bio-Ontology Research Group (BORG)" in data_source
    assert "Mathematical Sciences and Engineering (CEMSE) Division" in data_source
    assert "Physical Science and Engineering (PSE) Division" in data_source
    assert "Biological and Environmental Science and Engineering (BESE)" in data_source
    assert "Institute of Data Science, Department of Advanced Computing" in data_source


def test_active_manuscript_prose_follows_robert_forbidden_word_list(
    main_tex, supplement_tex
):
    forbidden = re.compile(
        r"\b(?:thus|fortunately|unfortunately|interestingly|surprisingly|"
        r"clearly|obviously|very|quite|really|basic|basically|comprehensive|"
        r"unique|uniquely|rigorous|robust|utili[sz](?:e|ed|es|ing)|"
        r"demonstrat(?:e|ed|es|ing))\b",
        re.I,
    )
    # "Robust covariance" is the technical sandwich-estimator term.
    prose = (main_tex + supplement_tex).replace("robust covariance", "sandwich covariance")
    assert forbidden.search(prose) is None


def test_main_source_contains_its_prose_and_has_no_tex_fragment_includes(
    main_tex, supplement_tex
):
    # The preamble legitimately defines the \todo macro (\newcommand) and the
    # arabtex \let lines; the document BODY must contain no prose fragment
    # includes.  \includegraphics and \bibliography stay allowed.
    preamble, body = main_tex.split(r"\begin{document}", 1)
    assert re.search(r"\\(?:input|include|subfile)\s*\{", body) is None
    assert re.search(r"\\(?:input|include|subfile)\s*\{", preamble) is None
    assert r"\newcommand" not in body
    assert r"\PHEcologyMethods" not in main_tex
    assert r"\PHEcologyResults" not in main_tex
    assert r"\PHEcologyDiscussion" not in main_tex
    assert re.search(r"\\PH[A-Za-z]+", main_tex) is None
    flat_main = _without_value_math(_flat(main_tex))
    for visible_prose in (
        "The archived workbook contained 767 pH measurements from all campaigns",
        "In total, 712 measurements passed",
        "Of 712 accepted measurements, 709 linked to ecology groups",
        "Each percentage uses the variance of its own residual response as the denominator",
    ):
        assert visible_prose in flat_main
    # Rows and columns are valid terms for actual tables/matrix operations.
    assert "702 profiles" not in flat_main
    assert "563" in flat_main and "560" in flat_main and "1,086" in flat_main


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


def test_control_method_explains_training_scope(main_tex, supplement_tex):
    flat_main = _without_value_math(_flat(main_tex))
    flat_supplement = _without_value_math(_flat(supplement_tex))
    assert "17 sequenced extraction blanks" in flat_main
    assert "could be linked by extraction day to 217 Trip~5 profiles" in flat_main
    assert "95 Trip~4 profiles" in flat_main and "7 candidate ASVs" in flat_main
    assert "they were kept separate because their biological batch links were incomplete" in flat_main
    assert "Positive standards were used only to assess" in flat_main
    assert "one extraction blank per extraction kit" in flat_supplement
    assert "maps EB1--EB17 to dates and 220 Trip~5 biological profiles" in flat_supplement
    assert "217 occur in the canonical table" in flat_supplement
    assert "six \\texttt{Negative}-labelled profiles" in flat_supplement
    assert "no-template PCR controls (water in place of template)" in flat_supplement
    assert "not sequenced" in flat_supplement


def test_geography_first_results_order_and_regional_novelty_opener(main_tex):
    flat = _flat(main_tex)
    introduction = flat.split(r"\section*{Introduction}", 1)[1].split(
        r"\section*{Results}", 1
    )[0]
    assert "world's largest continuous sand desert" in introduction
    assert "regional context for a repeated survey of the open Empty Quarter" in introduction
    # Naming policy: Rub' al-Khali is introduced once, then "Empty Quarter".
    assert "Sampling the same sites and soil compartments across campaigns" in introduction
    assert "The only direct microbial study" not in introduction
    assert "support more than 38\\,\\% of the world's population" not in introduction

    geography = flat.index(
        r"\subsection*{Bacterial communities change across the landscape}"
    )
    paired = flat.index(
        r"\subsection*{Compartment differences in composition and normalized Shannon diversity}"
    )
    environment = flat.index(
        r"\subsection*{Climate gradients track bacterial diversity and composition}"
    )
    function = flat.index(
        r"\subsection*{Predicted metabolic pathways follow geography and compartment}"
    )
    assert geography < paired < environment < function
    assert (
        r"\subsection*{Relic-DNA and low-biomass controls preserve the main patterns}"
        not in flat
    )
    controls = flat.index(r"\subsection*{Assay controls and sensitivity analyses}")
    methods = flat.index(r"\section*{Materials and Methods}")
    assert methods < controls
    headings = re.findall(r"\\subsection\*\{([^}]+)\}", flat)
    assert all("baseline" not in heading.lower() for heading in headings)
    assert r"\subsection*{No " not in flat


def test_abstract_leads_with_science_and_keeps_resource_subordinate(main_tex):
    abstract = re.search(
        r"\\begin\{abstract\}(.*?)\\end\{abstract\}", main_tex, re.S
    ).group(1)
    flat = _flat(abstract)
    assert "distant sites contained more different communities" in flat
    assert "taxon replacement accounted for most incidence dissimilarity" in flat
    assert "root-adjacent communities differed from bulk soil" in flat
    assert "regional biological baseline" in flat
    assert "knowledge graph" not in flat.lower()
    assert flat.index("distant sites") < flat.index("linked resource")


def test_pH_attenuation_is_bounded_and_negative_diagnostics_are_supplementary(
    main_tex, supplement_tex
):
    flat = _without_value_math(_flat(main_tex))
    flat_supplement = _without_value_math(_flat(supplement_tex))
    summary = _json("analysis/v3/ph_group_linkage_20260909/summary.json")
    counts = summary["counts"]
    assert counts["group_linked_ph_measurements"] == 709
    assert counts["site_campaign_position_groups"] == 563
    assert counts["composition_groups"] == 560
    geo = summary["geographic_same_cohort"]
    for key in ("without_ph_adjustment_partial_r2", "with_ph_adjustment_partial_r2"):
        assert f"{100 * geo[key]:.1f}" in flat
    assert "Each percentage uses the variance of its own residual response as the denominator" in flat
    assert "no stable direct composition association" not in flat
    assert "group-linked" in flat_supplement
    assert "fig:environment-spatial" not in main_tex


def test_paired_composition_claims_match_the_canonical_verdict(
    main_tex, supplement_tex
):
    verdict = _json("analysis/v3/compartment_composition/claim_verdict.json")
    omnibus = verdict["omnibus"]
    contrasts = verdict["contrasts"]

    assert omnibus["n_sites"] == 60
    assert omnibus["n_blocks"] == 170
    assert round(omnibus["pseudo_f"], 2) == 17.28
    assert omnibus["permutation_p"] == 0.001
    for text in (main_tex, supplement_tex):
        assert "170 complete" in _without_value_math(text)
        assert "pseudo-$F" in text
    assert f"{omnibus['pseudo_f']:.2f}" in main_tex
    assert f"{omnibus['pseudo_f']:.4f}" in supplement_tex

    expected_displacements = {
        "Deep-Surface": "0.331",
        "Rhizosphere-Surface": "0.433",
        "Rhizosphere-Deep": "0.478",
    }
    for contrast, printed in expected_displacements.items():
        result = contrasts[contrast]
        assert f"{result['standardized_displacement']:.3f}" == printed
        assert result["permutation_p"] == 0.001
        assert result["q_within_primary_family"] == 0.001
        assert printed in supplement_tex
    assert "All three comparisons remained supported after correction for three tests" in _without_value_math(_flat(main_tex))
    assert "fig2_soil_position.pdf" in main_tex

    assert "sign-flipped the mean within-block CLR difference vector" in supplement_tex
    assert "centrally symmetric paired-difference null" in _flat(supplement_tex)
    assert "root distance" in main_tex


def test_evenness_decomposition_is_numerically_and_semantically_bounded(
    main_tex, supplement_tex
):
    verdict = _json("analysis/v3/evenness_decomposition/claim_verdict.json")
    results = verdict["primary_results"]
    root_surface = results["Rhizosphere-Surface"]["evenness_sensitivity"]
    root_shallow = results["Rhizosphere-Deep"]["evenness_sensitivity"]

    assert f"{root_surface['mean_difference']:.5f}" == "-0.03148"
    assert f"{root_shallow['mean_difference']:.5f}" == "-0.04398"
    for value in ("-0.03148", "-0.04398"):
        assert value not in main_tex
        assert value in supplement_tex

    assert "$H/\\log(E[S_{25k}])$" in main_tex
    assert "$H/\\log(E[S_{25k}])$" in supplement_tex
    assert "custom mixed-depth index" in _flat(main_tex)
    assert "custom normalized index" in _flat(supplement_tex)
    assert "An additional analysis of normalized Shannon diversity" in _flat(main_tex)
    assert "617 of the 633" not in main_tex
    assert "617 of 633" in _without_value_math(supplement_tex)
    assert verdict["input"]["blocks_with_evenness_sensitivity"] == 617


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


def test_moran_claim_is_bounded_to_the_tested_neighbourhood_scale(
    main_tex, supplement_tex
):
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
    flat_supplement = _flat(supplement_tex)
    for current in ("$0.1101$", "$-0.0098$", "$p=0.243$", "$p=0.050$"):
        assert current in flat_supplement
    for stale in (
        "0.1104",
        "-0.0097",
        "$p=0.241$ at",
        "through $k=8$",
        "earlier coordinate version",
        "629",
    ):
        assert stale not in flat_supplement

    # Current inference propagates corrected coordinates through specified
    # row covariances.
    current = _json("analysis/v3/spatial_covariance_sensitivity_20260909/summary.json")
    assert current["model_count"] == 73
    assert current["original_descriptive_r2"] == pytest.approx(0.40073421629965933)
    assert current["family_p_min"] == 0.0001 and current["family_p_max"] == 1.0
    assert "73" in main_tex and "0.0001" in main_tex and "1.0000" in main_tex
    assert "specified spatial covariance" in _flat(main_tex)
    assert "Moran's $I=0.06464$" in supplement_tex
    assert "neighbour-count diagnostics" in _flat(supplement_tex)


def test_landscape_figure_contains_six_evidence_bearing_panels(main_tex):
    figure_script = (
        ROOT / "analysis/v3/make_submission_figures.py"
    ).read_text(encoding="utf-8")
    study_figure = figure_script.split(
        "def make_landscape_figure", 1
    )[1].split("def make_soil_position_figure", 1)[0]

    assert "plt.subplots(3, 2" in study_figure
    for content in (
        "Repeated 60-site desert transect",
        "Communities diverge with distance",
        "Climate–diversity associations",
        "Climate-associated genera",
    ):
        assert content in study_figure
    assert "Explicit analysis cohorts" not in study_figure
    assert "fig:overview-revised" not in main_tex
    assert "(c) The samples available" not in main_tex


def test_xrf_non_detection_is_not_written_as_evidence_of_absence(
    main_tex, supplement_tex
):
    flat_main = _flat(main_tex)
    flat_supplement = _flat(supplement_tex)
    assert "did not track within-site Shannon diversity" not in flat_main
    assert r"0.358\,\%" in flat_main and "within-site composition model" in flat_main
    assert "$p=0.990$" in flat_supplement and "Shannon" in flat_supplement
    assert "$p=0.990$" not in flat_main
    assert "Its adjusted association with Shannon diversity was null" not in flat_main
    assert "It had no conditional association with Shannon diversity" not in flat_main
    assert "no adjusted Shannon association;" not in flat_supplement


def test_depth_adjusted_claims_match_the_canonical_verdict(
    main_tex, supplement_tex
):
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
    assert "sequencing-depth and extraction-method adjustments" in _flat(main_tex)
    for value in ("0.00831", "0.17476", "$-0.273$"):
        assert value in supplement_tex
    assert "$+0.103$" in supplement_tex
    assert "laboratory batch were incompletely recorded" in _flat(main_tex)
    assert sensitivity_dependent["status"] == "sensitivity_dependent"
    assert sensitivity_dependent["direction_stable_across_models"] is True
    assert sum(
        sensitivity_dependent["interval_excludes_zero_by_model"].values()
    ) == 3
    assert "root-adjacent and surface soil" in _flat(main_tex)
    flat_supplement = _flat(supplement_tex)
    assert "root-adjacent minus surface" in flat_supplement or (
        "root-adjacent--surface" in flat_supplement
    )
    assert "three of the six" in flat_supplement or "only three fits" in flat_supplement
    assert "Extraction kit was recorded" in flat_supplement and "complete cases" in flat_supplement
    assert "expected richness retained a depth-adjusted interaction" not in (
        main_tex + supplement_tex
    )


def test_consolidation_removes_untraceable_between_site_xrf_numbers(main_tex):
    flat = _flat(main_tex)
    assert "PC1 and Shannon diversity had $\\rho=-0.68$" not in flat
    assert "partial $\\rho=-0.30$" not in flat
    assert "block size increased from 3 to 20 sites" not in flat
    assert r"0.358\,\%" in flat and "within-site composition model" in flat


def test_campaign_omission_is_bounded_as_an_influence_analysis(
    main_tex, supplement_tex
):
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
    main_flat = _flat(main_tex)
    assert "leave-one-campaign-out refits" in main_flat
    flat_supplement = _flat(supplement_tex)
    assert "These influence analyses omit one campaign at a time" in flat_supplement
    for value in ("25.2", "2.5", "26.8", "27.9", "17.5"):
        assert value in flat_supplement


def test_assay_aware_control_filter_is_bounded_and_headlines_are_stable(
    main_tex, supplement_tex
):
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

    assert "351,472" in _without_value_math(main_tex)
    for value in ("351", "351,472", "217"):
        assert value in _without_value_math(supplement_tex)
    for text in (main_tex, supplement_tex):
        assert "14,822" not in text
        assert "2.8684" not in text
    # The 25-verdict source above remains an archived sensitivity. Current
    # claims additionally use corrected coordinates and calibrated tests.
    assert "25 tracked metrics" in _without_value_math(_flat(supplement_tex))
    assert "312 filtered profiles" in _without_value_math(_flat(supplement_tex))
    assert "0.0256" in supplement_tex and "0.0312" in supplement_tex
    assert "candidate ASVs were removed only from profiles linked to their extraction blanks" in _flat(main_tex)
    flat = _flat(main_tex)
    assert "unfiltered biological table was the primary analysis input" in flat
    flat_supplement = _without_value_math(_flat(supplement_tex))
    assert "one extraction blank per extraction kit" in flat_supplement
    assert "Positive controls were excluded from training" in flat_supplement
    assert "Trip~5 also used D6300" in flat_supplement


def test_functional_top_k_sensitivity_matches_all_three_canonical_tables(
    supplement_tex,
):
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

    for value in observed:
        assert value in supplement_tex
    # The prose rounds the first null median to 0.0670.
    for value in ("0.0670", "0.0873", "0.1174"):
        assert value in supplement_tex
    assert "$p_{\\rm U}=0.001$ throughout" in supplement_tex


def test_distance_decay_is_surfaced_and_uses_whole_site_permutations(main_tex, supplement_tex):
    current = _json("analysis/v3/distance_decay_turnover/claim_verdict.json")
    assert current["site_pairs"] == 1770
    assert current["matched_sites"] == 60
    assert current["standardised_depth"] == 12865
    assert current["permutations"] == 9999
    assert current["omnibus_p"] == 0.0051
    slopes = pd.read_csv(ROOT / "analysis/v3/distance_decay_turnover/distance_decay_slopes.tsv", sep="\t")
    flat = _without_value_math(_flat(supplement_tex))
    main = _without_value_math(_flat(main_tex))
    assert "1,770 geographic pairs" in flat
    assert "Whole-site permutations were applied simultaneously" in flat
    assert "site labels defined the permutation unit" in flat
    for row in slopes[slopes.family == "aitchison"].itertuples():
        assert f"{row.slope_per_100km:.3f}" in flat
        assert f"{row.slope_per_100km:.3f}" in main
        assert row.two_sided_p == 0.0001
    assert "p=0.0051" in flat and "p=0.0051" in main
    contrast = slopes.query("family == 'contrast' and response == 'Deep-Rhizosphere'").iloc[0]
    assert contrast.max_t_adjusted_p == 0.0039
    assert r"p_{\mathrm{adj}}=0.0039" in flat
    for row in slopes[slopes.family.isin(["simpson_turnover", "nestedness"])].itertuples():
        assert f"{row.slope_per_100km:.4f}" in flat
    assert r"69--75\,\% of Sørensen dissimilarity" in flat
    assert "Replacement supplied most mean dissimilarity" in flat
    assert "randomly reassigned location labels to whole sites" in main
    assert "keeping all 1,770 pairwise distances attached to their site identities" in main


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
