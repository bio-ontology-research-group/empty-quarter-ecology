"""Numerical claim checks tying the ecology manuscript to canonical artifacts.

Current numerical assertions compare manuscript claims with their source
artifacts. Superseded analyses retain separate archival checks; their old
values are not required in the current paper.
"""

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from manuscript_paths import PAPER

ROOT = Path(__file__).resolve().parents[1]
MAIN = PAPER / "main.tex"
SUPPLEMENT = PAPER / "supplement.tex"
MANIFEST = PAPER / "figures" / "figure_manifest.tsv"
RESULTS = ROOT / "analysis/v3/results"


@pytest.fixture(scope="module")
def main_tex() -> str:
    return MAIN.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def supplement_tex() -> str:
    return SUPPLEMENT.read_text(encoding="utf-8")


def _without_value_math(text: str) -> str:
    """Normalize numeric LaTeX typesetting while retaining the claim wording."""
    return text.replace("$", "").replace("{,}", ",")


@pytest.fixture(scope="module")
def rain() -> dict:
    base = ROOT / "analysis/v3/rain_pulse_response"
    return {
        "decision": json.loads((base / "analysis_decision.json").read_text()),
        "cohort": pd.read_csv(base / "analysis_cohort.tsv", sep="\t"),
        "scan": pd.read_csv(base / "pulse_peak_scan.tsv", sep="\t"),
        "placebos": pd.read_csv(base / "future_rain_placebo_scan.tsv", sep="\t"),
    }


def test_reverse_primer_is_bakt_785r(main_tex):
    assert "Bakt\\_341F" in main_tex
    assert "Bakt\\_785R" in main_tex
    assert "GACTACHVGGGTATCTAATCC" in main_tex
    assert "CCTACGGGNGGCWGCAG" in main_tex
    assert "klindworth2013primers" in main_tex
    # The retired identity and its sequence must be gone.
    assert "806R" not in main_tex
    assert "GGACTACNVGGGTWTCTAAT" not in main_tex


def test_amplicon_filter_description_matches_the_canonical_branch(
    main_tex, supplement_tex
):
    build_cache = (
        ROOT / "analysis/v2/review/build_cache.py"
    ).read_text(encoding="utf-8")
    assert 'c.startswith("EB") or c.startswith("Negative")' in build_cache
    assert "depth >= 1000" in build_cache
    flat_main = _without_value_math(" ".join(main_tex.split()))
    flat_supplement = _without_value_math(" ".join(supplement_tex.split()))
    assert (
        "all 351,472 ASVs and 24 control profiles from the Trip~5 "
        "sequencing run (18 extraction blanks and 6 no-template PCR controls)"
        in flat_main
    )
    assert "profiles below 1,000 reads" in flat_main
    assert "applied no rare-feature or contaminant filter" in " ".join(
        flat_supplement.split()
    )
    assert "unfiltered canonical table was retained" in " ".join(
        supplement_tex.split()
    )
    assert "fewer than three samples" not in main_tex
    assert "threshold of 0.52" not in main_tex


def test_klindworth_entry_resolves(main_tex):
    bib = (PAPER / "sample.bib").read_text(encoding="utf-8")
    assert "@article{klindworth2013primers" in bib
    assert "10.1093/nar/gks808" in bib
    log = PAPER / "main.log"
    if log.exists():
        text = log.read_text(encoding="utf-8", errors="replace")
        assert "Citation" not in text or "undefined" not in text


def test_section_order_matches_isme_communications(main_tex):
    order = re.findall(r"^\\section\*\{([^}]+)\}", main_tex, re.M)
    assert order == [
        "Introduction",
        "Results",
        "Discussion",
        "Materials and Methods",
        # Added on Overleaf (Sep 2026); sits between Methods and the
        # bibliography as ISME Communications back matter.
        "Funding",
        "Competing interests",
    ]


def test_rainfall_denominators_are_reported_in_supplement(supplement_tex, rain):
    cohort = rain["cohort"]
    current = pd.read_csv(
        ROOT / "analysis/v3/rain_calendar_refit_20260909/analysis_cohort.tsv", sep="\t"
    )
    keys = ["Trip", "Site", "Type"]
    assert set(map(tuple, current[keys].to_numpy())) == set(map(tuple, cohort[keys].to_numpy()))
    assert current.loc[current.Trip == 4, "Date"].min() == "2024-08-29"
    assert current.loc[current.Trip == 4, "Date"].max() == "2024-08-31"
    supplement_flat = _without_value_math(" ".join(supplement_tex.split()))
    assert len(cohort) == 617
    assert cohort["Site"].nunique() == 60
    assert cohort["Type"].value_counts().to_dict() == {
        "Rhizosphere": 218,
        "Deep": 203,
        "Surface": 196,
    }
    assert "617 site--campaign--position groups" in supplement_flat
    assert "179 pairs" in supplement_flat
    assert "20 distinct rainfall series" in supplement_flat
    assert "Open-Meteo resolved 58 series" in supplement_flat
    assert "historical surface denominator" not in supplement_flat


def test_rainfall_family_extremes_match_the_artifact(main_tex, supplement_tex, rain):
    # Preserve the historical result test without requiring superseded values
    # in the active manuscript.
    decision = rain["decision"]
    assert decision["selected_peak_complete_days"] == 2.0
    assert decision["familywise_inference"]["conditional_lag_rotation_one_sided_p"] == pytest.approx(0.056)
    assert decision["familywise_inference"]["conditional_lag_rotation_two_sided_p"] == pytest.approx(0.07865)
    selected = rain["scan"].query("endpoint == 'richness_hurlbert_25000' and candidate_peak_days == 2.0").iloc[0]
    assert round(float(selected["estimate_per_mm_at_kernel_peak"]), 1) == 179.6
    current = ROOT / "analysis/v3/rain_calendar_refit_20260909"
    result = json.loads((current / "summary.json").read_text())
    intervals = pd.read_csv(current / "refitted_site_bootstrap_summary.tsv", sep="\t").set_index(["product", "quantity"])
    generated = (PAPER / "generated/rain_calendar_refit_tables.tex").read_text()
    for product, peak in (("nasa_power", 2.0), ("open_meteo", 1.0)):
        observed = result["observed"][product]
        assert observed["endpoint"] == "richness_hurlbert_25000"
        assert observed["peak_days"] == peak
        assert f"{observed['beta']:.1f}" in main_tex
        assert f"{observed['beta']:.1f}" in generated
        for column in ("lower_2_5", "upper_97_5"):
            assert f"{intervals.loc[(product, 'fixed_beta'), column]:.1f}" in generated
            assert f"{intervals.loc[(product, 'fixed_endpoint_selected_peak_days'), column]:.1f}" in generated
    assert result["year_orbit_size"] == 6
    assert result["exact_p"]["nasa_power_positive_max"] == pytest.approx(1 / 6)
    assert result["exact_p"]["joint_positive_max"] == pytest.approx(2 / 6)
    assert "$p=1/6$" in main_tex and "$p=2/6$" in main_tex
    assert "calendar-preserving" in supplement_tex
    for retired in ("$p=0.0560$", "$p=0.0787$", "179.6", "327.5"):
        assert retired not in main_tex + supplement_tex + generated


def test_campaign_interaction_reports_both_models(main_tex, supplement_tex):
    verdict = json.loads(
        (ROOT / "analysis/v3/depth_extraction/claim_verdict.json").read_text()
    )
    pair = verdict["campaign_by_position_interaction"]
    assert round(pair["unadjusted_wald_p"], 5) == 0.00831
    assert round(pair["depth_adjusted_wald_p"], 5) == 0.17476
    assert "$p=0.00831$" in supplement_tex or "Wald $p=0.00831$" in supplement_tex
    assert "0.17476" in supplement_tex
    assert "0.00831" not in main_tex and "0.17476" not in main_tex
    # The depth-adjusted model is null, so the unqualified interpretation
    # must not be reasserted.
    assert "campaign-dependent depth pattern, not a persistent refugium" not in main_tex


def test_rubisco_correlation_is_rounded_and_carries_its_p(main_tex, supplement_tex):
    summary = json.loads(
        (
            ROOT
            / "analysis/v3/measured_function_summary_results"
            / "encoded_function_summary.json"
        ).read_text()
    )
    marker = summary["rubisco_marker"]
    assert round(marker["across_sample_spearman"], 2) == 0.16
    assert round(marker["two_sided_p_value"], 3) == 0.067
    flat = " ".join(supplement_tex.split())
    assert "RuBisCO agreement was weaker ($r=0.16$, two-sided $p=0.067$)" in flat
    assert "two-sided $p=0.067$" in flat
    assert "RuBisCO $r=0.17$" not in flat
    assert "0.17" not in " ".join(main_tex.split()).split("RuBisCO")[-1][:80]


def test_paired_compartment_contrasts_match_the_artifact(
    main_tex, supplement_tex
):
    paired = pd.read_csv(RESULTS / "paired_compartment_effects.tsv", sep="\t")
    rows = paired[
        (paired["trip"].astype(str) == "all") & (paired["metric"] == "shannon")
    ]
    values = {
        row.comparison: (
            round(row.mean_difference, 3),
            round(row.ci_low, 3),
            round(row.ci_high, 3),
        )
        for row in rows.itertuples()
    }
    assert values["Rhizosphere-Surface"] == (-0.130, -0.321, 0.079)
    assert values["Rhizosphere-Deep"] == (-0.206, -0.398, -0.001)
    assert "$-0.130$" in supplement_tex and "$-0.206$" in main_tex
    assert "[$-0.321$,$0.079$]" in supplement_tex
    assert "[$-0.398$,$-0.001$]" in supplement_tex
    flat = _without_value_math(" ".join(main_tex.split()))
    assert "95\\,\\% site-bootstrap interval, -0.492 to -0.031" in flat
    assert "by 0.269 on average" in flat


def test_primary_and_evenness_shannon_bootstraps_are_identical():
    primary = pd.read_csv(
        RESULTS / "paired_compartment_effects.tsv", sep="\t"
    )
    primary = primary[
        (primary["trip"].astype(str) == "all")
        & (primary["metric"] == "shannon")
    ].set_index("comparison")
    decomposition = pd.read_csv(
        ROOT / "analysis/v3/evenness_decomposition/paired_contrasts.tsv",
        sep="\t",
    )
    decomposition = decomposition[
        decomposition["metric"] == "shannon"
    ].set_index("contrast")
    assert set(primary.index) == set(decomposition.index)
    assert set(primary["bootstrap_resamples"]) == {1_000_000}
    assert set(decomposition["bootstrap_resamples"]) == {1_000_000}
    assert set(primary["bootstrap_seed"]) == {20260723}
    assert set(decomposition["bootstrap_seed"]) == {20260723}
    for contrast in primary.index:
        assert (
            abs(
                primary.loc[contrast, "ci_low"]
                - decomposition.loc[contrast, "bootstrap_ci_low"]
            )
            < 1e-12
        )
        assert (
            abs(
                primary.loc[contrast, "ci_high"]
                - decomposition.loc[contrast, "bootstrap_ci_high"]
            )
            < 1e-12
        )


def test_no_ambiguous_negative_to_positive_ranges(main_tex, supplement_tex):
    pattern = re.compile(r"\$-[0-9.]+\$--\$")
    assert not pattern.search(main_tex)
    assert not pattern.search(supplement_tex)


def test_unimplemented_alpha_metrics_are_not_claimed(main_tex):
    alpha = pd.read_csv(
        ROOT / "analysis/v2/review/cache/alpha.tsv", sep="\t", index_col=0
    )
    assert "faith_pd" not in alpha.columns
    assert not any(column.startswith("hill") for column in alpha.columns)
    assert "Faith's" not in main_tex
    assert "Hill numbers" not in main_tex


def test_rarefaction_subset_is_disclosed_in_supplement(main_tex, supplement_tex):
    alpha = pd.read_csv(
        ROOT / "analysis/v2/review/cache/alpha.tsv", sep="\t", index_col=0
    )
    alpha = alpha[alpha["Trip"].between(1, 5)]
    below = int((alpha["depth"] < 25000).sum())
    core = alpha[alpha["Site"].astype(int) <= 60]
    below_core = int((core["depth"] < 25000).sum())
    assert (below, below_core) == (60, 55)
    flat_supplement = _without_value_math(supplement_tex)
    flat_main = _without_value_math(main_tex)
    assert "Sixty of the 1,237 profiles" in flat_supplement
    assert "55 of the 1,227" in flat_supplement
    assert "1,177" in flat_supplement
    assert "1,172" in flat_supplement
    assert "Sixty of the 1,237 profiles" not in flat_main


def test_primary_bootstrap_source_is_unambiguous(main_tex, supplement_tex):
    assert "analysis/v3/results/paired_compartment_effects.tsv" in supplement_tex
    flat_main = _without_value_math(" ".join(main_tex.split()))
    flat_supplement = _without_value_math(supplement_tex)
    assert "seed 20260723" in flat_main
    assert "1,000,000 whole-site bootstrap samples" in flat_main
    assert "normalized Shannon diversity" in main_tex
    assert "same 1,000,000-resample" in flat_supplement
    assert "their estimates, intervals and unadjusted $p$ values for the identical Shannon estimand agree" in " ".join(supplement_tex.split())


def test_pma_parent_tubes_and_equal_illumination_confirmed(main_tex, supplement_tex):
    combined = " ".join((main_tex + supplement_tex).split())
    for confirmed in (r"50\,\mu\mathrm{M}", "10-minute dark incubation", r"450\,\mathrm{nm}", "35 minutes"):
        assert confirmed in combined
    assert "two campsites" in combined
    assert "three" in combined and "site--compartment groups" in combined
    assert "three separately prepared slurries" in combined
    assert "one parent collection tube" in combined
    assert "Both treated and untreated aliquots" in combined
    assert "same 10-minute dark incubation and 35-minute light exposure" in combined
    assert "untreated-aliquot illumination unspecified" not in combined
    assert "biological replication unresolved" not in combined
    assert "executed replicate hierarchy remains unresolved" not in combined
    assert "descriptive" in combined
    endpoints = pd.read_csv(ROOT / "analysis/v3/pma_endpoint_results/pma_pair_endpoints.tsv", sep="\t")
    assert len(endpoints) == 9
    assert set(endpoints.pair_id.str[:-1]) == {"C1R", "C2R", "C2S"}
    assert set(endpoints.rarefaction_depth) == {123897}
    assert int((endpoints.rarefied_richness_difference_treated_minus_untreated < 0).sum()) == 8


def test_picrust_seed_and_multiplicity_provenance_is_explicit(
    main_tex, supplement_tex
):
    combined = _without_value_math(" ".join((main_tex + supplement_tex).split()))
    assert "separately executed Trips~1--4 and Trip~5 input sets" in combined
    assert "330,830 ASVs" in combined
    assert "combined canonical table used for the primary" in combined
    assert "351,472 ASVs" in combined
    assert "Multiplicity correction was applied separately to the stated" in combined
    for seed in ("20260723", "20260725", "20260728"):
        assert seed in _without_value_math(supplement_tex)
    assert "PMA endpoints are summarized descriptively within three groups" in combined
    assert "Encoded-function agreement used 9,999 whole-site" in combined
    assert "with seed 20260805" in combined


def test_confirmed_laboratory_volumes_and_rinsing_are_reported(main_tex, supplement_tex):
    combined = " ".join((main_tex + supplement_tex).split())
    for value in (r"4{,}000\,\mathrm{rpm}", "45$-second cycle",
                  r"8\,\mu\mathrm{L}", r"12\,\mu\mathrm{L}",
                  r"31.75\,\mu\mathrm{L}", r"28.75\,\mu\mathrm{L}",
                  r"24.75\,\mu\mathrm{L}", "$0.8$:$1$", "target molarity",
                  "room temperature", "Milli-Q", "Kimwipes"):
        assert value in combined
    assert "programme-to-library allocation before submission" not in combined
    assert "All Trip~4 pH measurements used field replicate~2" in combined
    assert "replicate~2 was requested" not in combined.lower()


def test_latest_confirmed_primary_pcr_and_ph_protocol(main_tex, supplement_tex):
    for text in (main_tex, supplement_tex):
        flat=" ".join(text.split())
        for fact in ("$35$ cycles", r"$95\,^{\circ}\mathrm{C}$ for $10\,\mathrm{s}$",
                     r"$65\,^{\circ}\mathrm{C}$ for $30\,\mathrm{s}$",
                     r"$75\,^{\circ}\mathrm{C}$", "13-620-112", "hula mixer",
                     "green stability indicator", "All Trip~4 pH measurements used field replicate~2"):
            assert fact in flat
        assert "$30$ cycles" not in flat
        assert "14~July, used replicate~1" not in flat
        assert "illumination unspecified" not in flat
    flat=" ".join(supplement_tex.split())
    for fact in (r"$10\times$ Mg$^{2+}$-containing buffer", r"$2.5\,\mathrm{mM}$ each",
                 r"$5\,\mathrm{U}/\mu\mathrm{L}$", "blank template",
                 r"V_i=V_{\mathrm{pool}}C_{\mathrm{pool}}/(N C_i)",
                 "indexing-PCR reaction mixture remains unspecified"):
        assert fact in flat


def test_final_ecology_cohort_is_reported_without_repair_history(supplement_tex):
    flat = _without_value_math(" ".join(supplement_tex.split()))
    assert "final canonical" in flat
    assert "ten biological profiles with fewer than 1,000 reads" in flat
    assert "1,237 quality-controlled biological profiles" in flat
    for process_history in (
        "multiset difference",
        "older \\texttt{feature-table.tsv}",
        "Source-column alignment was repaired",
        "not malformed records",
        "erroneous duplicate",
        "The replacement summed sequencing replicates",
    ):
        assert process_history not in supplement_tex


def test_reproducibility_boundary_does_not_promise_a_missing_status_table(
    main_tex,
):
    flat = " ".join(main_tex.split())
    assert "Supplementary Information records the current component" not in flat
    assert "The public repository provides the staged inputs and downstream analyses" in flat
    assert "outstanding archival and public-download checks" in flat


def test_figure_manifest_matches_the_current_rainfall_artifact():
    # The earlier figure_manifest.tsv describes the archived circular-lag
    # renderer. Current rainfall tables carry their own source custody.
    manifest = json.loads((PAPER / "generated/rain_calendar_refit_tables.manifest.json").read_text())
    assert manifest
    figure = PAPER / "figures/rain_calendar_refit.pdf"
    assert figure.read_bytes().startswith(b"%PDF")
    for relative, expected in manifest["outputs"].items():
        assert hashlib.sha256((PAPER / relative).read_bytes()).hexdigest() == expected
    assert r"{figures/rain_calendar_refit.pdf}" in SUPPLEMENT.read_text()
    source = ROOT / "analysis/v3/rain_calendar_refit_20260909"
    for item in ("summary.json", "calendar_year_orbit.tsv", "refitted_site_bootstrap_summary.tsv"):
        assert hashlib.sha256((source / item).read_bytes()).hexdigest() in json.dumps(manifest)


def test_xrf_geographic_scope_is_bounded(main_tex, supplement_tex):
    flat = " ".join(main_tex.split())
    assert "within-site composition model" in flat
    assert r"0.358\,\%" in flat and "$p=0.010$" in flat
    assert "residual variation" in flat
    assert "undocumented concentration units" in flat
    assert "site-level" in supplement_tex


def test_xrf_primary_status_sensitivity_matches_the_canonical_table(
    main_tex, supplement_tex
):
    sensitivity = pd.read_csv(
        ROOT / "analysis/xrf_audit/xrf_aggregation_sensitivity.tsv",
        sep="\t",
    )
    affected = sensitivity[
        (sensitivity["workflow"] == "lab_t5")
        & (sensitivity["groups_current_differs_from_primary"] > 0)
    ]
    assert len(affected) == 11
    minimum = float(affected["spearman_current_vs_primary"].min())
    assert 0.994 <= minimum < 0.995
    flat = _without_value_math(" ".join(supplement_tex.split()))
    assert "primary-status sensitivity for the 11 affected Trip~5" in flat
    assert "Spearman \\rho\\geq0.994" in flat
    assert "primary-status sensitivity" not in main_tex
    assert "maximum-positive sensitivity" not in flat


def test_all_supplementary_floats_are_cited(supplement_tex):
    labels = re.findall(r"\\label\{((?:tab|fig):[^}]+)\}", supplement_tex)
    assert labels
    for label in labels:
        assert f"\\ref{{{label}}}" in supplement_tex, label


def test_primary_aitchison_configuration_is_explicit(main_tex):
    verdict = json.loads(
        (
            ROOT / "analysis/v3/xrf_community_clr/claim_verdict.json"
        ).read_text()
    )
    assert verdict["input"]["primary_taxon_count"] == 200
    assert verdict["input"]["primary_zero_treatment"] == "pseudocount_0.5"
    flat = _without_value_math(" ".join(main_tex.split()))
    assert "selected the 200 most abundant genera among those detected in at least 20\\,\\% of groups" in flat
    assert "added 0.5 to every count to handle zeros" in flat


def test_phylogenetic_diagnostic_matches_the_artifact(main_tex, supplement_tex):
    summary = json.loads(
        (
            ROOT / "analysis/v3/phylo_signal_results/phylo_signal_summary.json"
        ).read_text()
    )
    diagnostics = summary["diagnostics"]
    assert round(diagnostics["short_distance_pearson_r"], 4) == 0.0262
    assert round(diagnostics["overall_mantel_r"], 4) == -0.0837
    assert round(diagnostics["closest_class_mean_niche_difference"], 3) == 0.808
    assert round(diagnostics["farthest_class_mean_niche_difference"], 3) == 0.515
    # The bounded prerequisite failure and its numbers are supplementary.
    assert "$r=0.0262$" not in main_tex
    assert "$r=0.0262$" in supplement_tex
    assert "$r=-0.0837$" in supplement_tex
    assert "0.808 versus 0.515" in _without_value_math(supplement_tex)
    assert "showed no positive association with niche difference" not in main_tex


def test_venue_structure_and_abstract_limit(main_tex, supplement_tex):
    keywords = re.search(r"\\keywords\{(.*?)\}", main_tex, re.S).group(1)
    count = len([k for k in keywords.replace("\n", " ").split(",") if k.strip()])
    assert 3 <= count <= 10
    assert main_tex.count("\\begin{figure}") <= 8
    assert "fig4_network_function.pdf" not in main_tex
    assert "fig4_network_function.pdf" not in supplement_tex
    bbl = PAPER / "main.bbl"
    if bbl.exists():
        assert bbl.read_text(encoding="utf-8").count("\\bibitem") <= 100

    texcount = shutil.which("texcount")
    assert texcount is not None, "texcount is required for venue word-limit validation"
    abstract = re.search(
        r"\\begin\{abstract\}(.*?)\\end\{abstract\}", main_tex, re.S
    ).group(1)
    counted = subprocess.run(
        [texcount, "-sum", "-q", "-"],
        input=abstract,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    abstract_words = int(re.search(r"Sum count: (\d+)", counted).group(1))
    # Keep copy-editing headroom below the formal 250-word ceiling.
    assert abstract_words <= 245


def test_rainfall_scope_is_consistent_across_the_whole_manuscript(main_tex, supplement_tex, rain):
    assert rain["decision"]["analysis_status"] == "temporally_localized_association_borderline"
    flat = _without_value_math(" ".join(main_tex.split()))
    si = _without_value_math(" ".join(supplement_tex.split()))
    assert "6 assignments" in flat and "2023, 2024 and 2025" in flat
    assert "exchangeability" in flat and "exchangeability" in si
    assert "refitting all nuisance coefficients in every draw" in flat
    assert "full endpoint--peak search separately within each product" in si
    assert "sampling window for event-triggered tests" in flat
    assert "rain caused" not in (flat + si).lower()
    assert "figS_campaign_rainfall.pdf" not in main_tex + supplement_tex
    assert "rain_calendar_refit.pdf" in supplement_tex


def test_every_main_text_cross_reference_resolves_in_main_text(main_tex):
    """Relocating a float must not leave a dangling \\ref behind."""
    labels = set(re.findall(r"\\label\{([^}]+)\}", main_tex))
    referenced = set(re.findall(r"\\ref\{([^}]+)\}", main_tex))
    assert referenced <= labels, f"dangling in main.tex: {referenced - labels}"
    assert "fig:network-function" not in labels
    assert "fig:network-function" not in referenced
    # Three evidence-bearing figures remain in the main paper.
    assert len(labels & {
        "fig:landscape",
        "fig:soil-position",
        "fig:function-controls",
    }) == 3
    assert main_tex.count("\\includegraphics") == 3
    assert main_tex.count("\\begin{table}") == 0


def test_relocated_detail_survives_in_the_supplement(main_tex, supplement_tex):
    """Every value removed from the main text must be findable elsewhere."""
    # The Supplement reports the positive-edge fractions as proportions
    # (0.6278) where the deleted main text used percentages (62.8%).
    network_values = ["0.05374", "0.0096", "0.6278", "0.05696", "180", "237"]
    function_values = [
        "0.82272",
        "0.13386",
        "0.08733",
        "0.13659",
        "0.16602",
        "28.5",
    ]
    phylogenetic_values = ["-0.0837", "0.808 versus 0.515"]
    normalized_supplement = _without_value_math(supplement_tex)
    for value in network_values + function_values + phylogenetic_values:
        assert value in normalized_supplement, f"{value} lost in relocation"
    # The detailed network and encoded-function methods left the main text.
    assert "Graphical Lasso" not in main_tex
    assert "intact genome annotation profiles among the fixed" not in main_tex
    assert "For KO-copy profiles, the observed median" not in main_tex
    # The main Methods deliberately keep a short shotgun-provenance paragraph
    # (added on Overleaf, Sep 2026) naming the companion metagenomic study
    # and its eggNOG-mapper annotation, without the encoded-function detail.
    flat_main = " ".join(main_tex.split())
    assert "companion metagenomic study of the same expeditions" in flat_main
    assert "eggNOG-mapper v2.1.12 annotation of the genome catalogue" in flat_main
    assert "CoverM genome-abundance profiles and KEGG Ortholog annotations" in flat_main
    # Provenance paragraph, trait-gene Methods and Software.
    assert main_tex.count("eggNOG-mapper") == 3
    # The Supplement carries the shotgun/genomic-potential methods instead.
    assert "990" in normalized_supplement and "2,000" in normalized_supplement


def test_network_calibration_diagnostic_is_named_correctly(supplement_tex):
    flat = " ".join(supplement_tex.split())
    assert "Expected false fraction" in flat
    assert "Expected false fractions summarize the permuted-taxon calibration" in flat
    assert "null-selection fraction" not in flat
    assert "0.05374" in flat


def test_companion_title_is_consistent(main_tex, supplement_tex):
    expected = "Landscape-scale bacterial biogeography across the Rub' al-Khali"
    flat_main = " ".join(main_tex.split())
    flat_supplement = " ".join(supplement_tex.split())
    assert "\\title{" + expected + "}" in flat_main
    assert expected in flat_supplement
    # The subtitle was dropped on 3 Sep 2026 at a co-author's suggestion.
    assert "reveals recurring" not in flat_main
    assert "reveals recurring" not in flat_supplement
    bibliography = (
        ROOT / "data-paper/sn-bibliography.bib"
    ).read_text(encoding="utf-8")
    assert expected.replace("Rub'", "{Rub'").replace("al-Khali", "al-Khali}") \
        in bibliography
    assert "compartment and genomic-potential structure" not in (
        main_tex + supplement_tex + bibliography
    )


def test_retired_ecology_tree_has_executable_build_guards():
    retired = ROOT / "ecology-paper"
    makefile = (retired / "Makefile").read_text(encoding="utf-8")
    latexmkrc = (retired / ".latexmkrc").read_text(encoding="utf-8")
    assert "ecology-paper is retired" in makefile
    assert "../empty-quarter-amplicon" in makefile
    assert "@false" in makefile
    assert "ecology-paper is retired" in latexmkrc


def test_abstract_interprets_the_supported_soil_position_result(main_tex):
    abstract = re.search(
        r"\\begin\{abstract\}(.*?)\\end\{abstract\}", main_tex, re.S
    ).group(1)
    flat = " ".join(abstract.split())
    dispersion = pd.read_csv(
        ROOT
        / "analysis/v3/compartment_composition/compartment_dispersion_results.tsv",
        sep="\t",
    )
    primary = dispersion[dispersion["analysis"] == "primary"].set_index(
        "contrast"
    )
    assert primary.loc["Rhizosphere-Surface", "q_primary_family"] < 0.05
    assert primary.loc["Rhizosphere-Deep", "q_primary_family"] > 0.05
    assert "root-adjacent communities differed from bulk soil in composition" in flat
    assert "lower normalized Shannon diversity" in flat
    assert "pseudo-$F" not in flat
    assert "partial $R" not in flat
