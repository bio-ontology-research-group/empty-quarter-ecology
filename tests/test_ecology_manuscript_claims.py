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
from manuscript_paths import DATA_REPO, PAPER, locked_data_bytes
import manuscript_text as mt
from manuscript_text import fmt, has_number, pct

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


def test_amplicon_filter_description_matches_the_canonical_branch():
    build_cache = (
        ROOT / "analysis/v2/review/build_cache.py"
    ).read_text(encoding="utf-8")
    assert 'c.startswith("EB") or c.startswith("Negative")' in build_cache
    assert "depth >= 1000" in build_cache
    audit = json.loads((ROOT / "analysis/v3/control_audit/summary.json").read_text())
    assert audit["canonical_profiles"] == 1271
    assert audit["canonical_features"] == 351472
    alpha = pd.read_csv(ROOT / "analysis/v2/review/cache/alpha.tsv", sep="\t", index_col=0)
    assert len(alpha) == 1237 and (alpha["depth"] >= 1000).all()
    # The 1,271 canonical profiles include the 24 assay controls:
    # 1,271 - 24 controls - 10 low-depth biological profiles = 1,237.
    assert audit["canonical_profiles"] - 24 - 10 == len(alpha)

    reads = mt.section("main", "Read processing", "Relic-DNA removal experiment")
    for number in ("1,271", "351,472", "24", "10", "1,000", "1,237"):
        assert has_number(reads, number), number
    # 1,271 counts every profile, controls included; it must not be presented
    # as 1,271 biological sample profiles plus 24 controls.
    assert "1,271 sample profiles" not in reads
    controls = mt.text("supplement")
    assert "24 assay-control profiles" in controls
    assert "18 extraction blanks and six PCR no-template controls" in controls
    # No rare-feature or contaminant filter was applied to the primary table.
    assert "retained all ASVs in downstream analyses" in mt.text("main")
    main = mt.text("main")
    assert "fewer than three samples" not in main
    assert "threshold of 0.52" not in main


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
        "Data availability",
        # Added on Overleaf (Sep 2026); sits between Methods and the
        # bibliography as ISME Communications back matter.
        "Funding",
        "Competing interests",
    ]


def test_rainfall_denominators_are_reported_in_supplement(rain):
    cohort = rain["cohort"]
    current = pd.read_csv(
        ROOT / "analysis/v3/rain_calendar_refit_20260909/analysis_cohort.tsv", sep="\t"
    )
    keys = ["Trip", "Site", "Type"]
    assert set(map(tuple, current[keys].to_numpy())) == set(map(tuple, cohort[keys].to_numpy()))
    assert current.loc[current.Trip == 4, "Date"].min() == "2024-08-29"
    assert current.loc[current.Trip == 4, "Date"].max() == "2024-08-31"
    assert len(cohort) == 617
    assert cohort["Site"].nunique() == 60
    assert cohort["Type"].value_counts().to_dict() == {
        "Rhizosphere": 218,
        "Deep": 203,
        "Surface": 196,
    }
    kernel = mt.section("supplement", "Rise-and-decay distributed lag kernel", "Complete calendar-year permutation calibration")
    assert "617 site × campaign × compartment combinations" in kernel
    assert "60 sites into 20 distinct rainfall series" in kernel
    assert "58 distinct series" in kernel
    assert "historical surface denominator" not in mt.text("supplement")


def test_rainfall_family_extremes_match_the_artifact(rain):
    # The archived circular-lag analysis keeps its own numerical record; its
    # superseded p values and estimates must not reappear in current text.
    decision = rain["decision"]
    assert decision["selected_peak_complete_days"] == 2.0
    # Archived circular-lag analysis, regenerated on corrected Site 52 inputs
    # (was p=0.056/0.07865 and 179.6 taxa per mm).
    assert decision["familywise_inference"]["conditional_lag_rotation_one_sided_p"] == pytest.approx(0.0545)
    assert decision["familywise_inference"]["conditional_lag_rotation_two_sided_p"] == pytest.approx(0.07705)
    selected = rain["scan"].query("endpoint == 'richness_hurlbert_25000' and candidate_peak_days == 2.0").iloc[0]
    assert round(float(selected["estimate_per_mm_at_kernel_peak"]), 1) == 179.4

    current = ROOT / "analysis/v3/rain_calendar_refit_20260909"
    result = json.loads((current / "summary.json").read_text())
    intervals = pd.read_csv(current / "refitted_site_bootstrap_summary.tsv", sep="\t").set_index(["product", "quantity"])
    generated = (PAPER / "generated/rain_calendar_refit_tables.tex").read_text()
    main_table = mt.section("main", r"\labeltab:rainfall", r"\endtable")
    main_text = mt.section("main", "Exploratory models suggest an early richness response", r"\begintable")
    supplement = mt.section("supplement", "Rise-and-decay distributed lag kernel", "Complete calendar-year permutation calibration")
    labels = {"nasa_power": "NASA POWER", "open_meteo": "Open-Meteo"}
    for product, peak in (("nasa_power", 2.0), ("open_meteo", 1.0)):
        observed = result["observed"][product]
        assert observed["endpoint"] == "richness_hurlbert_25000"
        assert observed["peak_days"] == peak
        beta = intervals.loc[(product, "fixed_beta")]
        r2 = intervals.loc[(product, "fixed_partial_r2")]
        timing = intervals.loc[(product, "fixed_endpoint_selected_peak_days")]
        selected = result["bootstrap_endpoint_frequencies"][product]["richness_hurlbert_25000"]
        # Supplement: one decimal, as in the generated table.
        for value in (observed["beta"], beta["lower_2_5"], beta["upper_97_5"]):
            assert f"{value:.1f}" in generated
            assert has_number(supplement, fmt(value, 1)), (product, value)
        assert has_number(supplement, pct(selected, 1))
        # Main table: integer effects, percent variance, timing and selection.
        row = main_table.split(labels[product], 1)[1].split(r"\\", 1)[0]
        assert f"{fmt(peak, 1)} [{fmt(timing['lower_2_5'], 1)}, {fmt(timing['upper_97_5'], 1)}]" in row
        assert f"{fmt(observed['beta'], 0)} [{fmt(beta['lower_2_5'], 0)}, {fmt(beta['upper_97_5'], 0)}]" in row
        assert f"{pct(observed['partial_r2'])} [{pct(r2['lower_2_5'])}, {pct(r2['upper_97_5'])}]" in row
        assert has_number(row, pct(selected, 0))
        assert has_number(main_text, fmt(observed["beta"], 0))
        assert has_number(main_text, pct(selected, 0))
        # The selected timing stays within the first week only if its upper
        # resampling bound is below seven days.
        assert timing["upper_97_5"] < 7
    assert "within the first week" in main_text
    assert result["year_orbit_size"] == 6
    assert result["exact_p"]["nasa_power_positive_max"] == pytest.approx(1 / 6)
    assert result["exact_p"]["open_meteo_positive_max"] == pytest.approx(2 / 6)
    assert result["exact_p"]["joint_positive_max"] == pytest.approx(2 / 6)
    assert "p=1/6" in main_text and "p=2/6" in main_text
    assert "1/6 & 2/6 & 2/6" in mt.text("supplement")
    combined = mt.combined() + " " + mt.normalize(generated)
    for retired in ("p=0.0560", "p=0.0545", "p=0.0787", "p=0.0771", "179.6", "327.5"):
        assert retired not in combined
    rainfall = main_text + " " + mt.section("supplement", "Short-term rainfall", "Beta diversity (between-sample)")
    assert not re.search(r"p ?= ?0\.05(?:6|45)\b", rainfall)


def test_campaign_interaction_reports_both_models():
    verdict = json.loads(
        (ROOT / "analysis/v3/depth_extraction/claim_verdict.json").read_text()
    )
    pair = verdict["campaign_by_position_interaction"]
    assert round(pair["unadjusted_wald_p"], 5) == 0.00831
    assert round(pair["depth_adjusted_wald_p"], 5) == 0.17476
    confounds = mt.section("supplement", "Potential confounds", "Climate associations")
    # Both models are reported at the supplement's three-decimal precision.
    assert has_number(confounds, fmt(pair["unadjusted_wald_p"], 3))
    assert has_number(confounds, fmt(pair["depth_adjusted_wald_p"], 3))
    assert "before adjustment for read depth" in confounds and "but not after" in confounds
    # The adjustment is bounded as a sensitivity analysis, not a batch control.
    assert "sensitivity analysis rather than definitive controls" in confounds
    main = mt.text("main")
    for value in ("0.00831", "0.17476", "0.008", "0.175"):
        assert not has_number(main, value)
    assert "campaign-dependent depth pattern, not a persistent refugium" not in main


def test_rubisco_correlation_is_rounded_and_carries_its_p():
    # The archived encoded-function RuBisCO marker (r=0.16) is a different
    # analysis from the trait-gene agreement reported now (rho=0.33).
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
    agreement = pd.read_csv(
        ROOT / "analysis/v3/trait_genes_20260909/source_agreement.tsv", sep="\t"
    ).set_index("trait")
    current = fmt(agreement.loc["rbcL_RuBisCO", "spearman_rho_genome_vs_picrust"], 2)
    assert current == "0.33"
    combined = mt.combined()
    agreement_tables = (
        mt.section("main", r"\labeltab:markers", r"\endtable"),
        mt.section("supplement", r"\labeltab:traits_picrust", r"\endtable"),
    )
    for table in agreement_tables:
        row = table.split("RuBisCO (rbcL/cbbL)", 1)[1].split(r"\\", 1)[0]
        assert row.strip().endswith(current), row
    # The superseded encoded-function value is not attributed to RuBisCO.
    assert not re.search(r"RuBisCO[^.]{0,80}\b0\.1[67]\b", combined)


def test_paired_compartment_contrasts_match_the_artifact():
    # Archived marginal estimate (now the supplementary alternative).
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

    estimands = _alpha_estimands()
    matched = "campaign_matched"
    labels = {
        "Rhizosphere-Surface": ("Root-adjacent - Surface", "Root-adjacent - surface"),
        "Rhizosphere-Deep": ("Root-adjacent - Subsurface", "Root-adjacent - shallow-subsurface"),
        "Deep-Surface": ("Subsurface - Surface", "Shallow-subsurface - surface"),
    }
    # Supplementary Table S (campaign-matched, common profiles): mean, CI and q.
    table_s = mt.section("supplement", r"\labeltab:matched-alpha", r"\endtable")
    blocks = {
        "shannon": table_s.split("Shannon diversity (H)", 1)[1].split("Rarefied richness", 1)[0],
        "richness_rare": table_s.split("Rarefied richness", 1)[1].split("Normalized Shannon", 1)[0],
        "normalized_shannon": table_s.split("Normalized Shannon", 1)[1],
    }
    decimals = {"shannon": 3, "richness_rare": 1, "normalized_shannon": 3}
    for metric, block in blocks.items():
        for contrast, (s_label, _) in labels.items():
            row = estimands.loc[("common_profiles", metric, matched, contrast)]
            line = block.split(s_label, 1)[1].split(r"\\", 1)[0]
            d = decimals[metric]
            assert has_number(line, fmt(row["mean"], d)), (metric, contrast, line)
            for bound in ("ci_low", "ci_high"):
                printed = {fmt(row[bound], d)} | ({fmt(row[bound], 4)} if abs(row[bound]) < 0.01 else set())
                assert any(has_number(line, value) for value in printed), (metric, contrast, bound, line)
            assert line.strip().endswith(_q_text(row["q_three_contrasts"])), (metric, contrast, line, row["q_three_contrasts"])
    # Main Table 2: two decimals for Shannon, integers for richness, three
    # decimals for normalized Shannon; stars mark q<0.05.
    main_table = mt.section("main", r"\labeltab:compartments", r"\endtable")
    main_decimals = {"shannon": 2, "richness_rare": 0, "normalized_shannon": 3}
    for contrast, (_, m_label) in labels.items():
        line = main_table.split(m_label + " &", 1)[1].split(r"\\", 1)[0]
        cells = [cell.strip() for cell in line.split("&")]
        for metric, cell in zip(("shannon", "richness_rare", "normalized_shannon"), cells[2:5]):
            row = estimands.loc[("common_profiles", metric, matched, contrast)]
            star = "*" if row["q_three_contrasts"] < 0.05 else ""
            assert cell == fmt(row["mean"], main_decimals[metric]) + star, (contrast, metric, cell)
    # Main text: the root-adjacent minus shallow-subsurface Shannon contrast.
    root_deep = estimands.loc[("common_profiles", "shannon", matched, "Rhizosphere-Deep")]
    text = mt.section("main", "Alpha diversity was lowest in root-adjacent soil", r"\begintable")
    assert (
        f"mean difference {fmt(root_deep['mean'], 2)}, 95% CI "
        f"{fmt(root_deep['ci_low'], 2)} to {fmt(root_deep['ci_high'], 2)}"
    ) in text
    for contrast in ("Rhizosphere-Surface", "Rhizosphere-Deep"):
        richness = estimands.loc[("common_profiles", "richness_rare", matched, contrast)]
        assert has_number(text, fmt(-richness["mean"], 0))
    # Marginal (endpoint-specific) richness q values quoted in the supplement.
    marginal = mt.section("supplement", "Marginal and paired compartment contrasts", r"\begintable")
    quoted = [
        fmt(estimands.loc[("endpoint_specific", "richness_rare", "marginal", c), "q_three_contrasts"], 2)
        for c in ("Rhizosphere-Surface", "Rhizosphere-Deep", "Deep-Surface")
    ]
    assert f"q = {', '.join(quoted)}" in marginal


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


def test_rarefaction_subset_is_disclosed_in_supplement():
    alpha = pd.read_csv(
        ROOT / "analysis/v2/review/cache/alpha.tsv", sep="\t", index_col=0
    )
    alpha = alpha[alpha["Trip"].between(1, 5)]
    below = int((alpha["depth"] < 25000).sum())
    core = alpha[alpha["Site"].astype(int) <= 60]
    below_core = int((core["depth"] < 25000).sum())
    assert (below, below_core) == (60, 55)
    retained = len(core) - below_core
    assert (len(core), retained) == (1227, 1172)
    indices = mt.section("supplement", "Diversity indices and rarefaction", "Marginal and paired compartment contrasts")
    assert has_number(indices, "55") and has_number(indices, "1,172")
    assert "below the selected read depth" in indices
    # Main Methods: the retained share and the reason for the 25,000-read depth.
    methods = mt.section("main", "Alpha diversity (within-sample)", "Compartment contrasts")
    assert int(100 * retained / len(core)) == 95
    assert "retain 95% of the samples" in methods
    assert has_number(methods, "25,000") and "5th percentile of library size" in methods
    assert has_number(methods, f"{int(alpha['depth'].median()):,}")
    assert "omitted from richness comparisons" in methods


def test_primary_bootstrap_source_is_unambiguous():
    manifest = json.loads(
        (ROOT / "analysis/v3/paired_alpha_sensitivity_20260909/manifest.json").read_text()
    )
    assert manifest["bootstrap_resamples"] == 1_000_000
    assert manifest["seed"] == 20260723
    assert manifest["bootstrap_unit"] == "site"
    assert manifest["source_sha256"] == hashlib.sha256(
        (ROOT / "analysis/v2/review/cache/alpha.tsv").read_bytes()
    ).hexdigest()
    table_s = mt.section("supplement", "Campaign-matched compartment differences", r"\endtable")
    assert "site-bootstrap intervals (10^6 draws)" in table_s
    assert "Wilcoxon signed-rank" in table_s
    caption = mt.section("main", "Pairwise differences between compartments", r"\labeltab:compartments")
    assert "common cohort with all three diversity measures" in caption
    assert "matched site–campaign blocks" in caption
    # Normalized Shannon is the reported evenness-sensitive endpoint.
    assert "normalized Shannon diversity" in mt.text("main")


def test_pma_parent_tubes_and_equal_illumination_confirmed():
    endpoints = pd.read_csv(ROOT / "analysis/v3/pma_endpoint_results/pma_pair_endpoints.tsv", sep="\t")
    assert len(endpoints) == 9
    assert set(endpoints.pair_id.str[:-1]) == {"C1R", "C2R", "C2S"}
    assert set(endpoints.rarefaction_depth) == {123897}
    lower = int((endpoints.rarefied_richness_difference_treated_minus_untreated < 0).sum())
    assert lower == 8
    summary = json.loads((ROOT / "analysis/v3/pma_endpoint_results/pma_summary.json").read_text())
    assert summary["provenance"]["seed"] == 20260805
    assert summary["rarefaction"]["depth"] == 123897

    methods = mt.section("main", "Relic-DNA removal experiment", "Environmental measurements")
    assert "three separate soil slurries from each of three parent soil tubes" in methods
    assert "two T5 campsites" in methods
    assert "50\\muM" in methods and "450nm" in methods
    assert "Both aliquots underwent a 10 min dark incubation" in methods
    assert "35 min" in methods
    assert has_number(methods, "123,897") and "minimum depth among the 18 libraries" in methods
    supplement = mt.section("supplement", "Relic DNA experiment", "Microbial community diversity")
    assert "both treated and untreated aliquots" in supplement
    assert "nine paired comparisons grouped within three parent tubes" in supplement
    # Results report the three group means and the eight-of-nine count.
    results = mt.section("main", "PMA-accessible DNA contributed to measured richness", r"\section*Discussion")
    assert "eight of the nine pairs" in results
    means = endpoints.assign(group=endpoints.pair_id.str[:-1]).groupby("group")
    richness = means["rarefied_richness_difference_treated_minus_untreated"].mean()
    shannon = means["shannon_difference_treated_minus_untreated"].mean()
    for group in ("C1R", "C2R", "C2S"):
        assert has_number(results, fmt(-richness[group], 0)), group
        value = fmt(shannon[group], 2)
        assert (("+" + value) if shannon[group] > 0 else value) in results, group
    table = mt.section("supplement", r"\labeltab:relicresults", r"\endtable")
    for group in ("C1R", "C2R", "C2S"):
        assert has_number(table, fmt(richness[group], 2))
        assert has_number(table, fmt(shannon[group], 3))
    # The design supports a descriptive group summary only.
    assert "describe the PMA response of those three parent samples" in results
    assert "descriptive comparisons" in mt.text("supplement")
    combined = mt.combined()
    for unresolved in ("untreated-aliquot illumination unspecified", "biological replication unresolved",
                       "executed replicate hierarchy remains unresolved"):
        assert unresolved not in combined


def test_picrust_seed_and_multiplicity_provenance_is_explicit():
    record = locked_data_bytes("metadata/functional/PICRUST2.md").decode("utf-8")
    assert "| Trip 5 |" in record and "330,830" in record
    trips_1_4 = next(line for line in record.splitlines() if line.startswith("| Trips 1–4"))
    assert "not recorded" in trips_1_4
    decision = json.loads((ROOT / "analysis/v3/picrust2_ecology/analysis_decision.json").read_text())
    assert decision["multiplicity"]["pathway_geography"].startswith("Benjamini-Hochberg across 200")
    assert decision["multiplicity"]["pathway_position"].startswith("Benjamini-Hochberg across 600")
    prediction = mt.section("supplement", "Prediction input and quality", "Comparison with metagenomes")
    assert "Trips 1–4 and Trip 5" in prediction and has_number(prediction, "1,227")
    # 330,830 is the Trip 5 input only; the Trips 1-4 ASV count was not recorded.
    if "330,830" in prediction:
        assert re.search(r"Trip 5[^.;]{0,60}330,830|330,830[^.;]{0,40}Trip 5", prediction), (
            "330,830 must be attributed to the Trip 5 PICRUSt2 input only"
        )
    combined = mt.combined()
    assert "Benjamini–Hochberg correction within each family" in combined
    assert "each set of 200 genus correlations" in combined
    # Seeds are recorded with each analysis; the manuscript points there.
    for path, seed in (
        ("analysis/v3/spatial_turnover_rescue/results/claim_verdict.json", 20260723),
        ("analysis/v3/xrf_community_clr/claim_verdict.json", 20260725),
        ("analysis/v3/distance_decay_turnover/claim_verdict.json", 20260728),
    ):
        assert json.loads((ROOT / path).read_text())["input"]["seed"] == seed
    assert "Random seeds are recorded in the repository" in mt.text("supplement")


def test_confirmed_laboratory_volumes_and_rinsing_are_reported():
    combined = mt.combined()
    for value in ("4,000rpm", "45-second cycle", "8\\muL", "12\\muL",
                  "31.75\\muL", "28.75\\muL", "24.75\\muL", "0.8:1", "target molarity",
                  "room temperature", "Milli-Q", "Kimwipes"):
        assert value in combined, value
    assert "programme-to-library allocation before submission" not in combined
    assert "replicate 2 was requested" not in combined.lower()


def test_latest_confirmed_primary_pcr_and_ph_protocol():
    combined = mt.combined()
    pcr = mt.section("supplement", "Amplification and indexing", "Library quantification")
    for fact in ("35 cycles", "95^\\circC for 10s", "65^\\circC for 30s", "75^\\circC for 5min",
                 "10× Mg^2+-containing buffer", "2.5mM each", "5U/\\muL"):
        assert fact in pcr, fact
    assert "V_i = V_poolC_pool/(N C_i)" in mt.text("supplement")
    for fact in ("13-620-112", "hula mixer"):
        assert fact in mt.text("main"), fact
    assert "green stability indicator" in mt.text("supplement")
    assert "30 cycles" not in combined
    assert "14 July, used replicate 1" not in combined
    assert "illumination unspecified" not in combined
    # The locked shared pH release (v1.0.1) records which field replicate was
    # assayed for Trip 4; the manuscript must report that specimen identity.
    confirmation = json.loads(
        locked_data_bytes("metadata/samples/ph/versions/EQ-PH-SHARED-v1.0.1/confirmation.json")
    )
    assert "replicate 2 was used for every Trip 4 pH assay" in confirmation["fact"]
    assert re.search(r"Trip 4[^.]{0,80}replicate 2|replicate 2[^.]{0,80}Trip 4", combined), (
        "Trip 4 pH assays used field replicate 2 (EQ-PH-SHARED-v1.0.1)"
    )


def test_final_ecology_cohort_is_reported_without_repair_history():
    main = mt.text("main")
    reads = mt.section("main", "Read processing", "Relic-DNA removal experiment")
    assert "10 samples that had <1,000 reads, leaving 1,237 samples" in reads
    # The abstract and Results report the 1,227 samples from sites 1-60.
    assert r"obtained $1{,}227$ $16$S rRNA profiles" in mt.source("main")
    assert "Four sites sampled only during T1 contributed 10 additional samples" in main
    assert has_number(main, "1,227")
    combined = mt.combined()
    for process_history in (
        "multiset difference",
        "older feature-table.tsv",
        "Source-column alignment was repaired",
        "not malformed records",
        "erroneous duplicate",
        "The replacement summed sequencing replicates",
    ):
        assert process_history not in combined


def test_reproducibility_boundary_does_not_promise_a_missing_status_table(
    main_tex,
):
    flat = " ".join(main_tex.split())
    assert "Supplementary Information records the current component" not in flat
    assert "The public repository provides the staged inputs and downstream analyses" in flat
    assert "records required for a complete raw-read replay" in flat
    assert "10.5281/zenodo.23134168" in flat


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


def test_xrf_geographic_scope_is_bounded():
    verdict = json.loads((ROOT / "analysis/v3/xrf_community_clr/claim_verdict.json").read_text())
    primary = verdict["primary"]
    assert primary["n_observations"] == 621 and primary["n_sites"] == 60
    assert verdict["all_partial_r2_below_0_01"] is True
    main = mt.section("main", "Soil chemical composition is a weak correlate", r"\begintable")
    # Within-site association, reported at two decimals in the main text and
    # three in the supplement.
    assert "Within sites" in main
    assert f"partial RDA, {pct(primary['partial_r2'], 2)}%, p={fmt(primary['permutation_p'], 3)}" in main
    supplement = mt.section("supplement", "Laboratory XRF elemental axis", "Landform")
    assert f"{pct(primary['partial_r2'], 3)}% of the remaining compositional variation (p={fmt(primary['permutation_p'], 3)})" in supplement
    assert f"Moran's I={fmt(primary['residual_moran_i'], 3)}, p={fmt(primary['residual_moran_p'], 3)}" in supplement
    assert has_number(supplement, "621")
    # No interval is reported because resampled partial R2 is biased upward.
    assert "biased upward" in supplement
    # Total elements are not salinity or soluble ions.
    discussion = mt.text("main")
    assert "XRF measures total elemental composition rather than soluble ions or salinity" in discussion
    assert "salinity gradient" not in discussion.lower()


def test_xrf_primary_status_sensitivity_matches_the_canonical_table():
    sensitivity = pd.read_csv(
        ROOT / "analysis/xrf_audit/xrf_aggregation_sensitivity.tsv",
        sep="\t",
    )
    affected = sensitivity[
        (sensitivity["workflow"] == "lab_t5")
        & (sensitivity["groups_current_differs_from_primary"] > 0)
    ]
    assert len(affected) == 11
    assert set(affected["current_rule"]) == {"last_reported"}
    minimum = float(affected["spearman_current_vs_primary"].min())
    assert 0.994 <= minimum < 0.995
    workflows = mt.section("supplement", "Measurement workflows", "XRF element axis")
    assert "T1–T4 records retained the largest positive reading and T5 retained the most recent reading" in workflows
    assert f"consistent for the {len(affected)} affected" in workflows
    assert f"ρ ≥ {fmt(minimum - 0.0005, 3)}" in workflows or f"ρ≥{fmt(minimum - 0.0005, 3)}" in workflows
    # The audit compared the T5 most-recent rule with the workbook's
    # primary-status rows, not with a largest-value rebuild.
    assert "using the ``largest'' rule" not in workflows, (
        "T5 sensitivity compared most-recent readings with primary-status readings"
    )
    assert "primary-status sensitivity" not in mt.text("main")
    assert "maximum-positive sensitivity" not in mt.text("supplement")


def test_all_supplementary_floats_are_cited():
    labels = {
        label for label in mt.labels("supplement")
        if label.startswith(("tab:", "fig:"))
    }
    assert labels
    cited = mt.refs("supplement") | {
        ref.removeprefix("si-") for ref in mt.refs("main") if ref.startswith("si-")
    }
    uncited = sorted(labels - cited)
    assert uncited == [], f"uncited supplementary floats: {uncited}"


def test_primary_aitchison_configuration_is_explicit():
    verdict = json.loads(
        (
            ROOT / "analysis/v3/xrf_community_clr/claim_verdict.json"
        ).read_text()
    )
    assert verdict["input"]["primary_taxon_count"] == 200
    assert verdict["input"]["primary_zero_treatment"] == "pseudocount_0.5"
    spatial = json.loads((ROOT / "analysis/v3/spatial_turnover_rescue/results/claim_verdict.json").read_text())["input"]
    assert spatial["prevalence_threshold"] == 0.2 and spatial["pseudocount"] == 0.5
    assert spatial["minimum_group_reads"] == 2000
    methods = mt.section("main", "Beta diversity (between-sample)", "Transect model")
    assert "discarded profiles with < 2,000 reads, yielding 630 profiles across 60 sites" in methods
    assert "200 genera with highest mean relative abundance present in at least 20% of profiles" in methods
    assert "adding a 0.5 pseudocount" in methods
    assert "Aitchison distance, the Euclidean distance between CLR profiles" in methods


def test_phylogenetic_diagnostic_matches_the_artifact():
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
    # The diagnostic failed its prerequisite and is not reported; the current
    # manuscript must not make a phylogenetic-signal or niche-conservatism claim.
    combined = mt.combined().lower()
    for claim in ("phylogenetic signal", "niche conservatism", "niche difference",
                  "phylogenetically conserved"):
        assert claim not in combined, claim


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


def test_rainfall_scope_is_consistent_across_the_whole_manuscript(rain):
    assert rain["decision"]["analysis_status"] == "temporally_localized_association_borderline"
    main = mt.text("main")
    si = mt.text("supplement")
    results = mt.section("main", "Exploratory models suggest an early richness response", r"\begintable")
    assert "only six such permutations were possible" in results
    assert "smallest attainable value = 1/6" in results
    assert "2023, 2024 and 2025" in si
    assert "exchangeability of the three observed annual rainfall fields" in si
    assert "including all nuisance coefficients" in si
    assert "separately within each product" in si
    assert "not distinguishable from the permuted assignments" in si
    # Rainfall is framed as exploratory and as a sampling window.
    assert "Exploratory models" in main
    assert "candidate window for additional sampling" in main or "window for event-based sampling" in main
    assert "rain caused" not in (main + si).lower()
    assert "figS_campaign_rainfall.pdf" not in mt.source("main") + mt.source("supplement")
    assert "rain_calendar_refit.pdf" in mt.source("supplement")


def test_every_main_text_cross_reference_resolves_in_main_text():
    """Each \\ref resolves in its own document or, via xr, in the other."""
    main_labels = mt.labels("main")
    supplement_labels = mt.labels("supplement")
    for ref in mt.refs("main"):
        if ref.startswith("si-"):
            assert ref.removeprefix("si-") in supplement_labels, ref
        else:
            assert ref in main_labels, ref
    for ref in mt.refs("supplement"):
        if ref.startswith("main-"):
            assert ref.removeprefix("main-") in main_labels, ref
        else:
            assert ref in supplement_labels, ref
    assert r"\externaldocument[si-]{supplement}" in mt.source("main")
    assert r"\externaldocument[main-]{main}" in mt.source("supplement")
    assert "fig:network-function" not in main_labels | mt.refs("main")
    # The five main-text figures after Rund Tawfiq's figure revision (Oct 2026).
    assert {"fig:landscape", "fig:composition-beta", "fig:compartment-differences",
            "fig:env-gradients", "fig:function-controls"} <= main_labels


def test_relocated_detail_survives_in_the_supplement():
    """Analyses outside the current manuscript keep their numbers; no unsupported claims remain."""
    network = json.loads((ROOT / "analysis/v3/network_rescue/results/claim_verdict.json").read_text())
    assert round(network["selected_combined_null_edge_ratio"], 5) == 0.05374
    assert network["original_claim_status"] == "retire"
    assert network["mechanistic_interpretation_permitted"] is False
    redundancy = pd.read_csv(
        ROOT / "analysis/v3/functional_redundancy_results/functional_redundancy_null.tsv", sep="\t"
    )
    every = redundancy[(redundancy.annotation_basis == "ko_copy_count") & (redundancy.group == "All")].iloc[0]
    assert [round(every[c], 5) for c in ("taxonomic_median_bray", "observed_functional_median_bray", "null_median")] == [
        0.82272, 0.13386, 0.08733
    ]
    combined = mt.combined()
    for retired in ("Graphical Lasso", "co-occurrence network", "hub taxa", "keystone",
                    "functional redundancy", "intact genome annotation profiles among the fixed",
                    "For KO-copy profiles, the observed median"):
        assert retired.lower() not in combined.lower(), retired
    # Shotgun provenance retained in the main Methods and Software.
    main = mt.text("main")
    assert "eggNOG-mapper v2.1.12" in main and "CoverM" in main
    assert "991 metagenome-assembled genomes" in main and "990 had KO annotations" in main
    # MAG read recruitment is recomputed from the pinned CoverM profiles.
    import tarfile

    with tarfile.open(DATA_REPO / "metadata/metagenome/coverm_profiles.tar.gz") as archive:
        recruited = []
        for member in archive.getmembers():
            if not member.name.endswith(".tsv"):
                continue
            profile = pd.read_csv(archive.extractfile(member), sep="\t")
            unmapped = profile.loc[profile.iloc[:, 0] == "unmapped"].iloc[0, 1]
            recruited.append(100 - float(unmapped))
    assert len(recruited) == 150
    mean = sum(recruited) / len(recruited)
    for name in ("main", "supplement"):
        assert f"{fmt(mean, 1)}% of shotgun reads" in mt.text(name), name
    supplement = mt.text("supplement")
    assert has_number(supplement, "990") and has_number(supplement, "2,000")


def test_network_calibration_diagnostic_is_named_correctly():
    alpha = pd.read_csv(ROOT / "analysis/v3/network_rescue/results/alpha_calibration.tsv", sep="\t")
    assert len(alpha) > 0
    verdict = json.loads((ROOT / "analysis/v3/network_rescue/results/claim_verdict.json").read_text())
    assert verdict["maximum_expected_false_fraction_gate"] == 0.1
    assert verdict["selected_combined_null_edge_ratio"] <= verdict["maximum_expected_false_fraction_gate"]
    assert verdict["comparative_density_contrasts_passing"] == 0
    # The network analysis is not part of the current manuscript; if it is
    # reintroduced the calibration must be named as an expected false fraction.
    combined = mt.combined()
    assert "null-selection fraction" not in combined
    if "network" in combined.lower():
        assert "expected false fraction" in combined.lower()


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


def test_abstract_interprets_the_supported_soil_position_result():
    abstract = mt.abstract()
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
    verdict = json.loads((ROOT / "analysis/v3/compartment_composition/claim_verdict.json").read_text())
    for contrast in ("Rhizosphere-Surface", "Rhizosphere-Deep"):
        assert verdict["contrasts"][contrast]["status"] == "supported"
    estimands = _alpha_estimands()
    for metric in ("richness_rare", "normalized_shannon"):
        for contrast in ("Rhizosphere-Surface", "Rhizosphere-Deep"):
            row = estimands.loc[("common_profiles", metric, "campaign_matched", contrast)]
            assert row["mean"] < 0 and row["q_three_contrasts"] < 0.05, (metric, contrast)
    assert "root-adjacent soil had the lowest alpha diversity" in abstract
    assert "communities distinct from bulk soil" in abstract
    assert "dissimilarity increased with geographic distance in all sampled compartments" in abstract
    assert "Exploratory models" in abstract
    assert "pseudo-F" not in abstract
    assert "partial R" not in abstract


def _alpha_estimands():
    table = pd.read_csv(
        ROOT / "analysis/v3/paired_alpha_sensitivity_20260909/paired_alpha_estimands.tsv",
        sep="\t",
    )
    return table.set_index(["cohort", "metric", "estimand", "contrast"])


def _q_text(q: float) -> str:
    """The supplement prints q to three decimals, or as m×10^-k when that rounds to zero."""
    if fmt(q, 3) != "0.000":
        return fmt(q, 3)
    exponent = int(f"{q:e}".split("e")[1])
    mantissa = q / 10 ** exponent
    return f"{fmt(mantissa, 0)}×10^{exponent}"
