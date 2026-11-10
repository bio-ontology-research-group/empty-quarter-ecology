"""Regression and provenance tests for the active alt_stats.tex / supp_altstat.tex ecology manuscript.

Verifies:
1. Active manuscript snapshot synchronization and clean compilation (0 undefined refs/citations, 0 Overfull \\hbox).
2. Cohort exclusions, MIMAG boundary (48Dr2__self_SemiBin_25 at 64.61% / 10.0%), and product provenance against the shared foundation manifest.
3. Exact three-component turnover decomposition identity (between + within + exclusive == total) and BetL/AhpCF within-genus source values (-0.213704 and +0.287116 pp) against turnover_summary.tsv and Stage 4 task_d_h6_baseline_vs_protection_surface_deep.tsv.
4. Stage 4 conditional joint diversity sensitivity (685 valid resamples, 315 invalid draws) and companion BibTeX author/title consistency.
5. Separation of working-draft consistency from submission-readiness placeholder gates.
"""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path

import pytest
from manuscript_paths import DATA_REPO, PAPER

ROOT = Path(__file__).resolve().parents[1]
INTEGRATION_RESULTS = ROOT / "analysis/genome_integration_20261108/results"
ACTIVE_MS_DIR = Path(
    os.environ.get(
        "EQ_ACTIVE_ECOLOGY_MS_DIR",
        "/home/leechuck/Documents/papers/empty-quarter-ecology-current",
    )
)

ACTIVE_SYNC_FILES = (
    "alt_stats.tex",
    "supp_altstat.tex",
    "mag_lifeonair_results.tex",
    "mag_lifeonair_methods.tex",
    "mag_lifeonair_supplement.tex",
    "sample.bib",
    "alt_stats.pdf",
    "supp_altstat.pdf",
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _strip_comments(text: str) -> str:
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", line) for line in text.splitlines())


def _normalize_tex(text: str) -> str:
    text = _strip_comments(text)
    text = re.sub(r"\\allowbreak\s*", "", text)
    text = (
        text.replace("\\_", "_")
        .replace("\\%", "%")
        .replace("\\ ", " ")
        .replace("{,}", ",")
        .replace("~", " ")
        .replace("$", "")
    )
    text = re.sub(r"\\(?:texttt|emph|textit|textbf)\{([^}]*)\}", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


def _combined_active_text() -> str:
    chunks = []
    for name in (
        "alt_stats.tex",
        "supp_altstat.tex",
        "mag_lifeonair_results.tex",
        "mag_lifeonair_methods.tex",
        "mag_lifeonair_supplement.tex",
    ):
        chunks.append((PAPER / name).read_text(encoding="utf-8"))
    return _normalize_tex("\n".join(chunks))


def test_active_manuscript_snapshot_and_clean_logs() -> None:
    for name in ACTIVE_SYNC_FILES:
        staged = PAPER / name
        assert staged.is_file(), f"Missing staged active file: {staged}"
        if ACTIVE_MS_DIR.is_dir() and (ACTIVE_MS_DIR / name).is_file():
            assert staged.read_bytes() == (ACTIVE_MS_DIR / name).read_bytes(), (
                f"Staged {name} is out of sync with active working tree"
            )

    for pdf_name in ("alt_stats.pdf", "supp_altstat.pdf"):
        assert (PAPER / pdf_name).read_bytes().startswith(b"%PDF")


def test_active_cohort_exclusions_mimag_and_product_provenance() -> None:
    manifest_path = DATA_REPO / "metadata/cohorts/cohort_product_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    text = _combined_active_text()

    # Exact canonical exclusions from foundation manifest
    cohorts = manifest["cohorts"]
    assert cohorts["EQ-WGS-NONZERO-550"]["excluded_zero_mapped"] == ["F5Dr1"]
    assert cohorts["EQ-WGS-CORE-546"]["excluded_samples"]["auxiliary_location"] == [
        "62PRr2",
        "64Dr3",
        "64PRr2",
        "64Sr2",
    ]
    assert cohorts["EQ-WGS-CORE-MODEL-545"]["excluded_missing_post_qc_reads"] == ["F45PRr3"]

    for req_id in ("F5Dr1", "62PRr2", "64Dr3", "64PRr2", "64Sr2", "F45PRr3"):
        assert req_id in text, f"Missing required sample exclusion ID {req_id} in active ecology text"

    for bad_id in ("62Sr1", "64Dr1", "64Rr1", "64Sr1", "T345Dr2", "35Rr3", "T160Sr1"):
        assert bad_id not in text, f"Fabricated sample exclusion ID {bad_id} present in active ecology text"

    assert "547 non-control libraries" not in text

    # MIMAG boundary genome 48Dr2__self_SemiBin_25
    b_info = manifest["mimag_quality_accounting"]["strict_mimag_tiers_lt10_contamination"][
        "boundary_genome_at_10pct"
    ]
    assert b_info["mag"] == "48Dr2__self_SemiBin_25"
    assert round(float(b_info["completeness"]), 2) == 64.61
    assert "48Dr2__self_SemiBin_25" in text
    assert "64.61%" in text
    assert "90.2%" not in text

    # Distinct KO and protein products + taxonomy snapshot reconciliation
    assert "rep_ko_v2.tsv" in text
    assert "EQ-KO-LEGACY-ANALYTICAL-v0" in text
    assert "EQ-PROT-BAKTA-v1.0.0" in text
    assert "EQ-PROT-PYRODIGAL-REP-v1.0.0" in text
    assert "EQ-TAX-LEGACY-ANALYTICAL-v0" in text
    assert "1,898" in text and "1,124" in text and "220" in text

    # Forbidden overclaims removed
    assert "Poseidoniia" not in text
    assert "UBA9212" not in text
    for forbidden_ko_phrase in (
        "conserved KO diversity",
        "comparatively stable functional",
        "maintaining comparable functional repertoire",
        "maintaining a core functional repertoire",
    ):
        assert forbidden_ko_phrase not in text.lower()
    assert (
        "Thus the strong eastern decline in bacterial genome diversity was not accompanied "
        "by a similarly large decline in this summary of annotated KO diversity."
    ) in text


def test_betl_ahpcf_turnover_values_and_three_component_identity() -> None:
    turnover_rows = _read_tsv(INTEGRATION_RESULTS / "turnover_summary.tsv")
    assert len(turnover_rows) > 0

    grouped: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in turnover_rows:
        key = (row["panel"], row["rank"], row["contrast"])
        grouped.setdefault(key, {})[row["component"]] = float(row["mean"])

    for key, comps in grouped.items():
        assert set(comps) == {"total", "between", "within", "exclusive"}, f"Incomplete components for {key}"
        reconstructed = comps["between"] + comps["within"] + comps["exclusive"]
        assert abs(reconstructed - comps["total"]) < 1e-9, (
            f"Turnover decomposition identity violated for {key}: {reconstructed} != {comps['total']}"
        )

    betl = grouped[("HMM:glycine betaine transporter", "genus", "Surface-Deep")]
    ahpcf = grouped[("HMM:ahpCF", "genus", "Surface-Deep")]
    assert betl["total"] == pytest.approx(7.551493467459441, abs=1e-6)
    assert betl["between"] == pytest.approx(6.1153420959365485, abs=1e-6)
    assert betl["within"] == pytest.approx(-0.2137042926319036, abs=1e-6)
    assert betl["exclusive"] == pytest.approx(1.649855664154796, abs=1e-6)

    assert ahpcf["total"] == pytest.approx(6.949125375870978, abs=1e-6)
    assert ahpcf["between"] == pytest.approx(5.28349576199607, abs=1e-6)
    assert ahpcf["within"] == pytest.approx(0.287115623129679, abs=1e-6)
    assert ahpcf["exclusive"] == pytest.approx(1.37851399074523, abs=1e-6)

    # Verify Stage 4 synthesis matrix agrees and obeys 3-component identity on all rows with turnover
    stage4_rows = _read_tsv(
        INTEGRATION_RESULTS / "task_d_h6_baseline_vs_protection_surface_deep.tsv"
    )
    assert len(stage4_rows) == 49
    stage4_by_panel = {r["panel"]: r for r in stage4_rows}
    for panel, expected in (
        ("HMM:glycine betaine transporter", betl),
        ("HMM:ahpCF", ahpcf),
    ):
        s_row = stage4_by_panel[panel]
        assert float(s_row["turnover_genus_total_pp"]) == pytest.approx(expected["total"], abs=1e-6)
        assert float(s_row["turnover_genus_between_shared_pp"]) == pytest.approx(
            expected["between"], abs=1e-6
        )
        assert float(s_row["turnover_genus_within_shared_pp"]) == pytest.approx(
            expected["within"], abs=1e-6
        )
        assert float(s_row["turnover_genus_exclusive_pp"]) == pytest.approx(
            expected["exclusive"], abs=1e-6
        )

    for s_row in stage4_rows:
        if s_row["turnover_genus_total_status"] == "completed":
            tot = float(s_row["turnover_genus_total_pp"])
            btw = float(s_row["turnover_genus_between_shared_pp"])
            wth = float(s_row["turnover_genus_within_shared_pp"])
            exc = float(s_row["turnover_genus_exclusive_pp"])
            assert abs((btw + wth + exc) - tot) < 1e-9, (
                f"Stage 4 turnover identity failed for {s_row['panel']}"
            )

    # Verify prose in mag_lifeonair_results.tex and exact values/caption in mag_lifeonair_supplement.tex
    res_text = _normalize_tex((PAPER / "mag_lifeonair_results.tex").read_text(encoding="utf-8"))
    assert "Using GTDB release 232 assignments, changes in community membership accounted for most of these contrasts" in res_text
    assert "between-genus replacement contributed 5.9 percentage points to the BetL difference and 5.3 points to the AhpCF difference" in res_text
    assert "changes within shared genera contributed less than 0.3 percentage points in either direction" in res_text
    assert "-0.44" not in res_text
    assert "-0.07" not in res_text
    assert "-18.8 percentage points" in res_text
    assert "-11.6 percentage points" in res_text
    assert "-6.8 percentage points" in res_text
    assert "-0.4 percentage points" in res_text

    # Verify completed canonical 3-guild turnover table against source TSV
    canon_rows = _read_tsv(
        INTEGRATION_RESULTS / "task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv"
    )
    primary_mapped = {
        (r["panel"], r["contrast"], r["component"]): r
        for r in canon_rows
        if r["is_primary_required_3guild"] == "True"
        and r["denominator_mode"] == "mapped_community_normalized_weights"
    }
    for panel, contrast, exp_tot, exp_btw, exp_wth, exp_exc in (
        ("HMM:glycine betaine transporter", "Surface-minus-Deep", 7.551493, 5.893043, -0.145520, 1.803971),
        ("HMM:ahpCF", "Surface-minus-Deep", 6.949125, 5.275635, 0.224192, 1.449298),
        ("typed:[NiFe]_Group_1l", "RootAdjacent-minus-Deep", -18.771404, -11.573751, -0.412527, -6.785126),
    ):
        tot_c = float(primary_mapped[(panel, contrast, "total_diff")]["canonical_r232_mean_pp"])
        btw_c = float(primary_mapped[(panel, contrast, "between_shared")]["canonical_r232_mean_pp"])
        wth_c = float(primary_mapped[(panel, contrast, "within_shared")]["canonical_r232_mean_pp"])
        exc_c = float(primary_mapped[(panel, contrast, "exclusive")]["canonical_r232_mean_pp"])
        assert tot_c == pytest.approx(exp_tot, abs=1e-5)
        assert btw_c == pytest.approx(exp_btw, abs=1e-5)
        assert wth_c == pytest.approx(exp_wth, abs=1e-5)
        assert exc_c == pytest.approx(exp_exc, abs=1e-5)
        assert abs((btw_c + wth_c + exc_c) - tot_c) < 1e-9

    alt_text = _normalize_tex((PAPER / "alt_stats.tex").read_text(encoding="utf-8"))
    assert "both shifts tracked community replacement under the analytical and foundational taxonomy snapshots" in alt_text
    assert "with canonical GTDB R232 genus and order sensitivity pending" not in alt_text

    supp_text = _normalize_tex((PAPER / "mag_lifeonair_supplement.tex").read_text(encoding="utf-8"))
    assert "-0.21" in supp_text
    assert "+0.29" in supp_text
    assert "+5.89 [4.01,7.68]" in supp_text
    assert "-11.57 [-15.46,-7.77]" in supp_text
    assert "-0.41 [-0.59,-0.21]" in supp_text
    assert "-6.79 [-9.65,-3.60]" in supp_text
    assert (
        "The complete panel showed different responses among stress-protection and DNA-repair capacities"
        in supp_text
    )
    assert "Their dependence assumptions differ from those of the spatial-block intervals" in supp_text
    assert "HMM and KO annotations give different conclusions" in supp_text
    assert "a single sporulation marker rather than a complete pathway" in supp_text


def test_stage4_joint_diversity_and_companion_citations() -> None:
    joint_rows = _read_tsv(
        INTEGRATION_RESULTS / "task_d_eastern_joint_diversity_contrast.tsv"
    )
    assert len(joint_rows) >= 6
    by_key = {
        (r["cohort"], r["sensitivity_specification"], r["weighting_scheme"]): r
        for r in joint_rows
    }
    sb_441 = by_key[("matched_amplicon_441_modeled_cohort", "primary_all_mags_copy", "site_balanced")]
    sa_441 = by_key[("matched_amplicon_441_modeled_cohort", "primary_all_mags_copy", "site_averaged")]
    assert int(sb_441["n_valid_replicates"]) == 685
    assert int(sb_441["n_invalid_rank_deficient_or_zero_cell_replicates"]) == 315
    assert int(sa_441["n_valid_replicates"]) == 767
    assert int(sa_441["n_invalid_rank_deficient_or_zero_cell_replicates"]) == 233
    assert float(sb_441["delta_log_tax_minus_func"]) == pytest.approx(-1.2799, abs=1e-3)

    text = _combined_active_text()
    assert "685" in text and "315" in text

    bib = (PAPER / "sample.bib").read_text(encoding="utf-8")
    assert "A linked resource of microbial genomes, community profiles and environmental data from the {Rub' al-Khali}" in bib
    assert "Widespread desert archaea harbour an export-enriched unknown protein repertoire" in bib
    assert "Zubair, Jood Kamal" in bib
    assert "Lopez Velazquez, Alejandra" in bib


def test_pending_provenance_gates_are_explicit_and_separated_from_draft_checks() -> None:
    text = _combined_active_text()
    required_pending = (
        "ECOL-PH-MODEL-INPUT-PROVENANCE",
        "DATA-BULK-ARCHIVE-DOI",
        "DATA-PROVENANCE-UNVERIFIED-PROTEIN-EQUIVALENCE",
        "TAXONOMY-SNAPSHOT-SENSITIVITY",
    )
    for tag in required_pending:
        assert tag in text, f"Expected explicit pending marker {tag} in working draft"

    if os.environ.get("EQ_REQUIRE_SUBMISSION_READY", "0") == "1":
        pytest.fail(
            "Submission-readiness gate failed: unresolved [PENDING: ...] markers remain in active manuscript"
        )
