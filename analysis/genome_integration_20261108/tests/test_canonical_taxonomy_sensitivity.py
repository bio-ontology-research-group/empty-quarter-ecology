#!/usr/bin/env python3
"""Unit and adversarial tests for canonical GTDB r232 taxonomy reconciliation, overlay, and sensitivity."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from canonical_taxonomy_sensitivity import (
    PRIMARY_REQUIRED_3GUILD_PANELS,
    build_10_contiguous_spatial_blocks,
    compute_site_blocked_genus_turnover,
    reconcile_and_validate_canonical_taxonomy,
    validate_3guild_turnover_table,
    verify_overlay_not_stale,
)


def _make_fixture_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    qual_df = pd.DataFrame(
        [
            {
                "genome": "15Dr1__self_SemiBin_98_sub",
                "gtdb_classification": "d__Archaea;p__Thermoplasmatota;c__SW-10-69-26;o__JACQPN01;f__JAIVGT01;g__JAIVGT01;s__",
                "checkm_completeness": 60.25,
                "contamination": 5.98,
                "phylum_group": "Thermoplasmatota",
            },
            {
                "genome": "59PRr1__single_SemiBin_134",
                "gtdb_classification": "d__Archaea;p__Thermoplasmatota;c__SW-10-69-26;o__SW-10-69-26;f__SW-10-69-26;g__;s__",
                "checkm_completeness": 91.43,
                "contamination": 2.35,
                "phylum_group": "Thermoplasmatota",
            },
            {
                "genome": "10Dr2__single_SemiBin_0",
                "gtdb_classification": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Micromonosporaceae;g__;s__",
                "checkm_completeness": 79.18,
                "contamination": 3.82,
                "phylum_group": "rest",
            },
        ]
    )
    thermo_df = qual_df.iloc[:2].copy().reset_index(drop=True)
    canon_df = pd.DataFrame(
        [
            {
                "mag": "15Dr1__self_SemiBin_98_sub",
                "Completeness": 60.25,
                "Contamination": 5.98,
                "gtdb_r232_classification": "d__Archaea;p__Thermoplasmatota;c__SW-10-69-26;o__JACQPN01;f__JAIVGT01;g__JAIVGT01;s__",
            },
            {
                "mag": "59PRr1__single_SemiBin_134",
                "Completeness": 91.43,
                "Contamination": 2.35,
                "gtdb_r232_classification": "d__Archaea;p__Thermoplasmatota;c__SW-10-69-26;o__SW-10-69-26;f__SW-10-69-26;g__CALIJB01;s__",
            },
            {
                "mag": "10Dr2__single_SemiBin_0",
                "Completeness": 79.18,
                "Contamination": 3.82,
                "gtdb_r232_classification": "d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Micromonosporaceae;g__Luedemannella;s__",
            },
        ]
    )
    return qual_df, thermo_df, canon_df


def test_reconcile_canonical_taxonomy_valid_fixture() -> None:
    qual_df, thermo_df, canon_df = _make_fixture_tables()
    res = reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, canon_df)
    summary = res["summary"]
    assert summary["n_total_genomes"] == 3
    assert summary["n_whole_taxonomy_strings_differ"] == 2
    assert summary["n_domain_differs"] == 0
    assert summary["n_genus_differs"] == 2
    assert summary["thermoplasmatota_237"]["n_genomes"] == 2
    assert summary["thermoplasmatota_237"]["n_whole_taxonomy_strings_differ"] == 1

    ov_qual = res["overlay_qual_df"]
    assert "legacy_gtdb_classification" in ov_qual.columns
    assert list(ov_qual["gtdb_classification"]) == list(canon_df["gtdb_r232_classification"])
    assert list(ov_qual["legacy_gtdb_classification"]) == list(qual_df["gtdb_classification"])


def test_same_count_wrong_genome_id_rejected() -> None:
    qual_df, thermo_df, canon_df = _make_fixture_tables()
    bad_canon = canon_df.copy()
    bad_canon.loc[2, "mag"] = "WRONG_GENOME_ID_SAME_ROW_COUNT"
    with pytest.raises(ValueError, match="1:1 genome membership mismatch"):
        reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, bad_canon)


def test_duplicate_genome_or_mag_ids_rejected() -> None:
    qual_df, thermo_df, canon_df = _make_fixture_tables()
    dup_qual = pd.concat([qual_df.iloc[:2], qual_df.iloc[[1]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        reconcile_and_validate_canonical_taxonomy(dup_qual, thermo_df, canon_df)

    dup_canon = pd.concat([canon_df.iloc[:2], canon_df.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, dup_canon)


def test_stale_taxonomy_label_substitution_rejected() -> None:
    qual_df, thermo_df, canon_df = _make_fixture_tables()
    res = reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, canon_df)
    stale_overlay = res["overlay_qual_df"].copy()
    stale_overlay["gtdb_classification"] = stale_overlay["legacy_gtdb_classification"]
    with pytest.raises(ValueError, match="Stale or mismatched taxonomy label substitution"):
        verify_overlay_not_stale(stale_overlay, canon_df, require_diff=True)

    missing_leg = res["overlay_qual_df"].drop(columns=["legacy_gtdb_classification"])
    with pytest.raises(ValueError, match="legacy_gtdb_classification"):
        verify_overlay_not_stale(missing_leg, canon_df, require_diff=True)


def test_domain_or_quality_mismatch_rejected() -> None:
    qual_df, thermo_df, canon_df = _make_fixture_tables()
    bad_dom = canon_df.copy()
    bad_dom.loc[2, "gtdb_r232_classification"] = "d__Archaea;p__Actinomycetota;c__Actinomycetes;o__;f__;g__;s__"
    with pytest.raises(ValueError, match="Domain labels differ"):
        reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, bad_dom)

    bad_comp = canon_df.copy()
    bad_comp.loc[0, "Completeness"] = 61.25
    with pytest.raises(ValueError, match="Completeness/contamination mismatch"):
        reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, bad_comp)


def test_spatial_block10_resampling_and_total_invariance_fixture() -> None:
    # 20 sites across 10 contiguous spatial blocks (2 sites per block)
    design_rows = [
        {"site": f"S{i:02d}", "transect_km": float(i * 15.0)}
        for i in range(20)
    ]
    site_to_block = build_10_contiguous_spatial_blocks(pd.DataFrame(design_rows), n_blocks=10)
    assert len(set(site_to_block.values())) == 10
    assert site_to_block["S00"] == site_to_block["S01"] == "spatial_block_00"

    genomes = ["G1", "G2", "G3", "G4"]
    t_vec = np.array([1.0, 0.0, 1.0, 0.0], dtype="float64")
    genera_leg = ["g__A", "g__A", "g__", "g__B"]
    genera_can = ["g__A", "g__C", "g__D", "g__B"]

    idx_tuples = []
    mat_rows = []
    rng = np.random.default_rng(42)
    for i in range(20):
        s = f"S{i:02d}"
        w_surf = rng.uniform(0.05, 0.35, size=4)
        w_deep = rng.uniform(0.05, 0.35, size=4)
        idx_tuples.append((s, "trip1", "Surface"))
        mat_rows.append(w_surf)
        idx_tuples.append((s, "trip1", "Deep"))
        mat_rows.append(w_deep)

    V_use = pd.DataFrame(
        np.vstack(mat_rows),
        index=pd.MultiIndex.from_tuples(idx_tuples, names=["site", "trip", "compartment"]),
    )

    can_res = compute_site_blocked_genus_turnover(
        V_use,
        t_vec,
        genomes,
        genera_can,
        c_target="Surface",
        normalize_to_mapped=True,
        site_to_block=site_to_block,
        n_boot=200,
        seed=20261108,
    )
    leg_res = compute_site_blocked_genus_turnover(
        V_use,
        t_vec,
        genomes,
        genera_leg,
        c_target="Surface",
        normalize_to_mapped=True,
        site_to_block=site_to_block,
        n_boot=200,
        seed=20261108,
    )

    assert can_res["total_diff"]["n_spatial_blocks"] == 10
    assert can_res["total_diff"]["n_paired_sites"] == 20
    assert can_res["_paired_sites"] == leg_res["_paired_sites"]
    assert can_res["_paired_visits"] == leg_res["_paired_visits"]
    assert abs(can_res["total_diff"]["mean_pp"] - leg_res["total_diff"]["mean_pp"]) < 1e-12
    assert abs(can_res["total_diff"]["ci95_low_pp"] - leg_res["total_diff"]["ci95_low_pp"]) < 1e-12

    # Missing site in site_to_block must raise ValueError
    incomplete_map = dict(site_to_block)
    incomplete_map.pop("S05")
    with pytest.raises(ValueError, match="missing from site_to_block"):
        compute_site_blocked_genus_turnover(
            V_use,
            t_vec,
            genomes,
            genera_can,
            c_target="Surface",
            normalize_to_mapped=True,
            site_to_block=incomplete_map,
            n_boot=50,
        )


def test_validate_3guild_turnover_table_gates() -> None:
    rows = []
    for p in PRIMARY_REQUIRED_3GUILD_PANELS:
        for ct in ("Surface-minus-Deep", "RootAdjacent-minus-Deep"):
            for dm in ("mapped_community_normalized_weights", "all_reads_unnormalized_weights"):
                comp_vals = [
                    ("total_diff", 10.0, 10.0),
                    ("between_shared", 6.0, 6.5),
                    ("within_shared", 2.0, 1.5),
                    ("exclusive", 2.0, 2.0),
                    ("unclassified_mass", 1.2, 0.8),
                ]
                for cname, l_val, c_val in comp_vals:
                    rows.append(
                        {
                            "panel": p,
                            "is_primary_required_3guild": True,
                            "panel_status": "completed",
                            "contrast": ct,
                            "denominator_mode": dm,
                            "component": cname,
                            "paired_membership_identical": True,
                            "legacy_mean_pp": l_val,
                            "legacy_ci95_low_pp": l_val - 1.0,
                            "legacy_ci95_high_pp": l_val + 1.0,
                            "canonical_r232_mean_pp": c_val,
                            "canonical_r232_ci95_low_pp": c_val - 1.0,
                            "canonical_r232_ci95_high_pp": c_val + 1.0,
                            "canonical_minus_legacy_pp": c_val - l_val,
                        }
                    )
    valid_df = pd.DataFrame(rows)
    validate_3guild_turnover_table(valid_df)

    # Missing primary panel must fail
    missing_primary_df = valid_df[valid_df["panel"] != "HMM:ahpCF"].copy()
    with pytest.raises(ValueError, match="Required primary panel 'HMM:ahpCF'"):
        validate_3guild_turnover_table(missing_primary_df)

    # Violated total_diff invariance must fail
    bad_tot_df = valid_df.copy()
    bad_tot_df.loc[bad_tot_df["component"] == "total_diff", "canonical_minus_legacy_pp"] = 0.25
    with pytest.raises(ValueError, match="Total contrast invariance violated"):
        validate_3guild_turnover_table(bad_tot_df)

    # Violated decomposition sum must fail
    bad_sum_df = valid_df.copy()
    bad_sum_df.loc[bad_sum_df["component"] == "between_shared", "canonical_r232_mean_pp"] = 9.0
    with pytest.raises(ValueError, match="Decomposition sum mismatch"):
        validate_3guild_turnover_table(bad_sum_df)


def test_canonical_real_outputs_if_present() -> None:
    out_root = Path(
        os.environ.get("EQ_CANONICAL_SCIENCE_OUT", str(Path(__file__).resolve().parents[1]))
    )
    rec_sum_p = out_root / "results/taxonomy_reconciliation_summary.json"
    if not rec_sum_p.exists():
        pytest.skip("PENDING_REAL_CANONICAL_OUTPUT_INTEGRATION: canonical reconciliation outputs not yet generated")

    rec_sum = json.loads(rec_sum_p.read_text())
    assert rec_sum["n_total_genomes"] == 5438
    assert rec_sum["n_whole_taxonomy_strings_differ"] == 1898
    assert rec_sum["n_domain_differs"] == 0
    assert rec_sum["n_phylum_differs"] == 247
    assert rec_sum["n_class_differs"] == 34
    assert rec_sum["n_order_differs"] == 220
    assert rec_sum["n_genus_differs"] == 1124

    th = rec_sum["thermoplasmatota_237"]
    assert th["n_genomes"] == 237
    assert th["n_whole_taxonomy_strings_differ"] == 42
    assert th["canonical_named_orders_n"] == 4
    assert th["canonical_unassigned_order_genomes_n"] == 2
    assert th["canonical_named_genera_n"] == 19
    assert th["canonical_unassigned_genus_genomes_n"] == 135

    turn_p = out_root / "results/task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv"
    if turn_p.exists():
        turn_df = pd.read_csv(turn_p, sep="\t")
        validate_3guild_turnover_table(turn_df)

    rg_p = out_root / "results/canonical_root_growth_taxonomy_status.json"
    if rg_p.exists():
        rg = json.loads(rg_p.read_text())
        assert rg["canonical_3guild_genus_turnover_status"] == "COMPLETED_EXACT_R232_DECOMPOSITION_10_SPATIAL_BLOCKS"
        assert "PENDING_CANONICAL_TAXONOMY_RERUN" in rg["secondary_grodon_order_genus_regression_status"]
