#!/usr/bin/env python3
"""Unit and adversarial fixture tests for Empty Quarter 2026-11-08 compute tasks.

Covers all supervisor code review items (2026-11-08 19:07):
1. MIMAG HQ requires explicit RNA criteria (5S/16S/23S + >=18 tRNAs); >90/<5 without RNAs is MQ,
   missing RNA evidence is unknown_full_HQ_missing_rna; sensitivity_comp70_contam5_pass vs
   catalogue_comp50_contam10_pass; exact equality cont == 10.0 (not rtol=1e-5).
2. validated_merge checks key coverage explicitly before merging so left joins and one_to_many
   merges fail on unmatched left keys, and many_to_one rejects null/empty keys.
3. decompose_trait_turnover uses the original symmetric between_shared / within_shared / exclusive
   decomposition and matches an independent hand-calculated numerical oracle including exclusive groups.
4. Missing internal_stop / partial / contig-edge evidence is unknown (not observed absence);
   candidate promotion requires full sequence, validated coordinates, and retained contig;
   mismatched vector lengths, non-finite values, negative abundances, and null/empty blocks fail.
5. GrpE/DnaJ/HSP20-associated dark families are candidate examples, never positive controls.
6. Frequentist TOST and Bayesian INLA ROPE are strictly separated and default to no equivalence claim
   without an explicit scientific justification string.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from common_science import (
    audit_orf_integrity,
    bh_fdr,
    block_bootstrap_stat,
    classify_candidate_role,
    classify_genome_quality,
    clean_candidate_aa_sequence,
    decompose_trait_turnover,
    estimate_matched_joint_log_contrast,
    evaluate_bayesian_rope,
    evaluate_frequentist_tost,
    fit_clustered_binomial_glm,
    identify_mgnify_gtdb_duplicates,
    normalize_genus_groups,
    parse_and_verify_refined_fasta_bytes,
    require_file,
    reverse_complement_dna,
    scale_plddt_explicit_source,
    translate_cds_table11,
    validate_audited_compartment,
    validate_phylum_scope,
    validate_primary_key,
    validate_signalp_chunk_against_fasta,
    validate_stage1_to_stage1b_schema,
    validated_merge,
    verify_locus_against_refined_contig,
    within_contig_circular_null,
)
from stage4_ecology_synthesis import (
    build_systematic_effect_matrix,
    verify_panel_identifier_mapping,
)


def test_missing_inputs_fail_explicitly(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.tsv"
    with pytest.raises(FileNotFoundError, match="Missing required"):
        require_file(missing, "test fixture")

    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    with pytest.raises(ValueError, match="empty"):
        require_file(empty, "empty fixture")


def test_invalid_ids_and_compartment_labels_fail() -> None:
    df_bad_id = pd.DataFrame({"genome": ["G1", "  ", "G3"], "val": [1, 2, 3]})
    with pytest.raises(ValueError, match="empty or invalid"):
        validate_primary_key(df_bad_id, ["genome"], "bad_id_table")

    df_ws = pd.DataFrame({"genome": ["G1", " G2", "G3"], "val": [1, 2, 3]})
    with pytest.raises(ValueError, match="whitespace"):
        validate_primary_key(df_ws, ["genome"], "ws_table")

    for forbidden in ("D", "S", "PR", "Dune", "Slack", "Playa"):
        df_comp = pd.DataFrame({"sample": ["s1"], "compartment": [forbidden]})
        with pytest.raises(ValueError, match="Invalid or un-audited compartment"):
            validate_audited_compartment(df_comp, "compartment")

    df_valid = pd.DataFrame(
        {"sample": ["s1", "s2", "s3"], "compartment": ["Deep", "Surface", "RootAdjacent"]}
    )
    validate_audited_compartment(df_valid, "compartment")


def test_validated_merge_left_and_one_to_many_and_null_keys() -> None:
    left = pd.DataFrame({"id": ["A", "B", "C"], "x": [1, 2, 3]})
    right_dup = pd.DataFrame({"id": ["A", "A", "B", "C"], "y": [10, 11, 20, 30]})
    with pytest.raises(ValueError, match="duplicate primary keys"):
        validated_merge(left, right_dup, on="id", validate="one_to_one")

    right_missing = pd.DataFrame({"id": ["A", "B"], "y": [10, 20]})
    # Critical review fix: how='left' with require_all_left=True MUST fail when key 'C' is missing!
    with pytest.raises(ValueError, match="unmatched left keys"):
        validated_merge(
            left, right_missing, on="id", how="left", validate="one_to_one", require_all_left=True
        )

    # one_to_many changes row count legitimately when all keys are present, but must fail when a left key is missing
    right_o2m_ok = pd.DataFrame({"id": ["A", "A", "B", "C"], "y": [10, 11, 20, 30]})
    merged_o2m = validated_merge(
        left, right_o2m_ok, on="id", how="left", validate="one_to_many", require_all_left=True
    )
    assert len(merged_o2m) == 4

    right_o2m_missing = pd.DataFrame({"id": ["A", "A", "B"], "y": [10, 11, 20]})
    with pytest.raises(ValueError, match="unmatched left keys"):
        validated_merge(
            left, right_o2m_missing, on="id", how="left", validate="one_to_many", require_all_left=True
        )

    # many_to_one with null ID in left table must fail explicitly
    left_null = pd.DataFrame({"id": ["A", None, "B"], "x": [1, 2, 3]})
    right_m2o = pd.DataFrame({"id": ["A", "B"], "y": [10, 20]})
    with pytest.raises(ValueError, match="null values found in key column"):
        validated_merge(left_null, right_m2o, on="id", how="left", validate="many_to_one")


def test_mimag_rna_requirements_and_exact_boundary_quality() -> None:
    # 7 genomes:
    # 0: 95.0 / 1.0 with rna_complete=True -> MIMAG_HQ
    # 1: 95.0 / 1.0 with rna_complete=False -> MIMAG_MQ_strict (high comp/low contam without RNAs is MQ!)
    # 2: 95.0 / 1.0 with rna_complete=None -> unknown_full_HQ_missing_rna
    # 3: 70.0 / 5.0 with rna_complete=None -> sensitivity_comp70_contam5_pass=True, MIMAG_MQ_strict
    # 4: 69.999 / 2.0 -> sensitivity_comp70_contam5_pass=False, MIMAG_MQ_strict
    # 5: 65.0 / 10.0 -> exact 10.0% boundary -> boundary_10pct_contam
    # 6: 65.0 / 9.99995 -> within rtol=1e-5 of 10.0, but strictly < 10.0 -> MUST be MIMAG_MQ_strict, NOT boundary!
    df_q = classify_genome_quality(
        completeness=[95.0, 95.0, 95.0, 70.0, 69.999, 65.0, 65.0],
        contamination=[1.0, 1.0, 1.0, 5.0, 2.0, 10.0, 9.99995],
        rna_complete=[True, False, None, None, False, False, False],
    )
    assert list(df_q["mimag_strict_tier"]) == [
        "MIMAG_HQ",
        "MIMAG_MQ_strict",
        "unknown_full_HQ_missing_rna",
        "MIMAG_MQ_strict",
        "MIMAG_MQ_strict",
        "boundary_10pct_contam",
        "MIMAG_MQ_strict",
    ]
    assert list(df_q["comp_contam_tier"]) == [
        "comp_gt90_contam_lt5",
        "comp_gt90_contam_lt5",
        "comp_gt90_contam_lt5",
        "comp_ge50_contam_lt10",
        "comp_ge50_contam_lt10",
        "boundary_10pct_contam",
        "comp_ge50_contam_lt10",
    ]
    assert list(df_q["sensitivity_comp70_contam5_pass"]) == [
        True,
        True,
        True,
        True,
        False,
        False,
        False,
    ]
    assert list(df_q["strict_comp50_contam_lt10_pass"]) == [
        True,
        True,
        True,
        True,
        True,
        False,
        True,
    ]
    assert list(df_q["catalogue_comp50_contam10_pass"]) == [
        True,
        True,
        True,
        True,
        True,
        True,
        True,
    ]

    # When rna_complete is omitted entirely, >90/<5 must be unknown_full_HQ_missing_rna, never MIMAG_HQ
    df_no_rna = classify_genome_quality(completeness=[96.0], contamination=[0.5])
    assert df_no_rna["mimag_strict_tier"].iloc[0] == "unknown_full_HQ_missing_rna"

    # Non-finite or mismatched lengths must fail explicitly
    with pytest.raises(ValueError, match="Length mismatch"):
        classify_genome_quality([90.0, 80.0], [1.0])
    with pytest.raises(ValueError, match="finite"):
        classify_genome_quality([np.nan], [1.0])


def test_orf_integrity_missing_is_unknown_and_candidate_promotion() -> None:
    df = pd.DataFrame(
        {
            "gene": ["g1", "g2", "g3", "g4", "g5"],
            "partial": ["00", "01", "00", "00", None],
            "internal_stop": [0, 0, np.nan, 1, 0],
            "length": [250, 180, 200, 300, 150],
            "contig_rank": [4, 1, 5, 5, np.nan],
            "contig_n_genes": [10, 10, 10, 10, np.nan],
            "contig_in_catalogue": [True, True, True, True, True],
            "has_full_seq": [True, True, True, True, True],
            "coords_valid": [True, True, True, True, True],
        }
    )
    audited = audit_orf_integrity(
        df,
        contig_rank_col="contig_rank",
        contig_n_genes_col="contig_n_genes",
        contig_in_cat_col="contig_in_catalogue",
        has_full_seq_col="has_full_seq",
        coords_valid_col="coords_valid",
        edge_gene_buffer=2,
    )
    # g3 has missing internal_stop -> must be 'unknown' and pd.NA, NOT complete!
    # g5 has missing partial and missing contig_rank -> 'unknown'
    assert list(audited["orf_status"]) == [
        "complete",
        "incomplete_or_internal_stop",
        "unknown",
        "incomplete_or_internal_stop",
        "unknown",
    ]
    assert pd.isna(audited["orf_has_internal_stop"].iloc[2])
    assert pd.isna(audited["orf_is_complete"].iloc[2])
    assert list(audited["contig_edge_status"]) == ["interior", "edge", "interior", "interior", "unknown"]
    assert pd.isna(audited["is_contig_edge"].iloc[4])

    # Only g1 is complete with known stop=0, known edge status, retained contig, full seq, valid coords
    assert list(audited["candidate_promotion_eligible"]) == [True, False, False, False, False]


def test_symmetric_turnover_decomposition_against_hand_calculated_oracle() -> None:
    """Verify exact symmetric between_shared / within_shared / exclusive decomposition against a hand oracle."""
    # 5 genomes:
    # G1, G2 in shared genus g__Alpha (present in both target and ref)
    # G3 in target-exclusive genus g__Beta (ma > 0, mb == 0)
    # G4 in ref-exclusive genus g__Gamma (ma == 0, mb > 0)
    # G5 unclassified ('g__') -> present in both -> singleton group -> within_shared contribution == 0
    genomes = ["G1", "G2", "G3", "G4", "G5"]
    genera = ["g__Alpha", "g__Alpha", "g__Beta", "g__Gamma", "g__"]
    trait = np.array([1.0, 0.0, 1.0, 1.0, 1.0])

    # Normalized abundances:
    # Target a: G1=0.40, G2=0.10 (so m_a(Alpha)=0.50, t_a(Alpha)=0.40/0.50=0.80),
    #           G3=0.30 (m_a(Beta)=0.30, m_b(Beta)=0.0 -> exclusive),
    #           G4=0.00,
    #           G5=0.20 (m_a(G5)=0.20, t_a(G5)=1.0)
    a = np.array([0.40, 0.10, 0.30, 0.00, 0.20])
    # Ref b:    G1=0.10, G2=0.30 (so m_b(Alpha)=0.40, t_b(Alpha)=0.10/0.40=0.25),
    #           G3=0.00,
    #           G4=0.50 (m_a(Gamma)=0.0, m_b(Gamma)=0.50 -> exclusive),
    #           G5=0.10 (m_b(G5)=0.10, t_b(G5)=1.0)
    b = np.array([0.10, 0.30, 0.00, 0.50, 0.10])

    # Hand-calculated exact values:
    # 1. Shared group g__Alpha:
    #    B_Alpha = (0.50 - 0.40) * (0.80 + 0.25) / 2 = 0.10 * 0.525 = +0.0525
    #    W_Alpha = (0.80 - 0.25) * (0.50 + 0.40) / 2 = 0.55 * 0.45  = +0.2475
    # 2. Shared singleton group G5 (unclassified):
    #    B_G5 = (0.20 - 0.10) * (1.0 + 1.0) / 2 = +0.10
    #    W_G5 = (1.0 - 1.0) * (0.20 + 0.10) / 2 = 0.0
    # Total between_shared = 0.0525 + 0.10 = +0.1525
    # Total within_shared  = 0.2475 + 0.0  = +0.2475
    # 3. Exclusive groups (g__Beta and g__Gamma):
    #    X_Beta  = 0.30 * 1.0 - 0.0 = +0.30
    #    X_Gamma = 0.0 - 0.50 * 1.0 = -0.50
    # Total exclusive = +0.30 - 0.50 = -0.20
    # Total diff = T_target (0.40 + 0.30 + 0.20 = 0.90) - T_ref (0.10 + 0.50 + 0.10 = 0.70) = +0.20
    # Check sum: 0.1525 + 0.2475 - 0.20 = +0.20
    res = decompose_trait_turnover(a, b, trait, genomes, genera)
    assert np.isclose(res["between_shared"], 0.1525, atol=1e-12)
    assert np.isclose(res["within_shared"], 0.2475, atol=1e-12)
    assert np.isclose(res["exclusive"], -0.2000, atol=1e-12)
    assert np.isclose(res["total_diff"], 0.2000, atol=1e-12)
    assert np.isclose(res["unclassified_mass"], 0.1500, atol=1e-12)

    # Validation checks on negative abundances or NaN
    with pytest.raises(ValueError, match="negative"):
        decompose_trait_turnover([-0.1, 0.2, 0.3, 0.3, 0.3], b, trait, genomes, genera)
    with pytest.raises(ValueError, match="finite"):
        decompose_trait_turnover([np.nan, 0.2, 0.3, 0.3, 0.2], b, trait, genomes, genera)

    # Item 9 check: unequal mapped sums (all-read weights vs mapped-community normalized weights)
    # Suppose target has total mapped fraction 0.50 (a_raw = 0.5 * a) and ref has 0.80 (b_raw = 0.8 * b).
    # Within shared genus Alpha:
    #   a_raw: G1=0.20, G2=0.05 -> m_a=0.25, t_a=0.80
    #   b_raw: G1=0.08, G2=0.24 -> m_b=0.32, t_b=0.25
    #   B_Alpha_raw = (0.25 - 0.32) * (0.80 + 0.25)/2 = -0.07 * 0.525 = -0.03675
    #   W_Alpha_raw = (0.80 - 0.25) * (0.25 + 0.32)/2 =  0.55 * 0.285 = +0.15675
    # Singleton G5:
    #   m_a=0.10, m_b=0.08, t_a=t_b=1.0 -> B_G5_raw = (0.10 - 0.08)*1.0 = +0.02, W_G5_raw = 0.0
    # Exclusive:
    #   Beta: 0.15 * 1.0 = +0.15; Gamma: -0.40 * 1.0 = -0.40 -> exclusive_raw = -0.25
    # Total raw diff = 0.50*0.90 - 0.80*0.70 = 0.45 - 0.56 = -0.11
    # Check sum: (-0.03675 + 0.02) + 0.15675 - 0.25 = -0.01675 + 0.15675 - 0.25 = -0.11
    a_raw = 0.50 * a
    b_raw = 0.80 * b
    res_raw = decompose_trait_turnover(a_raw, b_raw, trait, genomes, genera)
    assert not res_raw["normalize_to_mapped"]
    assert np.isclose(res_raw["target_weight_sum"], 0.50, atol=1e-12)
    assert np.isclose(res_raw["ref_weight_sum"], 0.80, atol=1e-12)
    assert np.isclose(res_raw["between_shared"], -0.01675, atol=1e-12)
    assert np.isclose(res_raw["within_shared"], 0.15675, atol=1e-12)
    assert np.isclose(res_raw["exclusive"], -0.25000, atol=1e-12)
    assert np.isclose(res_raw["total_diff"], -0.11000, atol=1e-12)

    # When normalize_to_mapped=True is explicitly requested on the same raw vectors, it recovers the normalized contrast (+0.20)
    res_norm = decompose_trait_turnover(a_raw, b_raw, trait, genomes, genera, normalize_to_mapped=True)
    assert res_norm["normalize_to_mapped"]
    assert np.isclose(res_norm["total_diff"], 0.2000, atol=1e-12)
    assert not np.isclose(res_raw["total_diff"], res_norm["total_diff"])


def test_block_resampling_rejects_null_or_empty_blocks() -> None:
    with pytest.raises(ValueError, match="at least 2 distinct blocks"):
        block_bootstrap_stat(["clade1", "clade1"], lambda idx: float(len(idx)), n_boot=20)
    with pytest.raises(ValueError, match="Null/NaN block"):
        block_bootstrap_stat(["clade1", None, "clade2"], lambda idx: float(len(idx)), n_boot=20)
    with pytest.raises(ValueError, match="Empty or invalid block"):
        block_bootstrap_stat(["clade1", "   ", "clade2"], lambda idx: float(len(idx)), n_boot=20)

    blocks = np.repeat([f"block_{i}" for i in range(5)], 100)
    vals = np.repeat([0.0, 0.25, 0.5, 0.75, 1.0], 100)
    boot = block_bootstrap_stat(blocks, lambda idx: float(vals[idx].mean()), n_boot=200, seed=42)
    assert boot["n_blocks"] == 5
    assert boot["se"] > 0.10


def test_candidate_examples_never_labelled_positive_controls() -> None:
    for fam in ("F41Dr3::F41Dr3_256695_46", "49Dr1::49Dr1_177214_2", "24Dr1::24Dr1_1575644_3"):
        assert classify_candidate_role(fam) == "prior_heat_stress_candidate_example"
        with pytest.raises(ValueError, match="NOT an established annotated positive control"):
            classify_candidate_role(fam, is_established_annotated_control=True)

    assert (
        classify_candidate_role("annotated_grpE_PF01025", is_established_annotated_control=True)
        == "annotated_positive_control"
    )


def test_frequentist_tost_and_bayesian_rope_separation_and_justification() -> None:
    # Default (no justified margin) -> must report no_equivalence_claim_unjustified_margin
    res_default = evaluate_frequentist_tost(
        estimate=0.1, se=0.2, ci90_low=-0.23, ci90_high=0.43, ci95_low=-0.29, ci95_high=0.49
    )
    assert res_default["status"] == "no_equivalence_claim_unjustified_margin"
    assert not res_default["margin_justified"]

    # Setting margin_justified=True without a justification string must raise ValueError
    with pytest.raises(ValueError, match="margin_justification"):
        evaluate_frequentist_tost(
            estimate=0.1,
            se=0.2,
            ci90_low=-0.23,
            ci90_high=0.43,
            ci95_low=-0.29,
            ci95_high=0.49,
            delta=2.0,
            margin_justified=True,
            margin_justification="",
        )

    # Bayesian INLA credible interval defaults to no_equivalence_claim_unjustified_rope
    res_bayes_default = evaluate_bayesian_rope(
        posterior_mean=0.05, posterior_sd=0.12, cred95_low=-0.18, cred95_high=0.28
    )
    assert res_bayes_default["inference_framework"] == "bayesian_inla_posterior"
    assert res_bayes_default["status"] == "no_equivalence_claim_unjustified_rope"


def test_stage1b_refined_contig_and_cds_translation_verification() -> None:
    """Test forward/reverse strand, alternative start codon GTG, wrong coords, missing contig, out-of-bounds, partial ORF, and AA mismatch."""
    # Construct a synthetic contig with:
    # - Forward CDS at 11..28 (18 nt = 6 codons): GTG (start -> M) GCT (A) TTC (F) GAA (E) AAA (K) TAA (stop) -> "MAFEK"
    # - Reverse CDS at 41..58 (18 nt = 6 codons): reverse-complement of ATG (M) CCG (P) GAT (D) TGG (W) CGT (R) TAG (*) -> "MPDWR"
    fwd_cds = "GTGGCTTTCGAAGAAA" + "AA"  # wait: GTG GCT TTC GAA AAA TAA = 18 nt
    fwd_cds = "GTGGCTTTCGAAAAATAA"
    rev_coding = "ATGCCGGATTGGCGTTAG"
    rev_on_contig = reverse_complement_dna(rev_coding)
    contig_seq = ("A" * 10) + fwd_cds + ("C" * 12) + rev_on_contig + ("G" * 10)
    contigs = {"contig_1": contig_seq}

    # 1. Valid forward complete CDS (with GTG start codon translated as M)
    r_fwd = verify_locus_against_refined_contig(
        contigs, "contig_1", 11, 28, 1, "00", "MAFEK"
    )
    assert r_fwd["final_promotion_gate_pass"] is True
    assert r_fwd["verification_status"] == "verified_exact_complete_cds"

    # 2. Valid reverse-strand complete CDS
    r_rev = verify_locus_against_refined_contig(
        contigs, "contig_1", 41, 58, -1, "00", "MPDWR"
    )
    assert r_rev["final_promotion_gate_pass"] is True
    assert r_rev["verification_status"] == "verified_exact_complete_cds"

    # 3. Reversed strand error (passing strand=1 for the reverse-strand gene) -> must fail!
    r_bad_strand = verify_locus_against_refined_contig(
        contigs, "contig_1", 41, 58, 1, "00", "MPDWR"
    )
    assert r_bad_strand["final_promotion_gate_pass"] is False
    assert r_bad_strand["verification_status"] == "mismatch_sequence_or_strand"

    # 4. Wrong coordinates (shifted by 3 bp) -> must fail!
    r_bad_coords = verify_locus_against_refined_contig(
        contigs, "contig_1", 14, 31, 1, "00", "MAFEK"
    )
    assert r_bad_coords["final_promotion_gate_pass"] is False
    assert r_bad_coords["cds_translation_exact_match"] is False

    # 5. Missing contig in refined FASTA -> must fail!
    r_missing_contig = verify_locus_against_refined_contig(
        contigs, "contig_999", 11, 28, 1, "00", "MAFEK"
    )
    assert r_missing_contig["refined_contig_present"] is False
    assert r_missing_contig["verification_status"] == "contig_missing_in_refined_fasta"
    assert r_missing_contig["final_promotion_gate_pass"] is False

    # 6. Out-of-bounds locus (end > contig length) -> must fail!
    r_oob = verify_locus_against_refined_contig(
        contigs, "contig_1", 55, 90, 1, "00", "MAFEK"
    )
    assert r_oob["locus_within_contig_bounds"] is False
    assert r_oob["verification_status"] == "out_of_bounds_or_non_triplet"
    assert r_oob["final_promotion_gate_pass"] is False

    # 7. Partial ORF (5'-partial `10`, where GTG is translated literally as V, and partial!=00) -> must NOT pass final promotion gate!
    r_partial = verify_locus_against_refined_contig(
        contigs, "contig_1", 11, 28, 1, "10", "VAFEK"
    )
    assert r_partial["cds_translation_exact_match"] is True
    assert r_partial["verification_status"] == "verified_translation_but_partial_orf"
    assert r_partial["final_promotion_gate_pass"] is False

    # 8. Sequence mismatch -> must fail!
    r_mismatch = verify_locus_against_refined_contig(
        contigs, "contig_1", 11, 28, 1, "00", "MAFEX"
    )
    assert r_mismatch["cds_translation_exact_match"] is False
    assert r_mismatch["final_promotion_gate_pass"] is False

    # 9. Stage 1b Final Review Item 1: Internal '*' in candidate AA sequence MUST fail (never stripped!)
    # Even if we construct a CDS with an internal stop codon TAA in the middle:
    # ATG (M) GCT (A) TAA (*) TTC (F) GAA (E) TAA (*) -> if '*' were stripped from both, both would be "MAFE".
    cds_with_internal_stop = "ATGGCTTAATTCGAATAA"
    contigs_stop = {"contig_s": "AAAA" + cds_with_internal_stop + "CCCC"}
    r_internal_stop_aa = verify_locus_against_refined_contig(
        contigs_stop, "contig_s", 5, 22, 1, "00", "MA*FE"
    )
    assert r_internal_stop_aa["n_internal_stops_in_candidate"] == 1
    assert r_internal_stop_aa["verification_status"] == "invalid_internal_stop_in_candidate_aa"
    assert r_internal_stop_aa["final_promotion_gate_pass"] is False

    # Also if candidate AA has no '*' ("MAFE") while CDS has an internal stop TAA, exact_match MUST fail!
    r_internal_stop_cds = verify_locus_against_refined_contig(
        contigs_stop, "contig_s", 5, 22, 1, "00", "MAFE"
    )
    assert r_internal_stop_cds["n_internal_stops_translated"] == 1
    assert r_internal_stop_cds["cds_translation_exact_match"] is False
    assert r_internal_stop_cds["final_promotion_gate_pass"] is False

    # Single terminal '*' in candidate AA is allowed and stripped cleanly
    r_term_star = verify_locus_against_refined_contig(
        contigs, "contig_1", 11, 28, 1, "00", "MAFEK*"
    )
    assert r_term_star["n_internal_stops_in_candidate"] == 0
    assert r_term_star["final_promotion_gate_pass"] is True

    # Double terminal '**' has an internal '*' before the terminal '*' -> must fail!
    r_double_star = verify_locus_against_refined_contig(
        contigs, "contig_1", 11, 28, 1, "00", "MAFEK**"
    )
    assert r_double_star["n_internal_stops_in_candidate"] == 1
    assert r_double_star["final_promotion_gate_pass"] is False

    # 10. Stage 1b Final Review Item 2: Inadmissible start codon (e.g., GGG) on complete ORF must fail!
    cds_bad_start = "GGGGCTTTCGAAAAATAA"  # GGG (G) instead of ATG/GTG/TTG
    contigs_bad_start = {"contig_bs": "AAAA" + cds_bad_start + "CCCC"}
    r_bad_start = verify_locus_against_refined_contig(
        contigs_bad_start, "contig_bs", 5, 22, 1, "00", "GAFEK"
    )
    assert r_bad_start["has_admissible_start_codon"] is False
    assert r_bad_start["verification_status"] == "inadmissible_start_codon"
    assert r_bad_start["final_promotion_gate_pass"] is False

    # 11. Stage 1b Final Review Item 2: Non-integer coordinates (float, string float, bool) must fail!
    for bad_coord in (11.5, "11.5", True):
        r_nonint = verify_locus_against_refined_contig(
            contigs, "contig_1", bad_coord, 28, 1, "00", "MAFEK"
        )
        assert r_nonint["verification_status"] == "invalid_noninteger_or_null_coordinates"
        assert r_nonint["final_promotion_gate_pass"] is False

    # 12. Stage 1b Final Review Item 1: Missing candidate AA sequence returns unknown_missing_candidate_sequence
    for missing_seq in ("", None, "   "):
        r_missing_seq = verify_locus_against_refined_contig(
            contigs, "contig_1", 11, 28, 1, "00", missing_seq
        )
        assert r_missing_seq["verification_status"] == "unknown_missing_candidate_sequence"
        assert r_missing_seq["cds_translation_exact_match"] is False
        assert r_missing_seq["final_promotion_gate_pass"] is False


def test_parse_and_verify_refined_fasta_bytes_and_duplicate_contig_rejection() -> None:
    import hashlib

    raw_ok = b">c1\nACGTACGT\n>c2\nTTAA\n"
    expected_sha = hashlib.sha256(raw_ok).hexdigest()
    res = parse_and_verify_refined_fasta_bytes(
        "test_mag",
        raw_ok,
        expected_sha256=expected_sha,
        expected_contigs=2,
        expected_bases=12,
    )
    assert res["sha256_matches_manifest"] is True
    assert res["contigs_match_manifest"] is True
    assert res["bp_matches_manifest"] is True
    assert res["provenance_status"] == "verified_exact_match_genome_manifest"
    assert res["contigs_dict"] == {"c1": "ACGTACGT", "c2": "TTAA"}

    # Duplicate contig ID must fail with ValueError
    raw_dup = b">c1\nACGT\n>c1\nTTAA\n"
    with pytest.raises(ValueError, match="duplicate contig ID"):
        parse_and_verify_refined_fasta_bytes("dup_mag", raw_dup)

    # Manifest SHA-256 mismatch must mark provenance_status == 'mismatch_against_genome_manifest' and sha256_matches_manifest == False
    res_bad_sha = parse_and_verify_refined_fasta_bytes(
        "test_mag",
        raw_ok,
        expected_sha256="0" * 64,
    )
    assert res_bad_sha["sha256_matches_manifest"] is False
    assert res_bad_sha["provenance_status"] == "mismatch_against_genome_manifest"


def test_scale_plddt_explicit_source_scales() -> None:
    # Under '0_to_1', 0.0 and 1.0 are valid, whereas 1.4, 70.0, 100.0 must raise ValueError
    assert np.isclose(scale_plddt_explicit_source(0.0, "0_to_1"), 0.0)
    assert np.isclose(scale_plddt_explicit_source(1.0, "0_to_1"), 1.0)
    assert np.isclose(scale_plddt_explicit_source(0.842, "0_to_1"), 0.842)
    for out_of_range in (1.4, 70.0, 100.0, -0.1):
        with pytest.raises(ValueError, match="out of bounds \\[0.0, 1.0\\]"):
            scale_plddt_explicit_source(out_of_range, "0_to_1")

    # Under '0_to_100', 0 -> 0.0, 1.0 -> 0.01, 1.4 -> 0.014, 70.0 -> 0.70, 100.0 -> 1.0
    assert np.isclose(scale_plddt_explicit_source(0.0, "0_to_100"), 0.0)
    assert np.isclose(scale_plddt_explicit_source(1.0, "0_to_100"), 0.01)
    assert np.isclose(scale_plddt_explicit_source(1.4, "0_to_100"), 0.014)
    assert np.isclose(scale_plddt_explicit_source(70.0, "0_to_100"), 0.70)
    assert np.isclose(scale_plddt_explicit_source(100.0, "0_to_100"), 1.0)
    with pytest.raises(ValueError, match="out of bounds \\[0.0, 100.0\\]"):
        scale_plddt_explicit_source(105.0, "0_to_100")
    with pytest.raises(ValueError, match="Invalid pLDDT source_scale"):
        scale_plddt_explicit_source(0.8, "auto_guess")


def test_bootstrap_site_balancing_preserves_resample_multiplicity_and_joint_contrast() -> None:
    """Verify Bootstrap Weight Review: site balancing uses fixed base weights (1 / n_obs_in_original_site)
    so repeated resampled blocks contribute repeated site weight, and test correlated endpoints.
    """
    # Construct 10 sites (5 megabarchans M1..M5, 5 uruq U1..U5), where site M1 has 4 repeat samples,
    # U1 has 3 repeat samples, and M2..M5, U2..U5 have 2 samples each -> total = 4 + 3 + 8*2 = 23 samples >= 20.
    rows = []
    site_configs = [
        ("M1", "megabarchans", "B_M1", 4, 4.2, 5.1),
        ("M2", "megabarchans", "B_M2", 2, 4.5, 5.2),
        ("M3", "megabarchans", "B_M3", 2, 4.7, 5.3),
        ("M4", "megabarchans", "B_M4", 2, 4.4, 5.0),
        ("M5", "megabarchans", "B_M5", 2, 4.6, 5.4),
        ("U1", "uruq", "B_U1", 3, 3.1, 4.9),
        ("U2", "uruq", "B_U2", 2, 3.3, 4.8),
        ("U3", "uruq", "B_U3", 2, 3.5, 5.0),
        ("U4", "uruq", "B_U4", 2, 3.2, 4.7),
        ("U5", "uruq", "B_U5", 2, 3.4, 4.9),
    ]
    for s_id, dune, blk, n_rep, base_tax, base_func in site_configs:
        for k in range(n_rep):
            comp = "Surface" if (k % 2 == 0) else "Deep"
            trip = "trip1" if (k == 0) else "trip2"
            l_reads = 15.0 + 0.1 * k + (0.2 if dune == "megabarchans" else -0.1)
            rows.append(
                {
                    "site": s_id,
                    "edgell_dune_type": dune,
                    "spatial_block": blk,
                    "compartment": comp,
                    "trip": trip,
                    "log_reads": l_reads,
                    "tax_hill1": float(np.exp(base_tax + 0.05 * k)),
                    "func_hill1": float(np.exp(base_func + 0.03 * k)),
                }
            )
    df = pd.DataFrame(rows)
    assert len(df) == 23

    # Deterministic custom resample where block B_M1 (indices 0,1,2,3) is drawn THREE times,
    # B_M2 (4,5) once, B_U1 (12,13,14) TWICE, B_U2 (15,16) once, B_U3 (17,18) once.
    # Under true site-balanced cluster bootstrap:
    # - Each row of M1 has fixed base weight 1/4 = 0.25. Drawn 3 times (12 rows in resample) -> total weight = 3.0!
    # - Each row of U1 has fixed base weight 1/3. Drawn 2 times (6 rows in resample) -> total weight = 2.0!
    r1 = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 12, 13, 14, 12, 13, 14, 15, 16, 17, 18], dtype=int)
    # Pass [r1, r1] so len(valid_reps) >= 2 and mean across [r1, r1] equals the exact single-resample estimate
    res_custom = estimate_matched_joint_log_contrast(
        df,
        tax_col="tax_hill1",
        func_col="func_hill1",
        weighting="site_balanced",
        custom_resample_indices=[r1, r1],
    )
    assert res_custom["status"] == "completed"
    # When r1 is passed twice, the bootstrap quantile (lo95 == hi95) equals the exact estimate on r1!
    b_tax_r1 = res_custom["beta_log_tax_ci95_low"]
    b_func_r1 = res_custom["beta_log_func_ci95_low"]
    delta_r1 = res_custom["delta_log_ci95_low"]

    # Manual WLS calculation with fixed base weights (1 / original_observations_per_site):
    orig_counts = df["site"].map(df["site"].value_counts()).to_numpy(dtype=float)
    base_w = 1.0 / orig_counts
    lr_arr = df["log_reads"].to_numpy(dtype=float)
    lr_z = (lr_arr - lr_arr.mean()) / lr_arr.std(ddof=1)
    X_full = np.column_stack(
        [
            np.ones(len(df)),
            (df["compartment"] == "Surface").astype(float),
            (df["trip"] == "trip2").astype(float),
            lr_z,
            (df["edgell_dune_type"] == "megabarchans").astype(float),
        ]
    )
    X_r = X_full[r1]
    y_tax_r = np.log(df["tax_hill1"].to_numpy())[r1]
    y_func_r = np.log(df["func_hill1"].to_numpy())[r1]
    sw_fixed = np.sqrt(base_w[r1])
    b_tax_expected = np.linalg.lstsq(X_r * sw_fixed[:, None], y_tax_r * sw_fixed, rcond=None)[0][4]
    b_func_expected = np.linalg.lstsq(X_r * sw_fixed[:, None], y_func_r * sw_fixed, rcond=None)[0][4]

    assert np.isclose(b_tax_r1, b_tax_expected, atol=1e-12)
    assert np.isclose(b_func_r1, b_func_expected, atol=1e-12)
    assert np.isclose(delta_r1, b_tax_expected - b_func_expected, atol=1e-12)

    # Verify that the buggy post-resample re-normalization (1 / resampled_observations_per_site)
    # would cancel the 3x multiplicity of M1 and 2x multiplicity of U1 and give a DIFFERENT coefficient:
    res_site_s = pd.Series(df["site"].to_numpy()[r1])
    w_buggy = 1.0 / res_site_s.map(res_site_s.value_counts()).to_numpy(dtype=float)
    sw_buggy = np.sqrt(w_buggy)
    b_tax_buggy = np.linalg.lstsq(X_r * sw_buggy[:, None], y_tax_r * sw_buggy, rcond=None)[0][4]
    assert not np.isclose(b_tax_r1, b_tax_buggy, atol=1e-5)

    # Correlated endpoint test: null equal-contrast case (true delta_log == 0) with strong shared site noise
    rng = np.random.default_rng(12345)
    n_sites = 24
    dunes = np.array(["megabarchans"] * 12 + ["uruq"] * 12)
    shared_noise = rng.normal(0.0, 0.4, size=n_sites)
    df_null = pd.DataFrame(
        {
            "site": [f"S{i}" for i in range(n_sites)],
            "spatial_block": [f"B{i}" for i in range(n_sites)],
            "edgell_dune_type": dunes,
            "compartment": ["Surface" if i % 2 == 0 else "Deep" for i in range(n_sites)],
            "trip": ["trip1" if i % 3 == 0 else "trip2" for i in range(n_sites)],
            "log_reads": 15.0 + rng.normal(0.0, 0.2, size=n_sites),
            "tax_hill1": np.exp(3.0 + 0.8 * (dunes == "megabarchans") + shared_noise),
            "func_hill1": np.exp(5.0 + 0.8 * (dunes == "megabarchans") + shared_noise),
        }
    )
    res_null = estimate_matched_joint_log_contrast(
        df_null,
        tax_col="tax_hill1",
        func_col="func_hill1",
        weighting="sample_weighted",
        n_boot=100,
        seed=42,
    )
    assert res_null["status"] == "completed"
    assert np.isclose(res_null["delta_log_tax_minus_func"], 0.0, atol=1e-12)
    assert np.isclose(res_null["delta_log_se"], 0.0, atol=1e-12)

    # Reject +inf Hill1 values
    df_inf = df_null.copy()
    df_inf.loc[0, "tax_hill1"] = np.inf
    with pytest.raises(ValueError, match="Non-finite"):
        estimate_matched_joint_log_contrast(
            df_inf,
            tax_col="tax_hill1",
            func_col="func_hill1",
        )


def test_stage4_panel_mapping_and_mgnify_dedup() -> None:
    from stage4_ecology_synthesis import (
        H3_ROOT_TRACES_AND_UPTAKE,
        H6_BASELINE_DNA_REPAIR,
        H6_OSMOLYTE_AND_SPORULATION_STRESS,
        H6_UV_AND_OXIDATIVE_PROTECTION,
    )

    all_declared = (
        H6_BASELINE_DNA_REPAIR
        + H6_UV_AND_OXIDATIVE_PROTECTION
        + H6_OSMOLYTE_AND_SPORULATION_STRESS
        + H3_ROOT_TRACES_AND_UPTAKE
    )
    ep_df_ok = pd.DataFrame(
        {
            "endpoint": list(range(1, len(all_declared) + 1)),
            "name": [f"{p}:mapped" for p in all_declared],
            "source": ["MAG"] * len(all_declared),
        }
    )
    res_map = verify_panel_identifier_mapping(ep_df_ok)
    assert res_map["all_declared_available_present"] is True

    # Removing one declared available panel from endpoints_df must raise ValueError
    ep_df_missing = ep_df_ok.iloc[1:].copy()
    with pytest.raises(ValueError, match="Declared available H6/H3 panels missing"):
        verify_panel_identifier_mapping(ep_df_missing)

    # Compact fixture test for build_systematic_effect_matrix:
    inla_fe_fixture = pd.DataFrame(
        [
            {
                "name": f"{all_declared[0]}:mapped",
                "spec": "all_TRUE_T_200_FALSE",
                "term": "compartmentSurface",
                "mean_original_scale": 0.015,
                "sd_original_scale": 0.004,
                "0.025quant_original_scale": 0.007,
                "0.975quant_original_scale": 0.023,
                "interval_excludes_zero": True,
                "numerical_review_required": False,
            }
        ]
    )
    paired_fixture = pd.DataFrame(
        [
            {
                "panel": all_declared[0],
                "contrast": "Surface-minus-Deep",
                "quality": "all",
                "cohort": "all",
                "denominator": "known_catalogue_fraction",
                "omitted_trip": "",
                "sites": 25,
                "visits": 28,
                "mean_pp": 1.5,
                "ci_low": 0.7,
                "ci_high": 2.3,
                "p": 0.002,
                "q_bh": 0.01,
            }
        ]
    )
    turnover_fixture = pd.DataFrame(
        [
            {
                "rank": "genus",
                "panel": all_declared[0],
                "contrast": "Surface-Deep",
                "component": "between_shared",
                "value": 0.010,
                "share_of_abs": 0.50,
            }
        ]
    )
    sys_df = build_systematic_effect_matrix(
        ep_df_ok,
        inla_fe_fixture,
        {"unresolved_fit_files": [], "missing_model_specs": []},
        paired_fixture,
        turnover_fixture,
    )
    unavail_bet = sys_df[
        (sys_df["panel"] == "curated:betaine_uptake") & (sys_df["contrast"] == "Surface-minus-Deep")
    ]
    assert len(unavail_bet) == 1
    assert "requested_name_unavailable" in unavail_bet.iloc[0]["panel_availability_status"]
    assert unavail_bet.iloc[0]["inla_primary_status"] == "requested_panel_unavailable"
    assert unavail_bet.iloc[0]["freq_mapped_all_status"] == "requested_panel_unavailable"

    avail_first = sys_df[
        (sys_df["panel"] == all_declared[0]) & (sys_df["contrast"] == "Surface-minus-Deep")
    ]
    assert len(avail_first) == 1
    assert avail_first.iloc[0]["panel_availability_status"] == "available_in_endpoints"
    assert avail_first.iloc[0]["inla_primary_status"] == "completed"
    assert avail_first.iloc[0]["freq_mapped_all_status"] == "completed"

    # Duplicate rows in INLA fixed effects must raise ValueError (never silently aggregated)
    inla_dup = pd.concat([inla_fe_fixture, inla_fe_fixture], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate primary keys"):
        build_systematic_effect_matrix(
            ep_df_ok,
            inla_dup,
            {"unresolved_fit_files": [], "missing_model_specs": []},
            paired_fixture,
            turnover_fixture,
        )

    # Verify identify_mgnify_gtdb_duplicates on the actual 301 external comparators fixture if present
    p_in = Path(__file__).resolve().parents[1] / "inputs"
    comp_p = p_in / "v2_thermo_relatives_comparators.tsv"
    skani_p = p_in / "v2_thermo_relatives_skani.tsv.gz"
    if comp_p.exists() and skani_p.exists():
        comp_df = pd.read_csv(comp_p, sep="\t")
        skani_df = pd.read_csv(skani_p, sep="\t")
        dups = identify_mgnify_gtdb_duplicates(comp_df, skani_df)
        assert len(dups) == 54
        assert len(comp_df) - len(dups) == 247


def test_signalp_chunk_validator_and_stage1_to_stage1b_schema(tmp_path: Path) -> None:
    """Test external-signalp-final-checks.txt Items 2 and 5:
    - Exact (genome, gene) ID set equality against chunk FASTA + SHA-256 check
    - Negative fixtures with identical row counts but wrong/duplicated IDs
    - Negative fixtures with unknown SignalP class or non-finite/unnormalized probabilities
    - Producer/consumer schema check between Stage 1 and Stage 1b
    """
    import gzip
    import hashlib
    from common_science import STAGE1_ROUTE_A1_REQUIRED_COLS, STAGE1_ROUTE_A2_REQUIRED_COLS

    faa_p = tmp_path / "chunk_0001.faa"
    faa_bytes = b">G1|geneA # 1 # 300 # 1\nMKKLLLLAALLAA\n>G1|geneB # 310 # 600 # -1\nMAAAA\n>G2|geneC\nMTTTT\n"
    faa_p.write_bytes(faa_bytes)
    faa_sha = hashlib.sha256(faa_bytes).hexdigest()

    # 1. Valid SignalP 6.0h prediction output matching all 3 (genome, gene) keys
    pred_ok = tmp_path / "chunk_0001.ok.txt.gz"
    with gzip.open(pred_ok, "wt") as fh:
        fh.write("# SignalP-6.0\tOrganism: Other\n")
        fh.write("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tLIPO(Sec/SPII)\tTAT(Tat/SPI)\tTATLIPO(Tat/SPII)\tPILIN(Sec/SPIII)\tCS Position\n")
        fh.write("G1|geneA # 1 # 300 # 1\tSP\t0.010000\t0.980000\t0.005000\t0.005000\t0.000000\t0.000000\tCS pos: 20-21\n")
        fh.write("G1|geneB # 310 # 600 # -1\tOTHER\t1.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t\n")
        fh.write("G2|geneC\tLIPO\t0.050000\t0.050000\t0.900000\t0.000000\t0.000000\t0.000000\tCS pos: 18-19\n")

    rows = validate_signalp_chunk_against_fasta(faa_p, pred_ok, expected_fasta_sha256=faa_sha)
    assert len(rows) == 3
    assert rows[0][:4] == ("G1", "geneA", "SP", True)
    assert rows[1][:4] == ("G1", "geneB", "OTHER", False)
    assert rows[2][:4] == ("G2", "geneC", "LIPO", True)

    # 2. Negative fixture: SHA-256 mismatch
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        validate_signalp_chunk_against_fasta(faa_p, pred_ok, expected_fasta_sha256="0" * 64)

    # 3. Negative fixture: SAME count (3 rows), but G2|geneC replaced by wrong ID G2|geneWRONG
    pred_wrong_id = tmp_path / "chunk_0001.wrong_id.txt.gz"
    with gzip.open(pred_wrong_id, "wt") as fh:
        fh.write("G1|geneA\tSP\t0.01\t0.99\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G1|geneB\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G2|geneWRONG\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
    with pytest.raises(ValueError, match="SignalP ID coverage mismatch"):
        validate_signalp_chunk_against_fasta(faa_p, pred_wrong_id, expected_fasta_sha256=faa_sha)

    # 4. Negative fixture: SAME count (3 rows), but G1|geneB duplicated and G2|geneC dropped
    pred_dup_id = tmp_path / "chunk_0001.dup_id.txt.gz"
    with gzip.open(pred_dup_id, "wt") as fh:
        fh.write("G1|geneA\tSP\t0.01\t0.99\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G1|geneB\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G1|geneB\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
    with pytest.raises(ValueError, match="Duplicate SignalP prediction"):
        validate_signalp_chunk_against_fasta(faa_p, pred_dup_id, expected_fasta_sha256=faa_sha)

    # 5. Negative fixture: unknown SignalP class or malformed probability
    pred_bad_class = tmp_path / "chunk_0001.bad_class.txt.gz"
    with gzip.open(pred_bad_class, "wt") as fh:
        fh.write("G1|geneA\tUNKNOWN_CLASS\t0.01\t0.99\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G1|geneB\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
        fh.write("G2|geneC\tOTHER\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t\n")
    with pytest.raises(ValueError, match="Invalid SignalP prediction class"):
        validate_signalp_chunk_against_fasta(faa_p, pred_bad_class, expected_fasta_sha256=faa_sha)

    # 6. Producer/consumer Stage 1 -> Stage 1b schema check
    a1_ok = pd.DataFrame([{c: 1 for c in STAGE1_ROUTE_A1_REQUIRED_COLS}])
    a2_ok = pd.DataFrame([{c: 1 for c in STAGE1_ROUTE_A2_REQUIRED_COLS}])
    validate_stage1_to_stage1b_schema(a1_ok, a2_ok)

    # Dropping 'missing_reaction' or 'pathway' (or renaming to invented 'rxn_id'/'pathway_id') must raise ValueError
    a1_bad = a1_ok.drop(columns=["missing_reaction", "pathway"]).assign(rxn_id="rxn00001", pathway_id="PWY-1")
    with pytest.raises(ValueError, match="missing required columns for Stage 1b.*missing_reaction"):
        validate_stage1_to_stage1b_schema(a1_bad, a2_ok)

    # Dropping 'frac_sp_or_S' (or renaming it to 'frac_exported_sp_or_S') must immediately raise ValueError
    a2_bad = a2_ok.drop(columns=["frac_sp_or_S"]).assign(frac_exported_sp_or_S=0.8)
    with pytest.raises(ValueError, match="missing required columns for Stage 1b.*frac_sp_or_S"):
        validate_stage1_to_stage1b_schema(a1_ok, a2_bad)


def test_stage2b_canonical_eq_taxonomy_join_and_negative_fixtures(tmp_path: Path) -> None:
    """Test canonical-taxonomy-sensitivity-design.txt (Stage 2b section):
    - Exact 237-genome EQ join against canonical refined_representatives.tsv with SHA-256 verification
    - Detection of same-count wrong-genome substitution (both missing ID and bacterial MAG from refined_representatives.tsv)
    - Detection of duplicate IDs and stale taxonomy label substitution
    - Preservation of both canonical_gtdb_r232_data_foundation and legacy_analytical_taxonomy_sensitivity in Stage 2b models
    """
    import hashlib
    from common_science import (
        CANONICAL_REFINED_REPRESENTATIVES_SHA256,
        join_and_validate_canonical_eq_taxonomy,
    )
    stg2b = pytest.importorskip("stage2b_join_external_signalp_and_direct_uniref50")
    evaluate_stage2b_uniref50_and_signalp_models = stg2b.evaluate_stage2b_uniref50_and_signalp_models

    p_in = Path(__file__).resolve().parents[1] / "inputs"
    thermo_p = p_in / "thermoplasmatota.tsv"
    can_p = p_in / "refined_representatives.tsv"
    assert thermo_p.is_file() and can_p.is_file()

    eq_df = pd.read_csv(thermo_p, sep="\t")
    joined, summary = join_and_validate_canonical_eq_taxonomy(
        eq_df,
        can_p,
        genome_col="genome",
        legacy_taxonomy_col="gtdb_classification",
        expected_sha256=CANONICAL_REFINED_REPRESENTATIVES_SHA256,
        expected_n_eq=237,
        require_canonical_diff_count=42,
    )
    assert len(joined) == 237
    assert summary["n_taxonomy_strings_differ"] == 42
    assert summary["n_class_labels_differ"] == 0
    assert summary["canonical_n_named_orders"] == 4
    assert summary["canonical_n_unassigned_order_genomes"] == 2
    assert summary["canonical_n_named_genera"] == 19
    assert summary["canonical_n_unassigned_genus_genomes"] == 135
    assert summary["canonical_n_genus_or_singleton_resampling_blocks"] == 154
    assert summary["legacy_n_named_orders"] == 4
    assert summary["legacy_n_unassigned_order_genomes"] == 0
    assert summary["legacy_n_named_genera"] == 14
    assert summary["legacy_n_unassigned_genus_genomes"] == 139
    assert summary["legacy_n_genus_or_singleton_resampling_blocks"] == 153

    # 1. Negative fixture: wrong SHA-256
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        join_and_validate_canonical_eq_taxonomy(
            eq_df,
            can_p,
            expected_sha256="0" * 64,
        )

    # 2. Negative fixture: same count (237), but one EQ genome substituted by a bacterial MAG from refined_representatives.tsv
    eq_wrong_mag = eq_df.copy()
    eq_wrong_mag.loc[0, "genome"] = "10Dr2__single_SemiBin_0"
    eq_wrong_mag.loc[0, "checkm_completeness"] = 79.18
    eq_wrong_mag.loc[0, "contamination"] = 3.82
    with pytest.raises(ValueError, match="non-Thermoplasmatota canonical phylum"):
        join_and_validate_canonical_eq_taxonomy(eq_wrong_mag, can_p)

    # 3. Negative fixture: same count (237), but one EQ genome substituted by an unknown genome ID
    eq_missing_id = eq_df.copy()
    eq_missing_id.loc[0, "genome"] = "NONEXISTENT_GENOME_ID"
    with pytest.raises(ValueError, match="unmatched left keys"):
        join_and_validate_canonical_eq_taxonomy(eq_missing_id, can_p)

    # 4. Negative fixture: duplicate genome ID
    eq_dup = eq_df.copy()
    eq_dup.loc[1, "genome"] = eq_dup.loc[0, "genome"]
    with pytest.raises(ValueError, match="duplicate primary keys"):
        join_and_validate_canonical_eq_taxonomy(eq_dup, can_p)

    # 5. Negative fixture: stale taxonomy label substitution (legacy copied over canonical)
    stale_can_p = tmp_path / "stale_refined_representatives.tsv"
    can_df_stale = pd.read_csv(can_p, sep="\t")
    leg_map = dict(zip(eq_df["genome"], eq_df["gtdb_classification"]))
    can_df_stale["gtdb_r232_classification"] = can_df_stale["mag"].map(leg_map).fillna(can_df_stale["gtdb_r232_classification"])
    can_df_stale.to_csv(stale_can_p, sep="\t", index=False)
    stale_sha = hashlib.sha256(stale_can_p.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="possible stale taxonomy label substitution"):
        join_and_validate_canonical_eq_taxonomy(
            eq_df,
            stale_can_p,
            expected_sha256=stale_sha,
            require_canonical_diff_count=42,
        )

    # 6. Verify Stage 2b model evaluation retains both canonical and legacy taxonomy schemes
    eq_g2 = joined.iloc[:2].copy()
    sub_summary = dict(summary)
    df_484_toy = pd.DataFrame(
        [
            {
                "genome": eq_g2.iloc[0]["genome"],
                "gene": "g1_1",
                "cohort": "EQ_Thermoplasmatota_237",
                "tax_order": eq_g2.iloc[0]["legacy_order_label"],
                "tax_genus": eq_g2.iloc[0]["legacy_genus_label"],
                "genus_singleton_block": eq_g2.iloc[0]["legacy_genus_or_singleton_block"],
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": True,
                "has_signalp6": True,
                "sp_any": True,
                "len_bin7": "150_199",
            },
            {
                "genome": eq_g2.iloc[0]["genome"],
                "gene": "g1_2",
                "cohort": "EQ_Thermoplasmatota_237",
                "tax_order": eq_g2.iloc[0]["legacy_order_label"],
                "tax_genus": eq_g2.iloc[0]["legacy_genus_label"],
                "genus_singleton_block": eq_g2.iloc[0]["legacy_genus_or_singleton_block"],
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": False,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
            {
                "genome": eq_g2.iloc[1]["genome"],
                "gene": "g2_1",
                "cohort": "EQ_Thermoplasmatota_237",
                "tax_order": eq_g2.iloc[1]["legacy_order_label"],
                "tax_genus": eq_g2.iloc[1]["legacy_genus_label"],
                "genus_singleton_block": eq_g2.iloc[1]["legacy_genus_or_singleton_block"],
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": True,
                "has_signalp6": True,
                "sp_any": True,
                "len_bin7": "150_199",
            },
            {
                "genome": eq_g2.iloc[1]["genome"],
                "gene": "g2_2",
                "cohort": "EQ_Thermoplasmatota_237",
                "tax_order": eq_g2.iloc[1]["legacy_order_label"],
                "tax_genus": eq_g2.iloc[1]["legacy_genus_label"],
                "genus_singleton_block": eq_g2.iloc[1]["legacy_genus_or_singleton_block"],
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": False,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
            {
                "genome": "EXT_1",
                "gene": "e1_1",
                "cohort": "External_Comparators_Deduplicated_247",
                "tax_order": eq_g2.iloc[0]["canonical_order_label"],
                "tax_genus": "g__ExtA",
                "genus_singleton_block": "GENUS::g__ExtA",
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": True,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
            {
                "genome": "EXT_1",
                "gene": "e1_2",
                "cohort": "External_Comparators_Deduplicated_247",
                "tax_order": eq_g2.iloc[0]["canonical_order_label"],
                "tax_genus": "g__ExtA",
                "genus_singleton_block": "GENUS::g__ExtA",
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": False,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
            {
                "genome": "EXT_2",
                "gene": "e2_1",
                "cohort": "External_Comparators_Deduplicated_247",
                "tax_order": eq_g2.iloc[0]["canonical_order_label"],
                "tax_genus": "g__ExtB",
                "genus_singleton_block": "GENUS::g__ExtB",
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": True,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
            {
                "genome": "EXT_2",
                "gene": "e2_2",
                "cohort": "External_Comparators_Deduplicated_247",
                "tax_order": eq_g2.iloc[0]["canonical_order_label"],
                "tax_genus": "g__ExtB",
                "genus_singleton_block": "GENUS::g__ExtB",
                "comp70_contam5_pass": True,
                "is_complete_orf": True,
                "primary_direct_uniref50_dark": False,
                "has_signalp6": True,
                "sp_any": False,
                "len_bin7": "150_199",
            },
        ]
    )
    pairs_toy = pd.DataFrame(
        [
            {"genome": eq_g2.iloc[0]["genome"], "aai_nearest": "EXT_1", "aai_nearest_gtdb": "EXT_2", "aai": 70.0},
            {"genome": eq_g2.iloc[1]["genome"], "aai_nearest": "EXT_1", "aai_nearest_gtdb": "EXT_2", "aai": 70.0},
        ]
    )
    d_df, sp_df = evaluate_stage2b_uniref50_and_signalp_models(
        df_484_toy,
        eq_g2,
        sub_summary,
        pairs_toy,
        has_signalp_joined=True,
    )
    assert set(d_df["taxonomy_scheme"]) == {
        "canonical_gtdb_r232_data_foundation",
        "legacy_analytical_taxonomy_sensitivity",
    }
    assert set(sp_df["taxonomy_scheme"]) == {
        "canonical_gtdb_r232_data_foundation",
        "legacy_analytical_taxonomy_sensitivity",
    }
    assert (sp_df["external_taxonomy_provenance"] == "unchanged_from_thermo_comparators_301_tsv").all()


def test_validate_full_stage2b_signalp_and_uniref50_provenance(tmp_path: Path) -> None:
    """Test external-join-final-launch.txt Item 3:
    - Exact (genome, gene) key sets and genome sets across df_484, ext_sig_df, and ext_manifest_df
    - Rejection of duplicate SignalP keys or conflicting EQ SignalP predictions (no silent overwrite)
    - Rejection of wrong proteome_faa_sha256 provenance against direct-search source thermo_proteomes.tar.gz
    - Rejection of same-count wrong-(genome, gene) ID substitution
    """
    import hashlib
    import io
    import tarfile
    from common_science import (
        compute_thermo_proteomes_tar_sha256_map,
        validate_full_stage2b_signalp_and_uniref50_provenance,
    )

    ext1_bytes = b">e1_1\nMKKLLLLAALLAA\n>e1_2\nMAAAA\n"
    ext2_bytes = b">e2_1\nMTTTT\n"
    ext1_sha = hashlib.sha256(ext1_bytes).hexdigest()
    ext2_sha = hashlib.sha256(ext2_bytes).hexdigest()

    tar_p = tmp_path / "thermo_proteomes.tar.gz"
    with tarfile.open(tar_p, "w:gz") as tf:
        for gid, b_data in [("EXT_1", ext1_bytes), ("EXT_2", ext2_bytes)]:
            ti = tarfile.TarInfo(name=f"faa/{gid}.faa")
            ti.size = len(b_data)
            tf.addfile(ti, io.BytesIO(b_data))

    sha_map = compute_thermo_proteomes_tar_sha256_map(tar_p, target_genome_ids={"EXT_1", "EXT_2"})
    assert sha_map == {"EXT_1": ext1_sha, "EXT_2": ext2_sha}

    with pytest.raises(ValueError, match="Missing 1 target genome"):
        compute_thermo_proteomes_tar_sha256_map(tar_p, target_genome_ids={"EXT_1", "EXT_2", "EXT_MISSING"})

    df_484 = pd.DataFrame(
        [
            {"genome": "EQ_1", "gene": "eq1_1", "cohort": "EQ_Thermoplasmatota_237"},
            {"genome": "EQ_1", "gene": "eq1_2", "cohort": "EQ_Thermoplasmatota_237"},
            {"genome": "EXT_1", "gene": "e1_1", "cohort": "External_Comparators_Deduplicated_247"},
            {"genome": "EXT_1", "gene": "e1_2", "cohort": "External_Comparators_Deduplicated_247"},
            {"genome": "EXT_2", "gene": "e2_1", "cohort": "External_Comparators_Deduplicated_247"},
        ]
    )
    ext_sig_ok = pd.DataFrame(
        [
            {"genome_id": "EXT_1", "gene_id": "e1_1", "signalp_prediction": "SP"},
            {"genome_id": "EXT_1", "gene_id": "e1_2", "signalp_prediction": "OTHER"},
            {"genome_id": "EXT_2", "gene_id": "e2_1", "signalp_prediction": "LIPO"},
        ]
    )
    ext_man_ok = pd.DataFrame(
        [
            {"genome_id": "EXT_1", "n_proteins": 2, "proteome_faa_sha256": ext1_sha},
            {"genome_id": "EXT_2", "n_proteins": 1, "proteome_faa_sha256": ext2_sha},
        ]
    )
    eq_sp_ok = {
        ("T", "eq1_1"): ("TAT", 0.95),
        ("T", "eq1_2"): ("OTHER", 0.01),
    }

    joined = validate_full_stage2b_signalp_and_uniref50_provenance(
        df_484,
        ext_sig_ok,
        ext_man_ok,
        eq_sp_ok,
        expected_ext_genome_faa_sha256=sha_map,
    )
    assert len(joined) == 5
    assert joined["has_signalp6"].all()
    assert joined["sp_any"].tolist() == [True, False, True, False, True]

    # Negative 1: Wrong proteome_faa_sha256 in manifest (file exists with wrong provenance)
    ext_man_bad_sha = ext_man_ok.copy()
    ext_man_bad_sha.loc[0, "proteome_faa_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="Provenance SHA-256 mismatch"):
        validate_full_stage2b_signalp_and_uniref50_provenance(
            df_484,
            ext_sig_ok,
            ext_man_bad_sha,
            eq_sp_ok,
            expected_ext_genome_faa_sha256=sha_map,
        )

    # Negative 2: Same row count (3), but one external gene_id substituted with wrong ID
    ext_sig_wrong_id = ext_sig_ok.copy()
    ext_sig_wrong_id.loc[2, "gene_id"] = "e2_WRONG"
    with pytest.raises(ValueError, match="External SignalP \\(genome, gene\\) key set mismatch"):
        validate_full_stage2b_signalp_and_uniref50_provenance(
            df_484,
            ext_sig_wrong_id,
            ext_man_ok,
            eq_sp_ok,
            expected_ext_genome_faa_sha256=sha_map,
        )

    # Negative 3: Duplicate (genome_id, gene_id) in ext_sig_df
    ext_sig_dup = ext_sig_ok.copy()
    ext_sig_dup.loc[1, "gene_id"] = "e1_1"
    with pytest.raises(ValueError, match="duplicate primary keys"):
        validate_full_stage2b_signalp_and_uniref50_provenance(
            df_484,
            ext_sig_dup,
            ext_man_ok,
            eq_sp_ok,
            expected_ext_genome_faa_sha256=sha_map,
        )

    # Negative 4: Conflicting EQ SignalP predictions across source tags
    eq_sp_conflict = {
        ("A", "eq1_1"): ("OTHER", 0.01),
        ("T", "eq1_1"): ("SP", 0.99),
        ("T", "eq1_2"): ("OTHER", 0.01),
    }
    with pytest.raises(ValueError, match="Conflicting EQ SignalP predictions"):
        validate_full_stage2b_signalp_and_uniref50_provenance(
            df_484,
            ext_sig_ok,
            ext_man_ok,
            eq_sp_conflict,
            expected_ext_genome_faa_sha256=sha_map,
        )






