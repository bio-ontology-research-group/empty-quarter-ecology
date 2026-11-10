#!/usr/bin/env python3
"""Stage 4 (Task D): Ecology Synthesis, Rund H6/H3 Effect & Uncertainty Matrix,
Explicit Pathway Flags, Frequentist/Bayesian Separation & Spatially Blocked Joint Diversity Contrasts.

Implements all requirements from:
- `coordinated-revision-20261108/compute-code-review.txt`
- `coordinated-revision-20261108/stage4-scientific-review.txt` (Items 1-5)
- `coordinated-revision-20261108/bootstrap-weight-review.txt`
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    bh_fdr,
    decompose_trait_turnover,
    estimate_matched_joint_log_contrast,
    require_file,
    validate_audited_compartment,
    validate_primary_key,
)

P_ECO_INT = Path("/scratch/project/k10018/eq-ecology-integration-20261108")
P_PROT_INT = Path("/scratch/project/k10018/eq-proteome-integration-20261106")
P_ECO_MAG = Path("/scratch/project/k10018/eq-ecology-mag-20261105")
P_ARCH = Path("/scratch/project/k10018/eq-archaeal-function-20261106")
P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_RES = P_OUT_ROOT / "results"

# Exploratory post-fit synthesis groupings (separated per stage4-scientific-review.txt Item 2:
# baseline DNA repair vs UV/oxidative protection vs osmolyte/sporulation stress)
H6_BASELINE_DNA_REPAIR = [
    "legacy:base_excision_repair",
    "legacy:mismatch_repair",
    "legacy:nucleotide_excision_repair",
    "legacy:SOS_translesion",
    "legacy:nhej",
]
H6_UV_AND_OXIDATIVE_PROTECTION = [
    "legacy:photolyase",
    "legacy:carotenoid",
    "legacy:dna_protection_oxidative",
    "HMM:ahpCF",
    "KO_projection:ahpCF",
    "HMM:sodA",
    "KO_projection:sodA",
    "HMM:lycopene biosynthesis",
    "KO_projection:lycopene biosynthesis",
    "HMM:TF_redox",
    "KO_projection:TF_redox",
]
H6_OSMOLYTE_AND_SPORULATION_STRESS = [
    "curated:trehalose_strict",
    "curated:ectoine_marker",
    "curated:glucosylglycerol_strict",
    "curated:mannosylglycerate_strict",
    "HMM:trehalose_biosynthesis",
    "KO_projection:trehalose_biosynthesis",
    "HMM:betaine_biosynthesis",
    "KO_projection:betaine_biosynthesis",
    "HMM:mannosylglycerate_biosynthesis",
    "KO_projection:mannosylglycerate_biosynthesis",
    "HMM:proline_biosynthesis",
    "KO_projection:proline_biosynthesis",
    "HMM:glycine betaine/proline transport system",
    "KO_projection:glycine betaine/proline transport system",
    "HMM:glycine betaine transporter",
    "KO_projection:glycine betaine transporter",
    "HMM:choline/glycine/proline betaine transport",
    "KO_projection:choline/glycine/proline betaine transport",
    "HMM:proline transport",
    "KO_projection:proline transport",
    "HMM:opu_and",
    "KO_projection:opu_and",
    "HMM:mechanosensitive channel",
    "KO_projection:mechanosensitive channel",
    "HMM:spoIIGA",
    "KO_projection:spoIIGA",
]
H3_ROOT_TRACES_AND_UPTAKE = [
    "typed:[NiFe]_Group_1l",
    "typed:[NiFe]_Group_1h",
    "typed:[NiFe]_Group_2a",
    "reads:hydB_1l",
    "reads_excess:hydB_1l",
    "reads:hydB_1h",
    "reads_excess:hydB_1h",
    "reads:hydB_2a",
    "reads_excess:hydB_2a",
    "reads:coxL_formI",
    "reads_excess:coxL_formI",
    "HMM:dpp-1_full",
    "KO_projection:dpp-1_full",
    "HMM:dpp-1_key",
    "KO_projection:dpp-1_key",
    "HMM:dpp-2_full",
    "KO_projection:dpp-2_full",
    "HMM:dpp-2_key",
    "KO_projection:dpp-2_key",
    "HMM:opp-1_full",
    "KO_projection:opp-1_full",
    "HMM:opp-1_key",
    "KO_projection:opp-1_key",
    "HMM:aap_full",
    "KO_projection:aap_full",
    "HMM:aap_key",
    "KO_projection:aap_key",
    "HMM:gln_full",
    "KO_projection:gln_full",
    "HMM:glt_full",
    "KO_projection:glt_full",
    "HMM:amt",
    "KO_projection:amt",
    "HMM:urease",
    "KO_projection:urease",
    "HMM:ammonia assimilation to glutamate",
    "KO_projection:ammonia assimilation to glutamate",
    "HMM:ammonia assimilation to glutamine-GS",
    "KO_projection:ammonia assimilation to glutamine-GS",
    "HMM:nitrogen fixation_key",
    "KO_projection:nitrogen fixation_key",
    "HMM:ammonia oxidation to nitrite_full",
    "KO_projection:ammonia oxidation to nitrite_full",
]

# Requested legacy/shorthand names that do not exist verbatim in `endpoints.tsv`
# Must be emitted explicitly as unavailable rows with cross-reference to actual endpoints (Item 2)
REQUESTED_UNAVAILABLE_ALIASES: dict[str, tuple[str, str]] = {
    "curated:betaine_uptake": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_see_HMM_and_KO_glycine_betaine_transport_and_opu_and",
    ),
    "curated:betaine_synthesis": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_see_HMM_and_KO_betaine_biosynthesis",
    ),
    "curated:ectoine_biosynthesis": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_actual_endpoint_is_curated:ectoine_marker",
    ),
    "HMM:proVWX": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_actual_endpoint_is_HMM:glycine betaine/proline transport system",
    ),
    "KO_projection:proVWX": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_actual_endpoint_is_KO_projection:glycine betaine/proline transport system",
    ),
    "HMM:opuD_betL": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_actual_endpoints_are_HMM:glycine betaine transporter_and_HMM:opu_and",
    ),
    "KO_projection:opuD_betL": (
        "H6_osmolyte_and_sporulation_stress",
        "requested_name_unavailable_in_endpoints_actual_endpoints_are_KO_projection:glycine betaine transporter_and_KO_projection:opu_and",
    ),
    "HMM:dpp": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_HMM:dpp-1_full_dpp-1_key_dpp-2_full_dpp-2_key",
    ),
    "KO_projection:dpp": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_KO_projection:dpp-1_full_dpp-1_key_dpp-2_full_dpp-2_key",
    ),
    "HMM:opp": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_HMM:opp-1_full_and_opp-1_key",
    ),
    "KO_projection:opp": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_KO_projection:opp-1_full_and_opp-1_key",
    ),
    "HMM:aap": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_HMM:aap_full_and_aap_key",
    ),
    "KO_projection:aap": (
        "H3_root_trace_gas_and_uptake",
        "requested_shorthand_unavailable_actual_endpoints_are_KO_projection:aap_full_and_aap_key",
    ),
}


def classify_exploratory_hypothesis_family(panel_name: str) -> str:
    """Assign each panel to its distinct exploratory synthesis family (never pooling unlike H6 classes)."""
    if panel_name in H6_BASELINE_DNA_REPAIR:
        return "H6_baseline_DNA_repair"
    if panel_name in H6_UV_AND_OXIDATIVE_PROTECTION:
        return "H6_UV_and_oxidative_protection"
    if panel_name in H6_OSMOLYTE_AND_SPORULATION_STRESS:
        return "H6_osmolyte_and_sporulation_stress"
    if panel_name in H3_ROOT_TRACES_AND_UPTAKE:
        return "H3_root_trace_gas_and_uptake"
    if panel_name in REQUESTED_UNAVAILABLE_ALIASES:
        return REQUESTED_UNAVAILABLE_ALIASES[panel_name][0]
    return "other_microtrait_or_community"


def pathway_interpretation_flag(panel_name: str) -> str:
    """Flag single-marker sporulation and HMM/KO trehalose route disagreement explicitly."""
    base = panel_name.replace(":mapped", "").replace(":whole", "")
    if "spoIIGA" in base:
        return "SINGLE_MARKER_NOT_FULL_SPORULATION_PATHWAY: spoIIGA encodes a single stage-II sporulation protease/sigma-processing factor, not a complete endospore formation pathway."
    if "trehalose" in base:
        return "HMM_KO_ROUTE_DISAGREEMENT_FLAGGED: microTrait HMM trehalose_biosynthesis and single-step KO projection disagree with curated two-step OtsAB/TreYZ strict trehalose synthesis."
    return "standard_audited_rule"


def compute_hill1_from_weights(W: np.ndarray) -> np.ndarray:
    """Compute Hill1 = exp(Shannon entropy in nats) per row of non-negative matrix W."""
    row_sums = W.sum(axis=1, keepdims=True)
    out = np.full(W.shape[0], np.nan, dtype="float64")
    valid = row_sums[:, 0] > 0
    if not valid.any():
        return out
    P = W[valid] / row_sums[valid]
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0.0)
    H = -np.sum(P * logP, axis=1)
    out[valid] = np.exp(H)
    return out


def verify_panel_identifier_mapping(endpoints_df: pd.DataFrame) -> dict[str, Any]:
    """Verify that every declared available panel in H6/H3 exists in `endpoints.tsv`."""
    ep_names = set(endpoints_df["name"].astype(str))
    base_ep_names = {
        nm[:-7] if nm.endswith(":mapped") else (nm[:-6] if nm.endswith(":whole") else nm)
        for nm in ep_names
    }
    declared_available = (
        H6_BASELINE_DNA_REPAIR
        + H6_UV_AND_OXIDATIVE_PROTECTION
        + H6_OSMOLYTE_AND_SPORULATION_STRESS
        + H3_ROOT_TRACES_AND_UPTAKE
    )
    missing_declared = [p for p in declared_available if p not in base_ep_names]
    if missing_declared:
        raise ValueError(f"Declared available H6/H3 panels missing from endpoints.tsv: {missing_declared}")
    return {
        "n_declared_available_h6_h3": len(declared_available),
        "n_requested_unavailable_aliases_tracked": len(REQUESTED_UNAVAILABLE_ALIASES),
        "all_declared_available_present": True,
    }


def build_systematic_effect_matrix(
    endpoints_df: pd.DataFrame,
    inla_fe_df: pd.DataFrame,
    inla_summary_json: dict[str, Any],
    paired_all_df: pd.DataFrame,
    turnover_sum_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build the systematic effect/uncertainty matrix with unique-join assertions and explicit interval bounds/status."""
    verify_panel_identifier_mapping(endpoints_df)

    # 1. Assert unique keys on frequentist paired tests (no silent duplicate aggregation)
    paired_main = paired_all_df[
        paired_all_df["omitted_trip"].isna() | (paired_all_df["omitted_trip"].astype(str) == "")
    ].copy()
    validate_primary_key(
        paired_main,
        ["panel", "contrast", "quality", "cohort", "denominator"],
        "paired_guild_tests_main",
    )
    paired_lookup = paired_main.set_index(["panel", "contrast", "quality", "cohort", "denominator"])

    # 2. Assert unique keys on INLA fixed effects
    inla_comp = inla_fe_df[inla_fe_df["term"].isin(["compartmentSurface", "compartmentRootAdjacent"])].copy()
    inla_comp["contrast"] = np.where(
        inla_comp["term"] == "compartmentSurface", "Surface-minus-Deep", "RootAdjacent-minus-Deep"
    )
    validate_primary_key(inla_comp, ["name", "spec", "contrast"], "spatial_fixed_effects_compartment")
    inla_lookup = inla_comp.set_index(["name", "spec", "contrast"])

    # Track failed/missing INLA specs from summary.json
    ep_name_to_id = dict(zip(endpoints_df["name"].astype(str), endpoints_df["endpoint"].astype(int)))
    ep_name_to_source = dict(zip(endpoints_df["name"].astype(str), endpoints_df["source"].astype(str)))
    failed_specs_set = set()
    for rel_p in inla_summary_json.get("unresolved_fit_files", []):
        # format: results/fits/075/all_FALSE_T_200_FALSE.FAILED
        parts = Path(rel_p).parts
        if len(parts) >= 4:
            ep_id = int(parts[-2])
            sp_tag = parts[-1].replace(".FAILED", "")
            failed_specs_set.add((ep_id, sp_tag))
    missing_specs_set = set()
    for m_str in inla_summary_json.get("missing_model_specs", []):
        ep_s, sp_tag = m_str.split("/", 1)
        missing_specs_set.add((int(ep_s), sp_tag))

    # 3. Assert unique keys on genus turnover summary (never use pivot_table(aggfunc='first'))
    turn_genus = turnover_sum_df[turnover_sum_df["rank"] == "genus"].copy()
    validate_primary_key(turn_genus, ["panel", "contrast", "component"], "turnover_summary_genus")
    turn_lookup = turn_genus.set_index(["panel", "contrast", "component"])

    base_panels: list[str] = []
    for nm in endpoints_df["name"].astype(str).tolist():
        bp = nm[:-7] if nm.endswith(":mapped") else (nm[:-6] if nm.endswith(":whole") else nm)
        if bp not in base_panels:
            base_panels.append(bp)

    # Append requested-but-unavailable legacy/shorthand names so none are silently omitted
    all_panels_ordered = base_panels + [
        alias for alias in REQUESTED_UNAVAILABLE_ALIASES if alias not in base_panels
    ]

    inla_spec_defs = [
        ("inla_primary", "all_TRUE_T_200_FALSE", False),
        ("inla_nodune", "all_FALSE_T_200_FALSE", False),
        ("inla_matched", "matched_TRUE_T_200_FALSE", False),
        ("inla_qual70_5", "all_TRUE_T_200_TRUE", False),
        ("inla_gaussian", "all_TRUE_gaussian_200_FALSE", False),
        ("inla_range50", "all_TRUE_T_50_FALSE", False),
        ("inla_range500", "all_TRUE_T_500_FALSE", False),
        ("inla_whole_read_primary", "all_TRUE_T_200_FALSE", True),
    ]

    freq_spec_defs = [
        ("freq_mapped_all", "all", "all", "known_catalogue_fraction"),
        ("freq_qual70_5_mapped", "completeness70_contamination5", "all", "known_catalogue_fraction"),
        ("freq_matched_cohort_mapped", "all", "matched_amplicon", "known_catalogue_fraction"),
        ("freq_whole_read_bound", "all", "all", "whole_read_bound"),
    ]

    matrix_rows: list[dict[str, Any]] = []
    for bp in all_panels_ordered:
        hyp_group = classify_exploratory_hypothesis_family(bp)
        flag_str = pathway_interpretation_flag(bp)
        is_unavailable_alias = bp in REQUESTED_UNAVAILABLE_ALIASES
        avail_status = (
            REQUESTED_UNAVAILABLE_ALIASES[bp][1]
            if is_unavailable_alias
            else "available_in_endpoints"
        )

        ep_mapped_name = f"{bp}:mapped" if f"{bp}:mapped" in ep_name_to_id else (bp if bp in ep_name_to_id else None)
        ep_whole_name = f"{bp}:whole" if f"{bp}:whole" in ep_name_to_id else None

        for contrast_name, turn_contrast in [
            ("Surface-minus-Deep", "Surface-Deep"),
            ("RootAdjacent-minus-Deep", "RootAdjacent-Deep"),
        ]:
            row_dict: dict[str, Any] = {
                "panel": bp,
                "exploratory_post_fit_hypothesis_family": hyp_group,
                "contrast": contrast_name,
                "panel_availability_status": avail_status,
                "pathway_interpretation_flag": flag_str,
            }

            # Frequentist paired site-block bootstrap columns (with actual bounds and status)
            for f_prefix, qual_k, coh_k, denom_k in freq_spec_defs:
                fkey = (bp, contrast_name, qual_k, coh_k, denom_k)
                if is_unavailable_alias:
                    row_dict[f"{f_prefix}_status"] = "requested_panel_unavailable"
                    row_dict[f"{f_prefix}_mean_pp"] = np.nan
                    row_dict[f"{f_prefix}_ci95_low"] = np.nan
                    row_dict[f"{f_prefix}_ci95_high"] = np.nan
                    if f_prefix == "freq_mapped_all":
                        row_dict["freq_paired_sites"] = np.nan
                        row_dict["freq_paired_visits"] = np.nan
                        row_dict["freq_mapped_wilcoxon_p"] = np.nan
                        row_dict["freq_mapped_q_bh_all41"] = np.nan
                elif fkey in paired_lookup.index:
                    pr = paired_lookup.loc[fkey]
                    row_dict[f"{f_prefix}_status"] = "completed"
                    row_dict[f"{f_prefix}_mean_pp"] = float(pr["mean_pp"])
                    row_dict[f"{f_prefix}_ci95_low"] = float(pr["ci_low"])
                    row_dict[f"{f_prefix}_ci95_high"] = float(pr["ci_high"])
                    if f_prefix == "freq_mapped_all":
                        row_dict["freq_paired_sites"] = int(pr["sites"])
                        row_dict["freq_paired_visits"] = int(pr["visits"])
                        row_dict["freq_mapped_wilcoxon_p"] = float(pr["p"])
                        row_dict["freq_mapped_q_bh_all41"] = float(pr["q_bh"])
                else:
                    row_dict[f"{f_prefix}_status"] = "not_applicable_non_mag_or_curated_endpoint"
                    row_dict[f"{f_prefix}_mean_pp"] = np.nan
                    row_dict[f"{f_prefix}_ci95_low"] = np.nan
                    row_dict[f"{f_prefix}_ci95_high"] = np.nan
                    if f_prefix == "freq_mapped_all":
                        row_dict["freq_paired_sites"] = np.nan
                        row_dict["freq_paired_visits"] = np.nan
                        row_dict["freq_mapped_wilcoxon_p"] = np.nan
                        row_dict["freq_mapped_q_bh_all41"] = np.nan

            # Never back-calculate SE/90% bounds from asymmetric 95% bootstrap percentiles (Item 4)
            row_dict["freq_equivalence_status_default"] = (
                "requested_panel_unavailable"
                if is_unavailable_alias
                else (
                    "no_equivalence_claim_unjustified_margin"
                    if row_dict["freq_mapped_all_status"] == "completed"
                    else "not_applicable_non_mag_endpoint"
                )
            )

            # Bayesian INLA spatial hierarchy columns (with actual bounds and explicit missing/failed/numerical status)
            n_fitted_map = 0
            n_excl_zero_map = 0
            n_num_flag_map = 0
            for i_prefix, sp_tag, is_whole in inla_spec_defs:
                target_ep = ep_whole_name if is_whole else ep_mapped_name
                if is_unavailable_alias or target_ep is None:
                    st_val = "requested_panel_unavailable" if is_unavailable_alias else "not_applicable_no_whole_endpoint"
                    row_dict[f"{i_prefix}_status"] = st_val
                    row_dict[f"{i_prefix}_mean_orig"] = np.nan
                    row_dict[f"{i_prefix}_sd_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_low_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_high_orig"] = np.nan
                    row_dict[f"{i_prefix}_excludes_zero"] = st_val
                    row_dict[f"{i_prefix}_numerical_review_flag"] = st_val
                    continue

                ep_id = ep_name_to_id[target_ep]
                ep_src = ep_name_to_source[target_ep]
                if sp_tag.endswith("_TRUE") and ep_src != "MAG":
                    st_val = "not_applicable_non_mag_endpoint"
                    row_dict[f"{i_prefix}_status"] = st_val
                    row_dict[f"{i_prefix}_mean_orig"] = np.nan
                    row_dict[f"{i_prefix}_sd_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_low_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_high_orig"] = np.nan
                    row_dict[f"{i_prefix}_excludes_zero"] = st_val
                    row_dict[f"{i_prefix}_numerical_review_flag"] = st_val
                    continue

                ikey = (target_ep, sp_tag, contrast_name)
                if ikey in inla_lookup.index:
                    ir = inla_lookup.loc[ikey]
                    num_req = bool(ir["numerical_review_required"])
                    excl_z = bool(ir["interval_excludes_zero"])
                    st_val = "completed_numerical_review_required" if num_req else "completed"
                    row_dict[f"{i_prefix}_status"] = st_val
                    row_dict[f"{i_prefix}_mean_orig"] = float(ir["mean_original_scale"])
                    row_dict[f"{i_prefix}_sd_orig"] = float(ir["sd_original_scale"])
                    row_dict[f"{i_prefix}_cred95_low_orig"] = float(ir["0.025quant_original_scale"])
                    row_dict[f"{i_prefix}_cred95_high_orig"] = float(ir["0.975quant_original_scale"])
                    row_dict[f"{i_prefix}_excludes_zero"] = "True" if excl_z else "False"
                    row_dict[f"{i_prefix}_numerical_review_flag"] = "True" if num_req else "False"
                    if not is_whole:
                        n_fitted_map += 1
                        if excl_z:
                            n_excl_zero_map += 1
                        if num_req:
                            n_num_flag_map += 1
                else:
                    if (ep_id, sp_tag) in failed_specs_set:
                        st_val = "failed_inla_fit"
                    elif (ep_id, sp_tag) in missing_specs_set:
                        st_val = "missing_inla_spec"
                    else:
                        st_val = "missing_inla_spec"
                    row_dict[f"{i_prefix}_status"] = st_val
                    row_dict[f"{i_prefix}_mean_orig"] = np.nan
                    row_dict[f"{i_prefix}_sd_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_low_orig"] = np.nan
                    row_dict[f"{i_prefix}_cred95_high_orig"] = np.nan
                    row_dict[f"{i_prefix}_excludes_zero"] = st_val
                    row_dict[f"{i_prefix}_numerical_review_flag"] = st_val

            row_dict["inla_n_mapped_specs_fitted"] = n_fitted_map if not is_unavailable_alias else np.nan
            row_dict["inla_n_mapped_specs_cred95_excludes_zero"] = n_excl_zero_map if n_fitted_map > 0 else np.nan
            row_dict["inla_n_mapped_specs_numerical_flag"] = n_num_flag_map if n_fitted_map > 0 else np.nan
            row_dict["bayesian_rope_status_default"] = (
                "requested_panel_unavailable"
                if is_unavailable_alias
                else (
                    "no_equivalence_claim_unjustified_rope"
                    if str(row_dict["inla_primary_status"]).startswith("completed")
                    else str(row_dict["inla_primary_status"])
                )
            )

            # Symmetric genus turnover decomposition (between_shared, within_shared, exclusive)
            for comp_k, col_label in [
                ("total", "total"),
                ("between", "between_shared"),
                ("within", "within_shared"),
                ("exclusive", "exclusive"),
            ]:
                tkey = (bp, turn_contrast, comp_k)
                if tkey in turn_lookup.index:
                    tr = turn_lookup.loc[tkey]
                    row_dict[f"turnover_genus_{col_label}_status"] = "completed"
                    row_dict[f"turnover_genus_{col_label}_pp"] = float(tr["mean"])
                    row_dict[f"turnover_genus_{col_label}_ci95_low"] = float(tr["low"])
                    row_dict[f"turnover_genus_{col_label}_ci95_high"] = float(tr["high"])
                else:
                    row_dict[f"turnover_genus_{col_label}_status"] = (
                        "requested_panel_unavailable" if is_unavailable_alias else "not_applicable_non_mag_trait"
                    )
                    row_dict[f"turnover_genus_{col_label}_pp"] = np.nan
                    row_dict[f"turnover_genus_{col_label}_ci95_low"] = np.nan
                    row_dict[f"turnover_genus_{col_label}_ci95_high"] = np.nan

            matrix_rows.append(row_dict)

    full_matrix_df = pd.DataFrame(matrix_rows)

    # Compute BH q-values within each separate exploratory post-fit family (`exploratory_post_fit_hypothesis_family x contrast`)
    full_matrix_df["bh_q_within_exploratory_family"] = np.nan
    for (hf, ct), grp in full_matrix_df.groupby(["exploratory_post_fit_hypothesis_family", "contrast"], sort=False):
        valid_p = grp["freq_mapped_wilcoxon_p"].notna()
        if valid_p.any():
            idx_v = grp.index[valid_p]
            full_matrix_df.loc[idx_v, "bh_q_within_exploratory_family"] = bh_fdr(
                full_matrix_df.loc[idx_v, "freq_mapped_wilcoxon_p"].to_numpy(dtype=float)
            )

    # Guard against typo creating all-NaN matrix with completed status
    n_completed_inla = int(full_matrix_df["inla_primary_status"].astype(str).str.startswith("completed").sum())
    if n_completed_inla == 0:
        raise ValueError("Systematic effect matrix has 0 completed primary INLA rows; check endpoint/spec identifiers")

    return full_matrix_df


def main() -> None:
    P_RES.mkdir(parents=True, exist_ok=True)
    print("[Stage 4] Loading design, 1,800 INLA fits, paired guild tests, and turnover summaries...", flush=True)

    design_df = pd.read_csv(require_file(P_ECO_INT / "results/design.tsv", "design.tsv"), sep="\t")
    validate_primary_key(design_df, ["sample_id"], "design.tsv")
    validate_audited_compartment(design_df, "compartment")
    assert len(design_df) == 546

    endpoints_df = pd.read_csv(require_file(P_ECO_INT / "results/endpoints.tsv", "endpoints.tsv"), sep="\t")
    inla_fe_df = pd.read_csv(require_file(P_ECO_INT / "results/spatial_fixed_effects.tsv", "spatial_fixed_effects.tsv"), sep="\t")
    inla_diag_df = pd.read_csv(require_file(P_ECO_INT / "results/spatial_fit_diagnostics.tsv", "spatial_fit_diagnostics.tsv"), sep="\t")
    inla_summary_json = json.loads(require_file(P_ECO_INT / "results/summary.json", "summary.json").read_text(encoding="utf-8"))
    turnover_sum_df = pd.read_csv(require_file(P_ECO_INT / "results/turnover_summary.tsv", "turnover_summary.tsv"), sep="\t")

    paired_hmm_df = pd.read_csv(require_file(P_PROT_INT / "results/paired_guild_tests_hmm.tsv", "paired_guild_tests_hmm.tsv"), sep="\t")
    paired_ko_df = pd.read_csv(require_file(P_PROT_INT / "results/paired_guild_tests_ko.tsv", "paired_guild_tests_ko.tsv"), sep="\t")
    paired_all_df = pd.concat([paired_hmm_df, paired_ko_df], ignore_index=True)

    audit_hmm_df = pd.read_csv(require_file(P_PROT_INT / "results/panel_audit_hmm.tsv", "panel_audit_hmm.tsv"), sep="\t")

    # =========================================================================
    # 1. Systematic Effect & Uncertainty Matrix across all 41 microTrait + Legacy/Curated/Typed/Read Panels
    # =========================================================================
    print("[Stage 4] Building systematic H6/H3 and full 41-rule effect/uncertainty matrix...", flush=True)
    full_matrix_df = build_systematic_effect_matrix(
        endpoints_df, inla_fe_df, inla_summary_json, paired_all_df, turnover_sum_df
    )
    full_matrix_df.to_csv(P_RES / "task_d_all_41_microtrait_and_curated_matrix.tsv", sep="\t", index=False)

    h6_h3_sub = full_matrix_df[
        (
            full_matrix_df["exploratory_post_fit_hypothesis_family"].isin(
                [
                    "H6_baseline_DNA_repair",
                    "H6_UV_and_oxidative_protection",
                    "H6_osmolyte_and_sporulation_stress",
                ]
            )
            & (full_matrix_df["contrast"] == "Surface-minus-Deep")
        )
        | (
            (full_matrix_df["exploratory_post_fit_hypothesis_family"] == "H3_root_trace_gas_and_uptake")
            & (full_matrix_df["contrast"] == "RootAdjacent-minus-Deep")
        )
    ].copy()
    h6_h3_sub.to_csv(P_RES / "task_d_h6_h3_systematic_matrix.tsv", sep="\t", index=False)

    h6_surf = full_matrix_df[
        full_matrix_df["exploratory_post_fit_hypothesis_family"].isin(
            [
                "H6_baseline_DNA_repair",
                "H6_UV_and_oxidative_protection",
                "H6_osmolyte_and_sporulation_stress",
            ]
        )
        & (full_matrix_df["contrast"] == "Surface-minus-Deep")
    ].copy()
    h6_surf.to_csv(P_RES / "task_d_h6_baseline_vs_protection_surface_deep.tsv", sep="\t", index=False)

    # =========================================================================
    # 2. All-Read vs Mapped-Community Turnover Verification (Review Item 3)
    # =========================================================================
    print("[Stage 4] Verifying CoverM fraction units and all-read vs mapped-community symmetric turnover...", flush=True)
    M_DIR = P_ECO_MAG / "inputs/ibex/rebuild_20260729/picrust2_current_20260803/mag_ko"
    mag_index_df = pd.read_csv(require_file(M_DIR / "mag_index.tsv", "mag_index.tsv"), sep="\t").sort_values("mag_index")
    qual_all_df = (
        pd.read_csv(require_file(P_ARCH / "inputs/genome_quality.tsv", "genome_quality.tsv"), sep="\t")
        .set_index("genome")
        .reindex(mag_index_df["mag_id"])
    )
    AA_raw = sparse.load_npz(require_file(M_DIR / "sample_by_mag_coverm_relative_abundance.npz", "coverm npz")).toarray()[
        design_df["sample_index"].to_numpy()
    ]

    # Verify explicitly that CoverM matrix stores fractions (0..1) rather than percentages (0..100) before multiplying by 100
    row_mapped_fractions = AA_raw.sum(axis=1)
    assert np.all((row_mapped_fractions > 0.005) & (row_mapped_fractions <= 1.0)), (
        f"CoverM per-row mapped sums outside declared fraction range (0.005, 1.0]: "
        f"min={row_mapped_fractions.min()}, max={row_mapped_fractions.max()}"
    )
    cat_mapped_pct = design_df["catalogue_mapped_percent"].to_numpy(dtype="float64")
    max_pct_diff = float(np.max(np.abs(100.0 * row_mapped_fractions - cat_mapped_pct)))
    assert max_pct_diff < 1e-6, f"100 * AA_raw.sum(1) differs from catalogue_mapped_percent by {max_pct_diff}"

    G_traits = (
        pd.read_csv(require_file(P_ECO_INT / "results/genome_traits.tsv", "genome_traits.tsv"), sep="\t")
        .set_index("genome")
        .reindex(mag_index_df["mag_id"])
    )
    genera_all = [
        next((x for x in str(tax).split(";") if x.startswith("g__")), "g__")
        for tax in qual_all_df["gtdb_classification"]
    ]
    genomes_all = mag_index_df["mag_id"].tolist()

    def _mean_visit_weights(W_mat: np.ndarray) -> pd.DataFrame:
        df_w = pd.DataFrame(W_mat)
        df_w[["site", "trip", "compartment"]] = design_df[["site", "trip", "compartment"]].to_numpy()
        return df_w.groupby(["site", "trip", "compartment"]).mean()

    V_allread = _mean_visit_weights(AA_raw)
    V_mapped = _mean_visit_weights(AA_raw / row_mapped_fractions[:, None])

    turn_genus_lookup = turnover_sum_df[turnover_sum_df["rank"] == "genus"].set_index(["panel", "contrast", "component"])
    key_turnover_traits = [
        "legacy:base_excision_repair",
        "legacy:mismatch_repair",
        "legacy:nucleotide_excision_repair",
        "legacy:SOS_translesion",
        "legacy:nhej",
        "legacy:photolyase",
        "legacy:carotenoid",
        "legacy:dna_protection_oxidative",
        "HMM:ahpCF",
        "KO_projection:ahpCF",
        "HMM:sodA",
        "KO_projection:sodA",
        "HMM:opu_and",
        "KO_projection:opu_and",
        "HMM:glycine betaine/proline transport system",
        "KO_projection:glycine betaine/proline transport system",
        "curated:ectoine_marker",
        "curated:trehalose_strict",
        "typed:[NiFe]_Group_1l",
        "HMM:dpp-1_full",
        "KO_projection:dpp-1_full",
    ]

    denom_turn_rows = []
    for t_col in key_turnover_traits:
        if t_col not in G_traits.columns:
            continue
        t_vec = G_traits[t_col].to_numpy(dtype="float64")
        if not np.isfinite(t_vec).all():
            continue
        for c_target, turn_ct in [("Surface", "Surface-Deep"), ("RootAdjacent", "RootAdjacent-Deep")]:
            for denom_mode, V_use, norm_flag in [
                ("all_reads_unnormalized_weights", V_allread, False),
                ("mapped_community_normalized_weights", V_mapped, True),
            ]:
                left_v = V_use.xs(c_target, level=2)
                right_v = V_use.xs("Deep", level=2)
                shared_visits = left_v.index.intersection(right_v.index)
                site_res: dict[str, list[dict[str, float]]] = {}
                for site, trip in shared_visits:
                    wa = left_v.loc[(site, trip)].to_numpy(dtype="float64")
                    wb = right_v.loc[(site, trip)].to_numpy(dtype="float64")
                    dec = decompose_trait_turnover(
                        wa, wb, t_vec, genomes_all, genera_all, normalize_to_mapped=norm_flag
                    )
                    site_res.setdefault(str(site), []).append(dec)
                site_means = [
                    {
                        k: 100.0 * float(np.mean([r[k] for r in r_list]))
                        for k in ("total_diff", "between_shared", "within_shared", "exclusive", "unclassified_mass")
                    }
                    for r_list in site_res.values()
                ]
                m_tot = float(np.mean([x["total_diff"] for x in site_means]))
                m_bet = float(np.mean([x["between_shared"] for x in site_means]))
                m_wit = float(np.mean([x["within_shared"] for x in site_means]))
                m_exc = float(np.mean([x["exclusive"] for x in site_means]))

                # Verify exact agreement with existing `turnover_summary.tsv` when `norm_flag == True`
                max_diff_existing = np.nan
                if norm_flag and (t_col, turn_ct, "total") in turn_genus_lookup.index:
                    ex_tot = float(turn_genus_lookup.loc[(t_col, turn_ct, "total"), "mean"])
                    ex_bet = float(turn_genus_lookup.loc[(t_col, turn_ct, "between"), "mean"])
                    ex_wit = float(turn_genus_lookup.loc[(t_col, turn_ct, "within"), "mean"])
                    ex_exc = float(turn_genus_lookup.loc[(t_col, turn_ct, "exclusive"), "mean"])
                    max_diff_existing = float(
                        max(abs(m_tot - ex_tot), abs(m_bet - ex_bet), abs(m_wit - ex_wit), abs(m_exc - ex_exc))
                    )
                    assert max_diff_existing < 1e-6, (
                        f"Mapped turnover mismatch vs existing turnover_summary.tsv for {t_col} {turn_ct}: {max_diff_existing}"
                    )

                denom_turn_rows.append(
                    {
                        "panel": t_col,
                        "contrast": f"{c_target}-minus-Deep",
                        "denominator_mode": denom_mode,
                        "coverm_input_unit_verified": "fraction_0_to_1_scaled_by_100_to_percentage_points",
                        "normalize_to_mapped": norm_flag,
                        "n_paired_sites": len(site_means),
                        "mean_total_diff_pp": m_tot,
                        "mean_between_shared_pp": m_bet,
                        "mean_within_shared_pp": m_wit,
                        "mean_exclusive_pp": m_exc,
                        "mean_unclassified_mass_pp": float(np.mean([x["unclassified_mass"] for x in site_means])),
                        "max_abs_diff_vs_existing_turnover_summary_pp": max_diff_existing,
                    }
                )
    pd.DataFrame(denom_turn_rows).to_csv(P_RES / "task_d_turnover_allread_vs_mapped_comparison.tsv", sep="\t", index=False)

    # =========================================================================
    # 3. Spatially Blocked Joint Contrast of `megabarchans` vs `uruq` (`bacteria_hill1` vs `KO_hq_mq_hill1`)
    #    Preserving Exact Covariate Design, Site-Balanced Fixed Base Weights, and Quality/Detection Sensitivities
    # =========================================================================
    print("[Stage 4] Computing matched joint log contrast of megabarchans vs uruq for bacteria_hill1 vs KO_hq_mq_hill1...", flush=True)
    ko_copy_mat = sparse.load_npz(require_file(M_DIR / "mag_by_ko_copy.npz", "mag_by_ko_copy.npz")).toarray()
    ko_pres_mat = (ko_copy_mat > 0).astype("float64")
    ko_sbest_mat = sparse.load_npz(
        require_file(M_DIR / "mag_by_ko_copy_single_best_sensitivity.npz", "mag_by_ko_copy_single_best_sensitivity.npz")
    ).toarray()

    is_bact = qual_all_df["gtdb_classification"].astype(str).str.startswith("d__Bacteria").to_numpy(dtype=bool)
    is_q70_5 = ((qual_all_df["checkm_completeness"] >= 70.0) & (qual_all_df["contamination"] <= 5.0)).to_numpy(dtype=bool)
    is_q90_5 = ((qual_all_df["checkm_completeness"] > 90.0) & (qual_all_df["contamination"] < 5.0)).to_numpy(dtype=bool)
    comp_frac = np.clip(qual_all_df["checkm_completeness"].to_numpy(dtype="float64") / 100.0, 0.50, 1.00)

    div_df = design_df[
        [
            "sample_id",
            "site",
            "trip",
            "compartment",
            "matched_amplicon",
            "edgell_dune_type",
            "transect_km",
            "post_qc_reads",
            "bacteria_hill1",
            "KO_hq_mq_hill1",
            "KO_hq_hill1",
            "KO_hq_mq_completeness_corrected_hill1",
            "KO_hq_mq_single_best_ko_hill1",
            "shannon_nats",
        ]
    ].copy()
    div_df["log_reads"] = np.log1p(pd.to_numeric(div_df["post_qc_reads"], errors="coerce"))
    div_df["bact_hill1_all_mags"] = pd.to_numeric(div_df["bacteria_hill1"], errors="coerce")
    div_df["bact_hill1_comp70_contam5"] = compute_hill1_from_weights(AA_raw[:, is_bact & is_q70_5])
    div_df["bact_hill1_comp90_contam5"] = compute_hill1_from_weights(AA_raw[:, is_bact & is_q90_5])
    div_df["amplicon_16S_hill1"] = np.exp(pd.to_numeric(div_df["shannon_nats"], errors="coerce"))

    div_df["ko_hill1_copy_all_mags"] = pd.to_numeric(div_df["KO_hq_mq_hill1"], errors="coerce")
    div_df["ko_hill1_presence_binary"] = compute_hill1_from_weights(AA_raw @ ko_pres_mat)
    div_df["ko_hill1_single_best_hit"] = compute_hill1_from_weights(AA_raw @ ko_sbest_mat)
    div_df["ko_hill1_comp70_contam5"] = compute_hill1_from_weights(AA_raw[:, is_q70_5] @ ko_copy_mat[is_q70_5])
    div_df["ko_hill1_completeness_corrected"] = compute_hill1_from_weights((AA_raw / comp_frac[None, :]) @ ko_copy_mat)

    # Spatial 10-block definition matching `ecology-integration-20261108/code/prepare.py`
    site_order = div_df.groupby("site")["transect_km"].mean().sort_values().index
    site_to_block = {s: f"spatial_block_{i:02d}" for i, grp in enumerate(np.array_split(site_order, 10)) for s in grp}
    div_df["spatial_block"] = div_df["site"].map(site_to_block)

    # Extract existing marginal INLA `edgell_dune_typemegabarchans` fixed effects for reference
    inla_dune_fe = inla_fe_df[
        inla_fe_df["name"].isin(["bacteria_hill1", "KO_hq_mq_hill1"])
        & (inla_fe_df["term"] == "edgell_dune_typemegabarchans")
    ].set_index(["name", "spec"])

    pair_sensitivities = [
        ("primary_all_mags_copy", "bact_hill1_all_mags", "ko_hill1_copy_all_mags", "all_TRUE_T_200_FALSE", "matched_TRUE_T_200_FALSE"),
        ("quality_comp70_contam5", "bact_hill1_comp70_contam5", "ko_hill1_comp70_contam5", None, None),
        ("quality_comp90_contam5_vs_comp70_ko", "bact_hill1_comp90_contam5", "ko_hill1_comp70_contam5", None, None),
        ("detection_binary_ko_presence", "bact_hill1_all_mags", "ko_hill1_presence_binary", None, None),
        ("detection_single_best_ko_hit", "bact_hill1_all_mags", "ko_hill1_single_best_hit", None, None),
        ("completeness_corrected_ko_weights", "bact_hill1_all_mags", "ko_hill1_completeness_corrected", None, None),
        ("amplicon_16S_vs_ko_copy", "amplicon_16S_hill1", "ko_hill1_copy_all_mags", None, None),
    ]

    joint_rows: list[dict[str, Any]] = []
    for cohort_label, coh_mask, inla_spec_col in [
        ("matched_amplicon_441_modeled_cohort", div_df["matched_amplicon"].astype(str).str.lower().isin(["true", "1"]).to_numpy(), 4),
        ("all_545_modeled_shotgun_cohort", np.ones(len(div_df), dtype=bool), 3),
    ]:
        sub_coh = div_df[coh_mask].copy()
        for sens_name, tax_col, func_col, sp_all, sp_mat in pair_sensitivities:
            inla_sp = sp_all if inla_spec_col == 3 else sp_mat
            inla_bact = inla_dune_fe.loc[("bacteria_hill1", inla_sp)] if (inla_sp and ("bacteria_hill1", inla_sp) in inla_dune_fe.index) else None
            inla_ko = inla_dune_fe.loc[("KO_hq_mq_hill1", inla_sp)] if (inla_sp and ("KO_hq_mq_hill1", inla_sp) in inla_dune_fe.index) else None

            for w_mode in ("sample_weighted", "site_balanced", "site_averaged"):
                res = estimate_matched_joint_log_contrast(
                    sub_coh,
                    tax_col=tax_col,
                    func_col=func_col,
                    dune_col="edgell_dune_type",
                    ref_dune="uruq",
                    target_dune="megabarchans",
                    compartment_col="compartment",
                    ref_compartment="Deep",
                    trip_col="trip",
                    ref_trip="trip1",
                    log_reads_col="log_reads",
                    site_col="site",
                    block_col="spatial_block",
                    weighting=w_mode,
                    n_boot=1000,
                    seed=20261108,
                )

                lo_d = float(res.get("delta_log_ci95_low", np.nan))
                hi_d = float(res.get("delta_log_ci95_high", np.nan))
                if res["status"] != "completed":
                    d_status = str(res["status"])
                elif lo_d > 0.0 or hi_d < 0.0:
                    d_status = "conditional_bootstrap_ci95_excludes_zero_not_replacement_for_INLA"
                else:
                    d_status = "conditional_bootstrap_ci95_crosses_zero_not_evidence_of_equivalence"

                joint_rows.append(
                    {
                        "cohort": cohort_label,
                        "estimand": "edgell_dune_type_megabarchans_minus_uruq_adjusted_for_compartment_trip_log_reads",
                        "sensitivity_specification": sens_name,
                        "weighting_scheme": w_mode,
                        "taxonomic_endpoint": tax_col,
                        "functional_endpoint": func_col,
                        "status": str(res["status"]),
                        "resampling_design": str(res.get("resampling_design", "")),
                        "n_matched_units": int(res.get("n_samples", 0)),
                        "n_sites": int(res.get("n_sites", 0)),
                        "n_megabarchan_units": int(res.get("n_megabarchan_samples", 0)),
                        "n_uruq_units": int(res.get("n_uruq_samples", 0)),
                        "n_spatial_blocks": int(res.get("n_blocks", 0)),
                        "n_boot_requested": int(res.get("n_boot_requested", 1000)),
                        "n_valid_replicates": int(res.get("n_valid_replicates", 0)),
                        "n_invalid_rank_deficient_or_zero_cell_replicates": int(res.get("n_invalid_replicates", 0)),
                        # Matched joint covariate-adjusted log contrasts (multiplicative natural-log scale)
                        "beta_log_tax_megabarchans_vs_uruq": float(res.get("beta_log_tax", np.nan)),
                        "beta_log_tax_se": float(res.get("beta_log_tax_se", np.nan)),
                        "beta_log_tax_ci95_low": float(res.get("beta_log_tax_ci95_low", np.nan)),
                        "beta_log_tax_ci95_high": float(res.get("beta_log_tax_ci95_high", np.nan)),
                        "beta_log_func_megabarchans_vs_uruq": float(res.get("beta_log_func", np.nan)),
                        "beta_log_func_se": float(res.get("beta_log_func_se", np.nan)),
                        "beta_log_func_ci95_low": float(res.get("beta_log_func_ci95_low", np.nan)),
                        "beta_log_func_ci95_high": float(res.get("beta_log_func_ci95_high", np.nan)),
                        "delta_log_tax_minus_func": float(res.get("delta_log_tax_minus_func", np.nan)),
                        "delta_log_se": float(res.get("delta_log_se", np.nan)),
                        "delta_log_ci95_low": lo_d,
                        "delta_log_ci95_high": hi_d,
                        "fold_ratio_tax_over_func": float(res.get("fold_ratio_tax_over_func", np.nan)),
                        "fold_ratio_ci95_low": float(res.get("fold_ratio_ci95_low", np.nan)),
                        "fold_ratio_ci95_high": float(res.get("fold_ratio_ci95_high", np.nan)),
                        "joint_contrast_inference_status": d_status,
                        # Reference marginal INLA SPDE fixed effects (where fitted)
                        "inla_marginal_log_tax_mean": float(inla_bact["mean_original_scale"]) if inla_bact is not None else np.nan,
                        "inla_marginal_log_tax_cred95_low": float(inla_bact["0.025quant_original_scale"]) if inla_bact is not None else np.nan,
                        "inla_marginal_log_tax_cred95_high": float(inla_bact["0.975quant_original_scale"]) if inla_bact is not None else np.nan,
                        "inla_marginal_log_func_mean": float(inla_ko["mean_original_scale"]) if inla_ko is not None else np.nan,
                        "inla_marginal_log_func_cred95_low": float(inla_ko["0.025quant_original_scale"]) if inla_ko is not None else np.nan,
                        "inla_marginal_log_func_cred95_high": float(inla_ko["0.975quant_original_scale"]) if inla_ko is not None else np.nan,
                        "inla_joint_posterior_contrast_status": "marginal_inla_fits_available_joint_bivariate_inla_not_fitted",
                    }
                )

    joint_df = pd.DataFrame(joint_rows)
    joint_df.to_csv(P_RES / "task_d_eastern_joint_diversity_contrast.tsv", sep="\t", index=False)

    # Save frozen panel definitions, multiplicity families, and chemistry-model distinction manifest
    panels_json_path = require_file(P_PROT_INT / "inputs/panels.json", "panels.json")
    panels_spec = json.loads(panels_json_path.read_text(encoding="utf-8"))

    manifest_payload = {
        "stage": "stage4_ecology_synthesis",
        "status": "completed",
        "date_frozen": "2026-11-09",
        "exploratory_synthesis_disclaimer": (
            "The H6 and H3 synthesis comparisons and panel groupings were structured after initial fits; "
            "they are an auditable exploratory synthesis, not a pre-specified or omnibus confirmatory test. "
            "H6 is separated into three distinct biological panels (baseline DNA repair, UV/oxidative protection, "
            "and osmolyte/sporulation stress) and never combined into an omnibus claim."
        ),
        "inla_vs_rund_chemistry_distinction": {
            "spatial_hierarchy_adaptation": (
                "Reuses the 1,800 INLA Matérn SPDE spatial + site random-effect fits across 207 endpoints "
                "(all vs matched cohort, dune-zone covariate FALSE/TRUE, Student-t vs Gaussian likelihood, "
                "prior range 50/200/500 km, all MAGs vs completeness>=70/contamination<=5)."
            ),
            "chemistry_models_status": "BLOCKED_UNAUDITED_JOINS",
            "chemistry_models_note": (
                "Rund's original soil-chemistry regressions are not run on shotgun functional endpoints because "
                "audited specimen-to-chemistry joins are unavailable for the full shotgun cohort, and XRF total "
                "elements are not measured porewater salinity. Environment/landform/dune-zone classifications "
                "are never treated as independent soil chemistry measurements."
            ),
        },
        "explicit_pathway_flags": {
            "spoIIGA": "SINGLE_MARKER_NOT_FULL_SPORULATION_PATHWAY",
            "trehalose_biosynthesis": "HMM_KO_ROUTE_DISAGREEMENT_FLAGGED",
        },
        "equivalence_and_interval_policy": (
            "Bayesian INLA posterior credible intervals and frequentist site-block bootstrap confidence intervals "
            "are reported in separate columns and never conflated. Asymmetric 95% bootstrap percentile intervals "
            "are never back-calculated into synthetic SEs or 90% TOST bounds. With no externally calibrated "
            "equivalence margin, primary equivalence status is 'no_equivalence_claim_unjustified_margin' / "
            "'no_equivalence_claim_unjustified_rope', and any interval crossing zero is never claimed as evidence of equivalence."
        ),
        "exploratory_post_fit_hypothesis_families": {
            "H6_baseline_DNA_repair": H6_BASELINE_DNA_REPAIR,
            "H6_UV_and_oxidative_protection": H6_UV_AND_OXIDATIVE_PROTECTION,
            "H6_osmolyte_and_sporulation_stress": H6_OSMOLYTE_AND_SPORULATION_STRESS,
            "H3_root_trace_gas_and_uptake": H3_ROOT_TRACES_AND_UPTAKE,
            "requested_unavailable_legacy_aliases": REQUESTED_UNAVAILABLE_ALIASES,
            "full_microtrait_41_panel": list(audit_hmm_df["panel"].astype(str)),
        },
        "panels_definition": panels_spec,
        "n_total_endpoints_inla": int(len(endpoints_df)),
        "n_total_inla_fits": int(len(inla_diag_df)),
        "n_numerical_review_inla_fits": int(inla_diag_df["numerical_review_required"].sum()),
    }
    (P_RES / "task_d_panel_and_multiplicity_manifest.json").write_text(json.dumps(manifest_payload, indent=2))
    print("[Stage 4] Completed successfully!", flush=True)


if __name__ == "__main__":
    main()
