#!/usr/bin/env python3
"""Canonical GTDB r232 Taxonomy Sensitivity Overlay, Reconciliation, Stage 1/1b Replay, and 3-Guild Genus Turnover Sensitivity.

Implements all requirements from `coordinated-revision-20261108/canonical-taxonomy-sensitivity-design.txt`:
1. 1:1 reconciliation of all 5,438 genomes in `genome_quality.tsv` and 237 genomes in `thermoplasmatota.tsv`
   against the data foundation `refined_representatives.tsv` (`genome == mag`), with strict validation of
   genome ID set equality, duplicate detection, domain invariance, completeness/contamination invariance,
   and stale taxonomy label substitution detection.
2. Read-only overlay input tree under `/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/overlay_arch`
   replacing `gtdb_classification` with canonical `gtdb_r232_classification` while retaining `legacy_gtdb_classification`
   as an explicit column and leaving all frozen source datasets untouched.
3. Exact reuse of `stage1_audit_and_export.py` and `stage1b_refined_fasta_and_corrections.py` algorithms by
   explicitly overriding module path globals (`P_ARCH`, `P_OUT_ROOT`, `P_IN`, `P_RES`), verifying that
   CDS/ORF coordinate and translation audits are invariant while recording canonical vs legacy lineage
   recurrence and Task B export OR sensitivity.
4. Exact recomputation of symmetric genus turnover decomposition (`total_diff`, `between_shared`,
   `within_shared`, `exclusive`, `unclassified_mass`) with 1,000-replicate site-block bootstrap CIs for
   BetL/betaine transport, AhpCF, and `[NiFe]` Group 1l under canonical vs legacy genus labels.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    block_bootstrap_stat,
    decompose_trait_turnover,
    normalize_genus_groups,
    require_file,
    sha256_file,
    validate_audited_compartment,
    validate_primary_key,
)
try:
    from stage1_audit_and_export import parse_gtdb_lineage
except ImportError:
    def parse_gtdb_lineage(tax_str: str) -> dict[str, str]:
        ranks = {"d": "", "p": "", "c": "", "o": "", "f": "", "g": "", "s": ""}
        if not isinstance(tax_str, str):
            return ranks
        for part in tax_str.split(";"):
            part = part.strip()
            if "__" in part:
                k, v = part.split("__", 1)
                if k in ranks:
                    ranks[k] = f"{k}__{v}" if v else f"{k}__"
        return ranks

DEFAULT_ORIG_ARCH = Path(os.environ.get("EQ_ORIG_ARCH_ROOT", "/scratch/project/k10018/eq-archaeal-function-20261106"))
DEFAULT_ORIG_SCI = Path(os.environ.get("EQ_ORIG_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
DEFAULT_CANON_OUT = Path(
    os.environ.get("EQ_CANONICAL_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261109-canonical-taxonomy")
)
DEFAULT_ECO_INT = Path(os.environ.get("EQ_ECO_INT_ROOT", "/scratch/project/k10018/eq-ecology-integration-20261108"))
DEFAULT_ECO_MAG = Path(os.environ.get("EQ_ECO_MAG_ROOT", "/scratch/project/k10018/eq-ecology-mag-20261105"))

RANK_KEYS = [
    ("d", "domain"),
    ("p", "phylum"),
    ("c", "class"),
    ("o", "order"),
    ("f", "family"),
    ("g", "genus"),
    ("s", "species"),
]

PRIMARY_REQUIRED_3GUILD_PANELS = [
    "HMM:glycine betaine transporter",
    "HMM:ahpCF",
    "typed:[NiFe]_Group_1l",
]

OPTIONAL_ADDITIONAL_TURNOVER_PANELS = [
    "KO_projection:glycine betaine transporter",
    "HMM:opu_and",
    "KO_projection:opu_and",
    "HMM:glycine betaine/proline transport system",
    "KO_projection:glycine betaine/proline transport system",
    "KO_projection:ahpCF",
]

THREE_GUILD_TURNOVER_PANELS = PRIMARY_REQUIRED_3GUILD_PANELS + OPTIONAL_ADDITIONAL_TURNOVER_PANELS


def reconcile_and_validate_canonical_taxonomy(
    qual_df: pd.DataFrame,
    thermo_df: pd.DataFrame,
    canonical_df: pd.DataFrame,
    *,
    require_changed_when_different_in_canonical: bool = True,
) -> dict[str, Any]:
    """Validate 1:1 join of `genome_quality.tsv` and `thermoplasmatota.tsv` to `refined_representatives.tsv`.

    Detects and raises `ValueError` on:
    - Duplicate or missing/empty genome/mag identifiers
    - Same row count with wrong/mismatched genome IDs
    - Completeness or contamination discrepancies (> 1e-6)
    - Domain label changes (`d__Bacteria` vs `d__Archaea`)
    - Changes in 237 EQ Thermoplasmatota phylum membership
    - Stale taxonomy label substitution in the generated overlay tables
    """
    validate_primary_key(qual_df, ["genome"], "genome_quality.tsv")
    validate_primary_key(thermo_df, ["genome"], "thermoplasmatota.tsv")
    validate_primary_key(canonical_df, ["mag"], "refined_representatives.tsv")

    req_qual = {"genome", "gtdb_classification", "checkm_completeness", "contamination"}
    req_canon = {"mag", "gtdb_r232_classification", "Completeness", "Contamination"}
    if not req_qual.issubset(qual_df.columns):
        raise ValueError(f"genome_quality.tsv missing required columns: {sorted(req_qual - set(qual_df.columns))}")
    if not req_canon.issubset(canonical_df.columns):
        raise ValueError(
            f"refined_representatives.tsv missing required columns: {sorted(req_canon - set(canonical_df.columns))}"
        )

    qual_ids = set(qual_df["genome"].astype(str))
    canon_ids = set(canonical_df["mag"].astype(str))
    thermo_ids = set(thermo_df["genome"].astype(str))

    if len(qual_df) != len(canonical_df) or qual_ids != canon_ids:
        only_qual = sorted(qual_ids - canon_ids)[:5]
        only_canon = sorted(canon_ids - qual_ids)[:5]
        raise ValueError(
            f"1:1 genome membership mismatch between genome_quality.tsv (n={len(qual_df)}) "
            f"and refined_representatives.tsv (n={len(canonical_df)}): "
            f"only_in_qual={only_qual}, only_in_canonical={only_canon}"
        )

    if not thermo_ids.issubset(qual_ids):
        missing_th = sorted(thermo_ids - qual_ids)[:5]
        raise ValueError(f"thermoplasmatota.tsv contains genomes absent from genome_quality.tsv: {missing_th}")

    canon_indexed = canonical_df.set_index("mag").reindex(qual_df["genome"].astype(str))
    if canon_indexed["gtdb_r232_classification"].isna().any():
        raise ValueError("Null canonical gtdb_r232_classification after 1:1 reindex")

    comp_diff = np.max(
        np.abs(
            pd.to_numeric(qual_df["checkm_completeness"], errors="raise").to_numpy(dtype="float64")
            - pd.to_numeric(canon_indexed["Completeness"], errors="raise").to_numpy(dtype="float64")
        )
    )
    cont_diff = np.max(
        np.abs(
            pd.to_numeric(qual_df["contamination"], errors="raise").to_numpy(dtype="float64")
            - pd.to_numeric(canon_indexed["Contamination"], errors="raise").to_numpy(dtype="float64")
        )
    )
    if float(comp_diff) > 1e-6 or float(cont_diff) > 1e-6:
        raise ValueError(
            f"Completeness/contamination mismatch between genome_quality.tsv and refined_representatives.tsv: "
            f"max_comp_diff={comp_diff}, max_contam_diff={cont_diff}"
        )

    leg_tax_list = qual_df["gtdb_classification"].astype(str).tolist()
    can_tax_list = canon_indexed["gtdb_r232_classification"].astype(str).tolist()

    leg_parsed = [parse_gtdb_lineage(x) for x in leg_tax_list]
    can_parsed = [parse_gtdb_lineage(x) for x in can_tax_list]

    rec_df = pd.DataFrame(
        {
            "genome": qual_df["genome"].astype(str).to_numpy(),
            "checkm_completeness": pd.to_numeric(qual_df["checkm_completeness"], errors="raise").to_numpy(dtype="float64"),
            "canonical_completeness": pd.to_numeric(canon_indexed["Completeness"], errors="raise").to_numpy(dtype="float64"),
            "contamination": pd.to_numeric(qual_df["contamination"], errors="raise").to_numpy(dtype="float64"),
            "canonical_contamination": pd.to_numeric(canon_indexed["Contamination"], errors="raise").to_numpy(dtype="float64"),
            "legacy_gtdb_classification": leg_tax_list,
            "canonical_gtdb_r232_classification": can_tax_list,
            "whole_taxonomy_differs": [l != c for l, c in zip(leg_tax_list, can_tax_list)],
        }
    )

    rank_diff_counts: dict[str, int] = {}
    for rk, rname in RANK_KEYS:
        l_col = [d[rk] for d in leg_parsed]
        c_col = [d[rk] for d in can_parsed]
        diff_col = [l != c for l, c in zip(l_col, c_col)]
        rec_df[f"legacy_{rname}"] = l_col
        rec_df[f"canonical_{rname}"] = c_col
        rec_df[f"{rname}_differs"] = diff_col
        rank_diff_counts[f"n_{rname}_differs"] = int(sum(diff_col))

    if rank_diff_counts["n_domain_differs"] != 0:
        raise ValueError(f"Domain labels differ for {rank_diff_counts['n_domain_differs']} genomes; expected 0")

    # Verify Thermoplasmatota phylum membership is identical between legacy and canonical
    canon_thermo_ids = set(rec_df.loc[rec_df["canonical_phylum"] == "p__Thermoplasmatota", "genome"])
    if canon_thermo_ids != thermo_ids:
        raise ValueError(
            f"Canonical Thermoplasmatota genome set (n={len(canon_thermo_ids)}) differs from thermoplasmatota.tsv (n={len(thermo_ids)})"
        )

    # Build overlay tables with legacy_gtdb_classification retained and gtdb_classification replaced
    overlay_qual_df = qual_df.copy()
    overlay_qual_df["legacy_gtdb_classification"] = overlay_qual_df["gtdb_classification"].astype(str)
    overlay_qual_df["gtdb_classification"] = can_tax_list

    canon_thermo_indexed = canonical_df.set_index("mag").reindex(thermo_df["genome"].astype(str))
    th_comp_diff = np.max(
        np.abs(
            pd.to_numeric(thermo_df["checkm_completeness"], errors="raise").to_numpy(dtype="float64")
            - pd.to_numeric(canon_thermo_indexed["Completeness"], errors="raise").to_numpy(dtype="float64")
        )
    )
    th_cont_diff = np.max(
        np.abs(
            pd.to_numeric(thermo_df["contamination"], errors="raise").to_numpy(dtype="float64")
            - pd.to_numeric(canon_thermo_indexed["Contamination"], errors="raise").to_numpy(dtype="float64")
        )
    )
    if float(th_comp_diff) > 1e-6 or float(th_cont_diff) > 1e-6:
        raise ValueError("Thermoplasmatota completeness/contamination mismatch vs canonical refined_representatives.tsv")

    overlay_thermo_df = thermo_df.copy()
    overlay_thermo_df["legacy_gtdb_classification"] = overlay_thermo_df["gtdb_classification"].astype(str)
    overlay_thermo_df["gtdb_classification"] = canon_thermo_indexed["gtdb_r232_classification"].astype(str).to_numpy()

    verify_overlay_not_stale(overlay_qual_df, canonical_df, require_diff=require_changed_when_different_in_canonical)
    verify_overlay_not_stale(overlay_thermo_df, canonical_df, require_diff=require_changed_when_different_in_canonical)

    # Subset reconciliation for the 237 Thermoplasmatota genomes
    rec_thermo_df = rec_df[rec_df["genome"].isin(thermo_ids)].copy().reset_index(drop=True)
    leg_g_blocks, leg_unclass = normalize_genus_groups(
        rec_thermo_df["genome"].tolist(), rec_thermo_df["legacy_genus"].tolist()
    )
    can_g_blocks, can_unclass = normalize_genus_groups(
        rec_thermo_df["genome"].tolist(), rec_thermo_df["canonical_genus"].tolist()
    )
    rec_thermo_df["legacy_genus_singleton_block"] = leg_g_blocks
    rec_thermo_df["legacy_is_unclassified_genus"] = leg_unclass
    rec_thermo_df["canonical_genus_singleton_block"] = can_g_blocks
    rec_thermo_df["canonical_is_unclassified_genus"] = can_unclass

    # Summarise Thermoplasmatota class/order/genus counts under legacy vs canonical
    can_orders = rec_thermo_df["canonical_order"]
    leg_orders = rec_thermo_df["legacy_order"]
    can_named_orders = sorted({x for x in can_orders if x not in ("", "o__")})
    leg_named_orders = sorted({x for x in leg_orders if x not in ("", "o__")})
    can_named_genera = sorted({x for x in rec_thermo_df["canonical_genus"] if x not in ("", "g__")})
    leg_named_genera = sorted({x for x in rec_thermo_df["legacy_genus"] if x not in ("", "g__")})

    summary = {
        "n_total_genomes": int(len(rec_df)),
        "n_whole_taxonomy_strings_differ": int(rec_df["whole_taxonomy_differs"].sum()),
        **rank_diff_counts,
        "max_completeness_abs_diff": float(comp_diff),
        "max_contamination_abs_diff": float(cont_diff),
        "thermoplasmatota_237": {
            "n_genomes": int(len(rec_thermo_df)),
            "n_whole_taxonomy_strings_differ": int(rec_thermo_df["whole_taxonomy_differs"].sum()),
            "n_class_differs": int(rec_thermo_df["class_differs"].sum()),
            "n_order_differs": int(rec_thermo_df["order_differs"].sum()),
            "n_family_differs": int(rec_thermo_df["family_differs"].sum()),
            "n_genus_differs": int(rec_thermo_df["genus_differs"].sum()),
            "canonical_class_counts": {
                str(k): int(v) for k, v in rec_thermo_df["canonical_class"].value_counts().items()
            },
            "legacy_named_orders_n": len(leg_named_orders),
            "legacy_named_orders": leg_named_orders,
            "legacy_unassigned_order_genomes_n": int(leg_orders.isin(["", "o__"]).sum()),
            "canonical_named_orders_n": len(can_named_orders),
            "canonical_named_orders": can_named_orders,
            "canonical_unassigned_order_genomes_n": int(can_orders.isin(["", "o__"]).sum()),
            "legacy_named_genera_n": len(leg_named_genera),
            "legacy_unassigned_genus_genomes_n": int(rec_thermo_df["legacy_is_unclassified_genus"].sum()),
            "legacy_genus_singleton_blocks_n": int(rec_thermo_df["legacy_genus_singleton_block"].nunique()),
            "canonical_named_genera_n": len(can_named_genera),
            "canonical_unassigned_genus_genomes_n": int(rec_thermo_df["canonical_is_unclassified_genus"].sum()),
            "canonical_genus_singleton_blocks_n": int(rec_thermo_df["canonical_genus_singleton_block"].nunique()),
        },
        "provenance_note": (
            "Reconciles analytical inputs/genome_quality.tsv and inputs/thermoplasmatota.tsv against data-paper "
            "canonical metadata/metagenome/canonical/refined_representatives.tsv (gtdb_r232_classification) by exact "
            "1:1 genome==mag join. Does not speculate on the release or cause of the legacy taxonomy strings."
        ),
    }

    return {
        "reconciled_5438_df": rec_df,
        "reconciled_thermo237_df": rec_thermo_df,
        "overlay_qual_df": overlay_qual_df,
        "overlay_thermo_df": overlay_thermo_df,
        "summary": summary,
    }


def verify_overlay_not_stale(
    overlay_df: pd.DataFrame,
    canonical_df: pd.DataFrame,
    *,
    require_diff: bool = True,
) -> None:
    """Verify that `overlay_df` has canonical `gtdb_r232_classification` in `gtdb_classification` and retains `legacy_gtdb_classification`."""
    if "legacy_gtdb_classification" not in overlay_df.columns or "gtdb_classification" not in overlay_df.columns:
        raise ValueError("Overlay table must contain both 'gtdb_classification' and 'legacy_gtdb_classification'")

    canon_map = dict(
        zip(canonical_df["mag"].astype(str), canonical_df["gtdb_r232_classification"].astype(str))
    )
    expected_canon = [canon_map.get(str(g)) for g in overlay_df["genome"].astype(str)]
    if any(x is None for x in expected_canon):
        raise ValueError("Overlay table contains genome IDs absent from canonical_df")

    actual_gtdb = overlay_df["gtdb_classification"].astype(str).tolist()
    if actual_gtdb != expected_canon:
        raise ValueError("Stale or mismatched taxonomy label substitution detected in overlay gtdb_classification")

    if require_diff:
        leg_gtdb = overlay_df["legacy_gtdb_classification"].astype(str).tolist()
        expected_diffs = sum(1 for l, c in zip(leg_gtdb, expected_canon) if l != c)
        actual_diffs = sum(1 for l, a in zip(leg_gtdb, actual_gtdb) if l != a)
        if expected_diffs > 0 and actual_diffs == 0:
            raise ValueError("Overlay gtdb_classification is still identical to legacy_gtdb_classification")


def build_canonical_overlay_tree(
    orig_arch_root: Path,
    orig_sci_in_root: Path,
    canonical_rep_tsv: Path,
    out_root: Path,
) -> dict[str, Any]:
    """Create the canonical output root and read-only overlay_arch tree without touching frozen source files."""
    require_file(canonical_rep_tsv, "canonical refined_representatives.tsv")
    qual_path = require_file(orig_arch_root / "inputs/genome_quality.tsv", "legacy genome_quality.tsv")
    thermo_path = require_file(orig_arch_root / "inputs/thermoplasmatota.tsv", "legacy thermoplasmatota.tsv")

    qual_df = pd.read_csv(qual_path, sep="\t")
    thermo_df = pd.read_csv(thermo_path, sep="\t")
    canon_df = pd.read_csv(canonical_rep_tsv, sep="\t")

    rec = reconcile_and_validate_canonical_taxonomy(qual_df, thermo_df, canon_df)

    in_dir = out_root / "inputs"
    res_dir = out_root / "results"
    logs_dir = out_root / "logs"
    overlay_arch = out_root / "overlay_arch"
    overlay_arch_in = overlay_arch / "inputs"

    for d in (in_dir, res_dir, logs_dir, overlay_arch_in):
        d.mkdir(parents=True, exist_ok=True)

    # Copy canonical refined_representatives.tsv into out_root/inputs/ and record SHA-256
    staged_canon_path = in_dir / "refined_representatives.tsv"
    if staged_canon_path.resolve() != canonical_rep_tsv.resolve():
        shutil.copy2(canonical_rep_tsv, staged_canon_path)
    canon_sha256 = sha256_file(staged_canon_path)

    # Symlink frozen science inputs into out_root/inputs/
    for item in sorted(orig_sci_in_root.iterdir()):
        dst = in_dir / item.name
        if dst.name == "refined_representatives.tsv":
            continue
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(item.resolve())

    # Symlink top-level directories/files in orig_arch_root (except inputs/) into overlay_arch
    for item in sorted(orig_arch_root.iterdir()):
        if item.name == "inputs":
            continue
        dst = overlay_arch / item.name
        if dst.exists() or dst.is_symlink():
            if dst.is_symlink() or dst.is_file():
                dst.unlink()
        if not dst.exists():
            dst.symlink_to(item.resolve())

    # Symlink non-overlaid files in orig_arch_root/inputs/ into overlay_arch/inputs/
    for item in sorted((orig_arch_root / "inputs").iterdir()):
        if item.name in ("genome_quality.tsv", "thermoplasmatota.tsv"):
            continue
        dst = overlay_arch_in / item.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(item.resolve())

    # Write canonical overlaid genome_quality.tsv and thermoplasmatota.tsv
    ov_qual_path = overlay_arch_in / "genome_quality.tsv"
    ov_thermo_path = overlay_arch_in / "thermoplasmatota.tsv"
    rec["overlay_qual_df"].to_csv(ov_qual_path, sep="\t", index=False)
    rec["overlay_thermo_df"].to_csv(ov_thermo_path, sep="\t", index=False)

    # Write reconciliation tables and summary JSON
    rec_5438_path = res_dir / "taxonomy_reconciliation_all5438.tsv.gz"
    rec_237_path = res_dir / "taxonomy_reconciliation_thermo237.tsv"
    rec_sum_path = res_dir / "taxonomy_reconciliation_summary.json"

    rec["reconciled_5438_df"].to_csv(rec_5438_path, sep="\t", index=False, compression="gzip")
    rec["reconciled_thermo237_df"].to_csv(rec_237_path, sep="\t", index=False)

    summary_payload = {
        **rec["summary"],
        "input_sha256": {
            "refined_representatives.tsv": canon_sha256,
            "legacy_genome_quality.tsv": sha256_file(qual_path),
            "legacy_thermoplasmatota.tsv": sha256_file(thermo_path),
            "overlay_genome_quality.tsv": sha256_file(ov_qual_path),
            "overlay_thermoplasmatota.tsv": sha256_file(ov_thermo_path),
        },
    }
    rec_sum_path.write_text(json.dumps(summary_payload, indent=2) + "\n")
    return summary_payload


def run_canonical_stage1_and_stage1b(
    overlay_arch_root: Path,
    out_root: Path,
    legacy_out_root: Path,
) -> dict[str, Any]:
    """Run existing Stage 1 and Stage 1b algorithms on the canonical taxonomy overlay by overriding module globals."""
    import stage1_audit_and_export as stg1
    import stage1b_refined_fasta_and_corrections as stg1b

    os.environ["EQ_SCIENCE_OUT"] = str(out_root)
    os.environ["EQ_ARCH_ROOT"] = str(overlay_arch_root)
    os.environ["EQ_SCIENCE_IN_ROOT"] = str(out_root / "inputs")

    # Explicitly override module-level path globals before calling main()
    stg1.P_ARCH = overlay_arch_root
    stg1.P_OUT_ROOT = out_root
    stg1.P_IN = out_root / "inputs"
    stg1.P_RES = out_root / "results"

    stg1b.P_ARCH = overlay_arch_root
    stg1b.P_OUT_ROOT = out_root
    stg1b.P_IN = out_root / "inputs"
    stg1b.P_RES = out_root / "results"

    print("[Canonical Wrapper] Running Stage 1 (stage1_audit_and_export.main) on canonical overlay...", flush=True)
    stg1.main()

    print("[Canonical Wrapper] Running Stage 1b (stage1b_refined_fasta_and_corrections.main) on canonical overlay...", flush=True)
    stg1b.main()

    return compare_canonical_vs_legacy_stage1_1b(out_root / "results", legacy_out_root / "results")


def compare_canonical_vs_legacy_stage1_1b(canon_res: Path, legacy_res: Path) -> dict[str, Any]:
    """Verify CDS audit invariance and compare canonical vs legacy lineage recurrence & Task B export models."""
    can_mem = pd.read_csv(
        require_file(canon_res / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz", "canonical member loci"),
        sep="\t",
    ).sort_values(["family_id", "gene"]).reset_index(drop=True)
    leg_mem = pd.read_csv(
        require_file(legacy_res / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz", "legacy member loci"),
        sep="\t",
    ).sort_values(["family_id", "gene"]).reset_index(drop=True)

    # CDS/ORF verification columns must be 100% invariant to taxonomy label updates
    invariant_cols = [
        "family_id",
        "gene",
        "genome",
        "contig",
        "start",
        "end",
        "strand",
        "cds_translation_exact_match",
        "has_admissible_start_codon",
        "has_terminal_stop_codon",
        "verification_status",
        "final_promotion_gate_pass",
        "verified_complete_interior_member",
    ]
    for col in invariant_cols:
        if not can_mem[col].equals(leg_mem[col]):
            raise ValueError(f"CDS audit invariance violated on column '{col}' between canonical and legacy Stage 1b")

    comp_rows: list[dict[str, Any]] = []

    # Compare Route A1 and Route A2 eligible family counts and recurrence
    for route_tag, fname in [
        ("Route_A1_gapfill", "task_a_stage1b_route1_gapfill_family_ranking_verified.tsv"),
        ("Route_A2_exported_dark", "task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv"),
    ]:
        c_df = pd.read_csv(require_file(canon_res / fname, f"canonical {fname}"), sep="\t")
        l_df = pd.read_csv(require_file(legacy_res / fname, f"legacy {fname}"), sep="\t")
        merged = c_df.merge(l_df, on="family_id", suffixes=("_canonical", "_legacy"))
        n_gblock_changed = int(
            (merged["n_verified_independent_genus_blocks_canonical"] != merged["n_verified_independent_genus_blocks_legacy"]).sum()
        )
        n_elig_can = int(c_df["final_candidate_promotion_eligible"].astype(bool).sum())
        n_elig_leg = int(l_df["final_candidate_promotion_eligible"].astype(bool).sum())
        comp_rows.append(
            {
                "comparison_category": route_tag,
                "metric": "n_final_candidate_promotion_eligible_families",
                "legacy_value": float(n_elig_leg),
                "canonical_value": float(n_elig_can),
                "canonical_minus_legacy": float(n_elig_can - n_elig_leg),
                "details": f"n_families_total={len(merged)};n_families_with_changed_verified_genus_blocks={n_gblock_changed}",
            }
        )

    # Compare Task B Stage 1b export models (within-genome, within-named-genus, within-genus-block, order-FE GLM)
    can_s1b = pd.read_csv(
        require_file(
            canon_res / "task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv",
            "canonical stage1b export sensitivity",
        ),
        sep="\t",
    )
    leg_s1b = pd.read_csv(
        require_file(
            legacy_res / "task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv",
            "legacy stage1b export sensitivity",
        ),
        sep="\t",
    )
    m_s1b = can_s1b.merge(
        leg_s1b,
        on=["cohort", "endpoint", "estimator"],
        suffixes=("_canonical", "_legacy"),
    )
    for r in m_s1b.itertuples(index=False):
        comp_rows.append(
            {
                "comparison_category": f"TaskB_Stage1b_{r.cohort}",
                "metric": f"{r.endpoint}::{r.estimator}",
                "legacy_value": float(r.or_estimate_legacy),
                "canonical_value": float(r.or_estimate_canonical),
                "canonical_minus_legacy": float(r.or_estimate_canonical - r.or_estimate_legacy),
                "details": (
                    f"legacy_CI=[{r.or_ci95_low_legacy:.3f},{r.or_ci95_high_legacy:.3f}](units={r.n_independent_lineage_units_legacy});"
                    f"canonical_CI=[{r.or_ci95_low_canonical:.3f},{r.or_ci95_high_canonical:.3f}](units={r.n_independent_lineage_units_canonical})"
                ),
            }
        )

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(canon_res / "canonical_vs_legacy_stage1_1b_comparison.tsv", sep="\t", index=False)
    return {
        "cds_audit_invariance_verified": True,
        "n_member_loci_checked": int(len(can_mem)),
        "n_comparison_rows": int(len(comp_df)),
    }


def build_10_contiguous_spatial_blocks(design_df: pd.DataFrame, n_blocks: int = 10) -> dict[str, str]:
    """Construct the exact 10 contiguous spatial blocks from mean site `transect_km` matching ecology-integration."""
    if "site" not in design_df.columns or "transect_km" not in design_df.columns:
        raise ValueError("design_df must contain 'site' and 'transect_km' columns to build spatial blocks")
    site_order = (
        design_df.assign(_tkm=pd.to_numeric(design_df["transect_km"], errors="raise"))
        .groupby("site")["_tkm"]
        .mean()
        .sort_values()
        .index
    )
    if len(site_order) < n_blocks:
        raise ValueError(f"Cannot split {len(site_order)} sites into {n_blocks} spatial blocks")
    return {
        str(s): f"spatial_block_{i:02d}"
        for i, grp in enumerate(np.array_split(site_order, n_blocks))
        for s in grp
    }


def compute_site_blocked_genus_turnover(
    V_use: pd.DataFrame,
    t_vec: np.ndarray,
    genomes_all: list[str],
    genera_labels: list[str],
    *,
    c_target: str,
    normalize_to_mapped: bool,
    site_to_block: dict[str, str],
    n_boot: int = 1000,
    seed: int = 20261108,
) -> dict[str, Any]:
    """Compute exact symmetric genus turnover decomposition with 10-contiguous-spatial-block bootstrap CIs.

    Retains equal site-weighting and resample multiplicity across the 10 spatial blocks (`site_to_block`),
    and additionally reports individual-site bootstrap CIs as an explicitly labelled sensitivity.
    Verifies exact decomposition sum (`between_shared + within_shared + exclusive == total_diff`) at
    every paired visit and in the site-mean estimate.
    """
    left_v = V_use.xs(c_target, level=2)
    right_v = V_use.xs("Deep", level=2)
    shared_visits = left_v.index.intersection(right_v.index)
    if len(shared_visits) == 0:
        raise ValueError(f"Zero paired visits between {c_target} and Deep")

    site_visit_records: dict[str, list[dict[str, float]]] = {}
    paired_visits_list: list[tuple[str, str]] = []
    for site, trip in shared_visits:
        wa = left_v.loc[(site, trip)].to_numpy(dtype="float64")
        wb = right_v.loc[(site, trip)].to_numpy(dtype="float64")
        dec = decompose_trait_turnover(
            wa, wb, t_vec, genomes_all, genera_labels, normalize_to_mapped=normalize_to_mapped
        )
        sum_err = abs((dec["between_shared"] + dec["within_shared"] + dec["exclusive"]) - dec["total_diff"])
        if sum_err > 1e-9:
            raise ValueError(
                f"Symmetric turnover decomposition identity violated at ({site}, {trip}): err={sum_err}"
            )
        site_visit_records.setdefault(str(site), []).append(dec)
        paired_visits_list.append((str(site), str(trip)))

    sites_ordered = sorted(site_visit_records.keys())
    missing_sites = [s for s in sites_ordered if s not in site_to_block or not str(site_to_block[s]).strip()]
    if missing_sites:
        raise ValueError(f"Paired sites missing from site_to_block mapping: {missing_sites}")
    site_blocks_ordered = [str(site_to_block[s]) for s in sites_ordered]

    comp_names = ("total_diff", "between_shared", "within_shared", "exclusive", "unclassified_mass")
    site_mat = np.zeros((len(sites_ordered), len(comp_names)), dtype="float64")
    for idx_s, s_id in enumerate(sites_ordered):
        r_list = site_visit_records[s_id]
        for idx_c, cname in enumerate(comp_names):
            site_mat[idx_s, idx_c] = 100.0 * float(np.mean([r[cname] for r in r_list]))

    def _mean_across_sites(idx: np.ndarray) -> np.ndarray:
        return np.mean(site_mat[idx], axis=0)

    # Primary bootstrap: resample 10 contiguous spatial blocks (retaining equal site weight & multiplicity)
    boot_spatial10 = block_bootstrap_stat(site_blocks_ordered, _mean_across_sites, n_boot=n_boot, seed=seed)
    # Sensitivity bootstrap: resample individual sites
    boot_site_only = block_bootstrap_stat(sites_ordered, _mean_across_sites, n_boot=n_boot, seed=seed)

    est_vec = boot_spatial10["estimate"]
    mean_sum_err = abs(float((est_vec[1] + est_vec[2] + est_vec[3]) - est_vec[0]))
    if mean_sum_err > 1e-9:
        raise ValueError(f"Site-mean decomposition sum identity violated: err={mean_sum_err}")

    out: dict[str, Any] = {
        "_paired_sites": tuple(sites_ordered),
        "_paired_visits": tuple(sorted(paired_visits_list)),
        "_n_spatial_blocks": int(boot_spatial10["n_blocks"]),
    }
    for idx_c, cname in enumerate(comp_names):
        out[cname] = {
            "mean_pp": float(boot_spatial10["estimate"][idx_c]),
            "se_pp": float(boot_spatial10["se"][idx_c]),
            "ci95_low_pp": float(boot_spatial10["ci95_low"][idx_c]),
            "ci95_high_pp": float(boot_spatial10["ci95_high"][idx_c]),
            "site_only_se_pp": float(boot_site_only["se"][idx_c]),
            "site_only_ci95_low_pp": float(boot_site_only["ci95_low"][idx_c]),
            "site_only_ci95_high_pp": float(boot_site_only["ci95_high"][idx_c]),
            "n_spatial_blocks": int(boot_spatial10["n_blocks"]),
            "n_paired_sites": int(len(sites_ordered)),
            "n_paired_visits": int(len(shared_visits)),
        }
    return out


def validate_3guild_turnover_table(df: pd.DataFrame) -> None:
    """Validate required primary panels, total contrast invariance, decomposition sum, and paired memberships."""
    comp_df = df[df["panel_status"] == "completed"].copy()
    for req_panel in PRIMARY_REQUIRED_3GUILD_PANELS:
        sub_p = comp_df[comp_df["panel"] == req_panel]
        if len(sub_p) != 20:
            raise ValueError(
                f"Required primary panel '{req_panel}' must have 20 completed rows "
                f"(2 contrasts x 2 denominators x 5 components), found {len(sub_p)}"
            )

    tot_rows = comp_df[comp_df["component"] == "total_diff"]
    if tot_rows.empty:
        raise ValueError("No completed total_diff rows in 3-guild turnover table")
    if not (tot_rows["canonical_minus_legacy_pp"].abs() < 1e-9).all():
        max_d = float(tot_rows["canonical_minus_legacy_pp"].abs().max())
        raise ValueError(f"Total contrast invariance violated between canonical and legacy: max_diff={max_d}")
    if not ((tot_rows["canonical_r232_ci95_low_pp"] - tot_rows["legacy_ci95_low_pp"]).abs() < 1e-9).all():
        raise ValueError("Total contrast 95% CI low bound differs between canonical and legacy")
    if not ((tot_rows["canonical_r232_ci95_high_pp"] - tot_rows["legacy_ci95_high_pp"]).abs() < 1e-9).all():
        raise ValueError("Total contrast 95% CI high bound differs between canonical and legacy")

    if not comp_df["paired_membership_identical"].astype(bool).all():
        raise ValueError("Paired site/visit membership mismatch between canonical and legacy turnover")

    for (p, ct, dm), grp in comp_df.groupby(["panel", "contrast", "denominator_mode"], sort=False):
        by_c = grp.set_index("component")
        for col in ("legacy_mean_pp", "canonical_r232_mean_pp"):
            tot_v = float(by_c.loc["total_diff", col])
            sum_parts = (
                float(by_c.loc["between_shared", col])
                + float(by_c.loc["within_shared", col])
                + float(by_c.loc["exclusive", col])
            )
            if abs(sum_parts - tot_v) > 1e-9:
                raise ValueError(
                    f"Decomposition sum mismatch for ({p}, {ct}, {dm}, {col}): sum={sum_parts} vs total={tot_v}"
                )


def run_3guild_canonical_genus_turnover_sensitivity(
    overlay_arch_root: Path,
    orig_arch_root: Path,
    out_root: Path,
    *,
    eco_int_root: Path = DEFAULT_ECO_INT,
    eco_mag_root: Path = DEFAULT_ECO_MAG,
    n_boot: int = 1000,
) -> pd.DataFrame:
    """Recompute exact genus turnover for BetL, AhpCF, and group1l under canonical vs legacy genus labels."""
    res_dir = out_root / "results"
    res_dir.mkdir(parents=True, exist_ok=True)

    design_df = pd.read_csv(require_file(eco_int_root / "results/design.tsv", "design.tsv"), sep="\t")
    validate_primary_key(design_df, ["sample_id"], "design.tsv")
    validate_audited_compartment(design_df, "compartment")
    site_to_block = build_10_contiguous_spatial_blocks(design_df, n_blocks=10)

    M_DIR = eco_mag_root / "inputs/ibex/rebuild_20260729/picrust2_current_20260803/mag_ko"
    mag_index_df = pd.read_csv(require_file(M_DIR / "mag_index.tsv", "mag_index.tsv"), sep="\t").sort_values("mag_index")
    genomes_all = mag_index_df["mag_id"].astype(str).tolist()

    qual_can_df = (
        pd.read_csv(require_file(overlay_arch_root / "inputs/genome_quality.tsv", "canonical genome_quality.tsv"), sep="\t")
        .set_index("genome")
        .reindex(genomes_all)
    )
    qual_leg_df = (
        pd.read_csv(require_file(orig_arch_root / "inputs/genome_quality.tsv", "legacy genome_quality.tsv"), sep="\t")
        .set_index("genome")
        .reindex(genomes_all)
    )

    genera_can = [
        next((x for x in str(tax).split(";") if x.startswith("g__")), "g__")
        for tax in qual_can_df["gtdb_classification"]
    ]
    genera_leg = [
        next((x for x in str(tax).split(";") if x.startswith("g__")), "g__")
        for tax in qual_leg_df["gtdb_classification"]
    ]

    AA_raw = sparse.load_npz(
        require_file(M_DIR / "sample_by_mag_coverm_relative_abundance.npz", "coverm npz")
    ).toarray()[design_df["sample_index"].to_numpy()]
    row_mapped_fractions = AA_raw.sum(axis=1)

    def _mean_visit_weights(W_mat: np.ndarray) -> pd.DataFrame:
        df_w = pd.DataFrame(W_mat)
        df_w[["site", "trip", "compartment"]] = design_df[["site", "trip", "compartment"]].to_numpy()
        return df_w.groupby(["site", "trip", "compartment"]).mean()

    V_allread = _mean_visit_weights(AA_raw)
    V_mapped = _mean_visit_weights(AA_raw / row_mapped_fractions[:, None])

    G_traits = (
        pd.read_csv(require_file(eco_int_root / "results/genome_traits.tsv", "genome_traits.tsv"), sep="\t")
        .set_index("genome")
        .reindex(genomes_all)
    )

    # Hard check on required primary 3-guild panels: missing or non-finite MUST fail immediately
    for req_panel in PRIMARY_REQUIRED_3GUILD_PANELS:
        if req_panel not in G_traits.columns:
            raise ValueError(f"Required primary 3-guild panel '{req_panel}' missing from genome_traits.tsv")
        req_vec = pd.to_numeric(G_traits[req_panel], errors="coerce").to_numpy(dtype="float64")
        if not np.isfinite(req_vec).all():
            raise ValueError(f"Required primary 3-guild panel '{req_panel}' contains non-finite values")

    rows: list[dict[str, Any]] = []
    for panel in THREE_GUILD_TURNOVER_PANELS:
        is_primary = panel in PRIMARY_REQUIRED_3GUILD_PANELS
        if panel not in G_traits.columns or not np.isfinite(pd.to_numeric(G_traits[panel], errors="coerce").to_numpy(dtype="float64")).all():
            if is_primary:
                raise ValueError(f"Primary panel '{panel}' unavailable or non-finite")
            rows.append(
                {
                    "panel": panel,
                    "is_primary_required_3guild": False,
                    "panel_status": "optional_panel_unavailable_or_nonfinite",
                    "contrast": "NA",
                    "denominator_mode": "NA",
                    "normalize_to_mapped": False,
                    "component": "NA",
                    "bootstrap_block_scheme": "ten_contiguous_spatial_blocks_transect_km",
                    "n_spatial_blocks": 0,
                    "n_paired_sites": 0,
                    "n_paired_visits": 0,
                    "paired_membership_identical": True,
                    "legacy_mean_pp": np.nan,
                    "legacy_se_pp": np.nan,
                    "legacy_ci95_low_pp": np.nan,
                    "legacy_ci95_high_pp": np.nan,
                    "canonical_r232_mean_pp": np.nan,
                    "canonical_r232_se_pp": np.nan,
                    "canonical_r232_ci95_low_pp": np.nan,
                    "canonical_r232_ci95_high_pp": np.nan,
                    "canonical_r232_site_only_ci95_low_pp": np.nan,
                    "canonical_r232_site_only_ci95_high_pp": np.nan,
                    "canonical_minus_legacy_pp": np.nan,
                    "unclassified_handling": "optional_panel_unavailable",
                }
            )
            continue

        t_vec = G_traits[panel].to_numpy(dtype="float64")
        for c_target in ("Surface", "RootAdjacent"):
            contrast_name = f"{c_target}-minus-Deep"
            for denom_mode, V_use, norm_flag in [
                ("mapped_community_normalized_weights", V_mapped, True),
                ("all_reads_unnormalized_weights", V_allread, False),
            ]:
                can_dec = compute_site_blocked_genus_turnover(
                    V_use,
                    t_vec,
                    genomes_all,
                    genera_can,
                    c_target=c_target,
                    normalize_to_mapped=norm_flag,
                    site_to_block=site_to_block,
                    n_boot=n_boot,
                )
                leg_dec = compute_site_blocked_genus_turnover(
                    V_use,
                    t_vec,
                    genomes_all,
                    genera_leg,
                    c_target=c_target,
                    normalize_to_mapped=norm_flag,
                    site_to_block=site_to_block,
                    n_boot=n_boot,
                )
                same_pairs = (
                    can_dec["_paired_sites"] == leg_dec["_paired_sites"]
                    and can_dec["_paired_visits"] == leg_dec["_paired_visits"]
                )
                for comp_name in ("total_diff", "between_shared", "within_shared", "exclusive", "unclassified_mass"):
                    c_info = can_dec[comp_name]
                    l_info = leg_dec[comp_name]
                    rows.append(
                        {
                            "panel": panel,
                            "is_primary_required_3guild": is_primary,
                            "panel_status": "completed",
                            "contrast": contrast_name,
                            "denominator_mode": denom_mode,
                            "normalize_to_mapped": norm_flag,
                            "component": comp_name,
                            "bootstrap_block_scheme": "ten_contiguous_spatial_blocks_transect_km",
                            "n_spatial_blocks": c_info["n_spatial_blocks"],
                            "n_paired_sites": c_info["n_paired_sites"],
                            "n_paired_visits": c_info["n_paired_visits"],
                            "paired_membership_identical": bool(same_pairs),
                            "legacy_mean_pp": l_info["mean_pp"],
                            "legacy_se_pp": l_info["se_pp"],
                            "legacy_ci95_low_pp": l_info["ci95_low_pp"],
                            "legacy_ci95_high_pp": l_info["ci95_high_pp"],
                            "canonical_r232_mean_pp": c_info["mean_pp"],
                            "canonical_r232_se_pp": c_info["se_pp"],
                            "canonical_r232_ci95_low_pp": c_info["ci95_low_pp"],
                            "canonical_r232_ci95_high_pp": c_info["ci95_high_pp"],
                            "canonical_r232_site_only_ci95_low_pp": c_info["site_only_ci95_low_pp"],
                            "canonical_r232_site_only_ci95_high_pp": c_info["site_only_ci95_high_pp"],
                            "canonical_minus_legacy_pp": float(c_info["mean_pp"] - l_info["mean_pp"]),
                            "unclassified_handling": "unclassified_g___genomes_kept_as_separate_singletons_never_pooled",
                        }
                    )

    out_df = pd.DataFrame(rows)
    validate_3guild_turnover_table(out_df)
    out_df.to_csv(res_dir / "task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv", sep="\t", index=False)

    # Explicitly distinguish completed canonical 3-guild genus turnover from pending secondary gRodon-order regression
    root_growth_status = {
        "canonical_3guild_genus_turnover_status": "COMPLETED_EXACT_R232_DECOMPOSITION_10_SPATIAL_BLOCKS",
        "primary_3guild_panels_verified": PRIMARY_REQUIRED_3GUILD_PANELS,
        "secondary_grodon_order_genus_regression_claim": "root_growth_order_or_genus_association_gRodon",
        "secondary_grodon_order_genus_regression_status": (
            "PENDING_CANONICAL_TAXONOMY_RERUN_VISIBLE_PLACEHOLDER_APPROVED_DO_NOT_ASSERT_LEGACY_AS_CURRENT"
        ),
        "note": (
            "Per canonical-taxonomy-sensitivity-design.txt Item 4 and canonical-runtime-final-review.txt Item 3, "
            "canonical 3-guild genus turnover decomposition (BetL: HMM:glycine betaine transporter, AhpCF: HMM:ahpCF, "
            "Group 1l: typed:[NiFe]_Group_1l) is recomputed directly under canonical gtdb_r232_classification using "
            "the 10 contiguous spatial blocks (transect_km). In contrast, the secondary gRodon/root-growth order/genus "
            "regression remains pending a canonical rerun and uses the supervisor-approved TAXONOMY-SNAPSHOT-SENSITIVITY "
            "placeholder rather than asserting legacy taxonomy labels as current."
        ),
    }
    (res_dir / "canonical_root_growth_taxonomy_status.json").write_text(json.dumps(root_growth_status, indent=2) + "\n")
    return out_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Canonical GTDB r232 taxonomy overlay and sensitivity runner.")
    parser.add_argument(
        "--mode",
        choices=[
            "all_stage1_1b_4",
            "overlay_only",
            "stage1_1b_only",
            "turnover_3guilds_only",
            "stage3_canonical",
            "stage3_summary_only",
        ],
        default="all_stage1_1b_4",
    )
    parser.add_argument("--orig-arch-root", type=Path, default=DEFAULT_ORIG_ARCH)
    parser.add_argument("--orig-sci-out", type=Path, default=DEFAULT_ORIG_SCI)
    parser.add_argument("--canon-out", type=Path, default=DEFAULT_CANON_OUT)
    parser.add_argument(
        "--canonical-rep-tsv",
        type=Path,
        default=DEFAULT_CANON_OUT / "inputs/refined_representatives.tsv",
    )
    args = parser.parse_args()

    overlay_arch = args.canon_out / "overlay_arch"

    if args.mode in ("all_stage1_1b_4", "overlay_only"):
        print("[Canonical Wrapper] Building canonical taxonomy overlay and 1:1 reconciliation tables...", flush=True)
        rec_sum = build_canonical_overlay_tree(
            args.orig_arch_root,
            args.orig_sci_out / "inputs",
            args.canonical_rep_tsv,
            args.canon_out,
        )
        print("[Canonical Wrapper] Reconciliation summary:", json.dumps(rec_sum, indent=2), flush=True)

    if args.mode in ("all_stage1_1b_4", "stage1_1b_only"):
        comp_sum = run_canonical_stage1_and_stage1b(overlay_arch, args.canon_out, args.orig_sci_out)
        print("[Canonical Wrapper] Stage 1/1b canonical comparison:", json.dumps(comp_sum, indent=2), flush=True)

    if args.mode in ("all_stage1_1b_4", "turnover_3guilds_only"):
        print("[Canonical Wrapper] Recomputing 3-guild genus turnover sensitivity (BetL, AhpCF, Group 1l)...", flush=True)
        turn_df = run_3guild_canonical_genus_turnover_sensitivity(
            overlay_arch, args.orig_arch_root, args.canon_out
        )
        print(f"[Canonical Wrapper] Wrote {len(turn_df)} rows to task_d_canonical_vs_legacy_genus_turnover_3guilds.tsv", flush=True)

    if args.mode in ("stage3_canonical", "stage3_summary_only"):
        import stage3_candidate_functional as stg3

        os.environ["EQ_SCIENCE_OUT"] = str(args.canon_out)
        os.environ["EQ_ARCH_ROOT"] = str(overlay_arch)
        os.environ["EQ_SCIENCE_IN_ROOT"] = str(args.canon_out / "inputs")
        stg3.P_ARCH = overlay_arch
        stg3.P_OUT_ROOT = args.canon_out
        stg3.P_IN = args.canon_out / "inputs"
        stg3.P_RES = args.canon_out / "results"
        stg3.P_DOSSIERS = stg3.P_RES / "stage3_family_dossiers"
        if args.mode == "stage3_summary_only":
            print("[Canonical Wrapper] Refreshing Stage 3 metabolic and GSPA summary tables on canonical overlay...", flush=True)
            stg3.refresh_stage3_metabolic_and_gspa_summaries_only()
        else:
            print("[Canonical Wrapper] Running Stage 3 on canonical Stage 1b outputs...", flush=True)
            stg3.main()


if __name__ == "__main__":
    main()
