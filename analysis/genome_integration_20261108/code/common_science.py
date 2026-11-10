#!/usr/bin/env python3
"""Core validated scientific utilities for Empty Quarter 2026-11-08 compute tasks.

Fixes incorporated from supervisor code review (2026-11-08 19:07):
1. Genome quality classification (`classify_genome_quality`):
   - Separates completeness/contamination-only tiers (`comp_contam_tier`, never labelled MIMAG)
     from strict MIMAG tiers (`mimag_strict_tier`) which require explicit RNA criteria
     (23S/16S/5S rRNAs + >=18 tRNAs) for `MIMAG_HQ`. High-completeness/low-contamination
     genomes lacking required RNAs are `MIMAG_MQ_strict`; missing RNA evidence is
     explicitly `unknown_full_HQ_missing_rna`.
   - Uses unambiguous `sensitivity_comp70_contam5_pass` (70/5 sensitivity filter) and
     `catalogue_comp50_contam10_pass` / `strict_comp50_contam_lt10_pass`.
   - Uses exact numeric equality `cont == 10.0` (no `np.isclose` default rtol=1e-5) for the
     boundary genome at 10.0% contamination.
2. Schema- and key-validated merge (`validated_merge`):
   - Explicitly checks key coverage (`left_keys - right_keys`) before merging whenever
     `require_all_left=True`, so `how='left'` and `validate='one_to_many'` also fail
     explicitly on unmatched keys, and validates non-null/non-empty join keys in both
     tables across all cardinality modes.
3. Ecological trait turnover decomposition (`decompose_trait_turnover`):
   - Uses the exact symmetric between-shared-genus (`between_shared`), within-shared-genus
     (`within_shared`), and compartment-exclusive (`exclusive`) decomposition matching
     `proteome-integration-20261106/code/ecology.py`. Unclassified genera are mapped to
     unique singleton groups so their within-group variance is identically zero.
4. ORF & candidate provenance integrity (`audit_orf_integrity`):
   - Missing internal-stop, partial, length, or contig-edge evidence is recorded as
     `'unknown'` (and `pd.NA`), never treated as an observed absence of stops or edge effects.
   - `candidate_promotion_eligible` requires known complete ORF, known full sequence matching
     coordinates, validated coordinates, retained refined contig, and known contig-edge status.
   - Validates vector lengths, finite values, non-negative abundances, and non-empty/non-null
     resampling blocks explicitly.
5. Candidate examples vs positive controls:
   - GrpE/DnaJ/HSP20-adjacent dark families (`F41Dr3::F41Dr3_256695_46`,
     `49Dr1::49Dr1_177214_2`, `24Dr1::24Dr1_1575644_3`) are labelled strictly as
     `prior_heat_stress_candidate_example` (not positive controls). Positive controls
     must be established annotated proteins (`annotated_positive_control`).
6. Frequentist TOST (`evaluate_frequentist_tost`) vs Bayesian ROPE (`evaluate_bayesian_rope`):
   - Never conflates Bayesian INLA credible intervals with frequentist TOST confidence
     intervals. Requires an explicitly justified margin/ROPE with a non-empty justification
     string (`margin_justified=True`, `justification=...`); defaults to no equivalence claim
     (`no_equivalence_claim_unjustified_margin` / `no_equivalence_claim_unjustified_rope`).
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats

ALLOWED_COMPARTMENTS = frozenset({"Deep", "Surface", "RootAdjacent", "Bulk"})
THERMO_PHYLUM_LABELS = frozenset({"Thermoplasmatota", "p__Thermoplasmatota"})
PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES = frozenset(
    {
        "F41Dr3::F41Dr3_256695_46",
        "49Dr1::49Dr1_177214_2",
        "24Dr1::24Dr1_1575644_3",
    }
)


def require_file(path: str | Path, description: str = "input file") -> Path:
    """Ensure a required input file exists and is non-empty, else fail explicitly."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing required {description}: {p}")
    if not p.is_file():
        raise ValueError(f"Expected regular file for {description}, got: {p}")
    if p.stat().st_size == 0:
        raise ValueError(f"Required {description} is empty (0 bytes): {p}")
    return p


def sha256_file(path: str | Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_key_columns_non_null(
    df: pd.DataFrame, key_cols: Sequence[str], table_name: str = "table"
) -> None:
    """Ensure key columns exist and contain no null, empty, or whitespace-padded IDs."""
    missing_cols = [c for c in key_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"{table_name}: missing required key columns {missing_cols}")
    if df.empty:
        raise ValueError(f"{table_name}: dataframe is empty")
    for col in key_cols:
        s = df[col]
        if s.isna().any():
            raise ValueError(f"{table_name}: null values found in key column '{col}'")
        if s.dtype == object or pd.api.types.is_string_dtype(s):
            raw = s.astype(str)
            stripped = raw.str.strip()
            if (stripped == "").any() or stripped.isin({"nan", "None", "<NA>"}).any():
                raise ValueError(f"{table_name}: empty or invalid string IDs in key column '{col}'")
            if (stripped != raw).any():
                raise ValueError(f"{table_name}: leading/trailing whitespace in key column '{col}'")


def validate_primary_key(
    df: pd.DataFrame, key_cols: Sequence[str], table_name: str = "table"
) -> None:
    """Fail explicitly on missing columns, null/empty IDs, or duplicate keys."""
    validate_key_columns_non_null(df, key_cols, table_name=table_name)
    dup_mask = df.duplicated(subset=list(key_cols), keep=False)
    if dup_mask.any():
        dup_examples = df.loc[dup_mask, list(key_cols)].head(5).to_dict(orient="records")
        raise ValueError(
            f"{table_name}: duplicate primary keys on {list(key_cols)} "
            f"(n_dup_rows={int(dup_mask.sum())}, examples={dup_examples})"
        )


def validated_merge(
    left: pd.DataFrame,
    right: pd.DataFrame,
    on: Sequence[str] | str,
    *,
    how: str = "inner",
    validate: str = "one_to_one",
    require_all_left: bool = False,
    left_name: str = "left",
    right_name: str = "right",
) -> pd.DataFrame:
    """Perform a schema-, key-validity-, and coverage-checked merge.

    Checks non-null/non-empty keys in both tables for all merge modes, enforces uniqueness
    for the declared cardinality, and explicitly verifies that every left key exists in
    `right` before merging when `require_all_left=True` (catching unmatched keys even for
    `how='left'` or `validate='one_to_many'`).
    """
    key_cols = [on] if isinstance(on, str) else list(on)
    validate_key_columns_non_null(left, key_cols, table_name=left_name)
    validate_key_columns_non_null(right, key_cols, table_name=right_name)

    if validate in ("one_to_one", "one_to_many"):
        validate_primary_key(left, key_cols, table_name=left_name)
    if validate in ("one_to_one", "many_to_one"):
        validate_primary_key(right, key_cols, table_name=right_name)

    if require_all_left:
        left_keys = set(map(tuple, left[key_cols].itertuples(index=False, name=None)))
        right_keys = set(map(tuple, right[key_cols].itertuples(index=False, name=None)))
        missing_keys = sorted(left_keys - right_keys)
        if missing_keys:
            raise ValueError(
                f"Join {left_name} -> {right_name} on {key_cols} has {len(missing_keys)} "
                f"unmatched left keys (examples: {missing_keys[:5]})"
            )

    return pd.merge(left, right, on=key_cols, how=how, validate=validate)


def validate_audited_compartment(df: pd.DataFrame, col: str = "compartment") -> None:
    """Ensure compartment labels come strictly from the audited vocabulary."""
    if col not in df.columns:
        raise ValueError(f"Missing audited compartment column '{col}'")
    if df[col].isna().any():
        raise ValueError(f"Null compartment labels present in '{col}'")
    vals = set(df[col].astype(str).unique())
    invalid = vals - ALLOWED_COMPARTMENTS
    if invalid:
        raise ValueError(
            f"Invalid or un-audited compartment labels {sorted(invalid)}; "
            f"allowed={sorted(ALLOWED_COMPARTMENTS)}. Never infer from sample IDs or transfer README."
        )


def validate_phylum_scope(
    df: pd.DataFrame,
    phylum_col: str = "phylum",
    allowed: Iterable[str] = THERMO_PHYLUM_LABELS,
) -> None:
    """Fail explicitly if any row leaks outside the target phylum scope."""
    if phylum_col not in df.columns:
        raise ValueError(f"Missing phylum column '{phylum_col}'")
    allowed_set = set(allowed)
    observed = set(df[phylum_col].dropna().astype(str).unique())
    leaked = observed - allowed_set
    if leaked or df[phylum_col].isna().any():
        raise ValueError(
            f"Phylum scope leakage detected in '{phylum_col}': "
            f"leaked={sorted(leaked)}, n_null={int(df[phylum_col].isna().sum())}"
        )


def classify_genome_quality(
    completeness: pd.Series | np.ndarray | Sequence[float],
    contamination: pd.Series | np.ndarray | Sequence[float],
    rna_complete: pd.Series | np.ndarray | Sequence[bool | None] | None = None,
) -> pd.DataFrame:
    """Classify genome quality separating completeness/contamination tiers from strict MIMAG tiers.

    - `sensitivity_comp70_contam5_pass`: completeness >= 70.0 and contamination <= 5.0
    - `strict_comp50_contam_lt10_pass`: completeness >= 50.0 and contamination < 10.0
    - `catalogue_comp50_contam10_pass`: completeness >= 50.0 and contamination <= 10.0
    - `comp_contam_tier` (completeness/contamination-only, never labelled MIMAG):
        * 'comp_gt90_contam_lt5': completeness > 90.0 and contamination < 5.0
        * 'comp_ge50_contam_lt10': not >90/<5, completeness >= 50.0 and contamination < 10.0
        * 'boundary_10pct_contam': completeness >= 50.0 and exact equality contamination == 10.0
        * 'below_50_10': otherwise
    - `mimag_strict_tier` (requires explicit `rna_complete` evidence for MIMAG HQ: 5S+16S+23S rRNAs and >=18 AA tRNAs):
        * For genomes with completeness > 90.0 and contamination < 5.0:
            - rna_complete is True  -> 'MIMAG_HQ'
            - rna_complete is False -> 'MIMAG_MQ_strict' (lacks RNA requirement for HQ)
            - rna_complete is None/NA -> 'unknown_full_HQ_missing_rna'
        * Otherwise:
            - completeness >= 50.0 and contamination < 10.0 -> 'MIMAG_MQ_strict'
            - completeness >= 50.0 and contamination == 10.0 -> 'boundary_10pct_contam'
            - else -> 'MIMAG_LQ_or_fail'
    """
    comp = pd.Series(completeness, dtype="float64").reset_index(drop=True)
    cont = pd.Series(contamination, dtype="float64").reset_index(drop=True)
    if len(comp) != len(cont):
        raise ValueError(f"Length mismatch: completeness ({len(comp)}) vs contamination ({len(cont)})")
    if not np.isfinite(comp.to_numpy()).all() or not np.isfinite(cont.to_numpy()).all():
        raise ValueError("Completeness and contamination must be finite numbers without NaN/Inf")
    if (comp < 0.0).any() or (comp > 100.0).any() or (cont < 0.0).any() or (cont > 100.0).any():
        raise ValueError("Completeness/contamination out of [0, 100] range")

    if rna_complete is None:
        rna_ser = pd.Series([pd.NA] * len(comp), dtype="boolean")
    else:
        if isinstance(rna_complete, (bool, np.bool_)):
            raise ValueError("rna_complete must be a sequence matching completeness length or None")
        rna_raw = list(rna_complete)
        if len(rna_raw) != len(comp):
            raise ValueError(f"Length mismatch: rna_complete ({len(rna_raw)}) vs completeness ({len(comp)})")
        rna_ser = pd.Series(rna_raw, dtype="boolean")

    sens_70_5 = (comp >= 70.0) & (cont <= 5.0)
    strict_50_lt10 = (comp >= 50.0) & (cont < 10.0)
    cat_50_le10 = (comp >= 50.0) & (cont <= 10.0)

    gt90_lt5 = (comp > 90.0) & (cont < 5.0)
    ge50_lt10 = (~gt90_lt5) & (comp >= 50.0) & (cont < 10.0)
    # Exact numeric equality on recorded contamination values (no rtol fuzzy matching)
    exact_10 = (~gt90_lt5) & (comp >= 50.0) & (cont == 10.0)

    comp_contam_tier = np.where(
        gt90_lt5,
        "comp_gt90_contam_lt5",
        np.where(
            ge50_lt10,
            "comp_ge50_contam_lt10",
            np.where(exact_10, "boundary_10pct_contam", "below_50_10"),
        ),
    )

    mimag_tier: list[str] = []
    for i in range(len(comp)):
        if bool(gt90_lt5.iloc[i]):
            r_val = rna_ser.iloc[i]
            if pd.isna(r_val):
                mimag_tier.append("unknown_full_HQ_missing_rna")
            elif bool(r_val):
                mimag_tier.append("MIMAG_HQ")
            else:
                mimag_tier.append("MIMAG_MQ_strict")
        elif bool(ge50_lt10.iloc[i]):
            mimag_tier.append("MIMAG_MQ_strict")
        elif bool(exact_10.iloc[i]):
            mimag_tier.append("boundary_10pct_contam")
        else:
            mimag_tier.append("MIMAG_LQ_or_fail")

    return pd.DataFrame(
        {
            "completeness": comp,
            "contamination": cont,
            "rna_complete": rna_ser,
            "sensitivity_comp70_contam5_pass": sens_70_5,
            "strict_comp50_contam_lt10_pass": strict_50_lt10,
            "catalogue_comp50_contam10_pass": cat_50_le10,
            "comp_contam_tier": comp_contam_tier,
            "mimag_strict_tier": mimag_tier,
        }
    )


def audit_orf_integrity(
    df: pd.DataFrame,
    *,
    partial_col: str = "partial",
    internal_stop_col: str = "internal_stop",
    length_col: str = "length",
    contig_rank_col: str | None = None,
    contig_n_genes_col: str | None = None,
    contig_in_cat_col: str | None = None,
    has_full_seq_col: str | None = None,
    coords_valid_col: str | None = None,
    edge_gene_buffer: int = 2,
) -> pd.DataFrame:
    """Audit ORF completeness, internal stops, contig-edge proximity, and candidate promotion eligibility.

    Missing internal-stop, partial, length, or contig-edge evidence is marked `'unknown'` / `pd.NA`,
    never treated as an observed absence of stops or edge effects.
    """
    for col in (partial_col, length_col):
        if col not in df.columns:
            raise ValueError(f"Missing required ORF column '{col}'")
    out = df.copy()

    # 1. Partial flag audit
    p_raw = out[partial_col]
    p_na = p_raw.isna() | p_raw.astype(str).str.strip().isin({"", "nan", "None", "<NA>", "missing", "unknown"})
    p_str = p_raw.astype(str).str.strip()
    p_complete = (~p_na) & p_str.isin({"00", "0", "False", "false", "0.0"})
    p_incomplete = (~p_na) & (~p_complete)

    # 2. Internal stop audit (missing column or NaN -> unknown, NOT 0 stops)
    if internal_stop_col in out.columns:
        stop_num = pd.to_numeric(out[internal_stop_col], errors="coerce")
        stop_na = stop_num.isna() | (stop_num < 0)
        stop_zero = (~stop_na) & (stop_num == 0)
        stop_pos = (~stop_na) & (stop_num > 0)
    else:
        stop_na = pd.Series(True, index=out.index)
        stop_zero = pd.Series(False, index=out.index)
        stop_pos = pd.Series(False, index=out.index)

    # 3. Length audit
    len_num = pd.to_numeric(out[length_col], errors="coerce")
    len_valid = len_num.notna() & np.isfinite(len_num) & (len_num > 0)

    orf_status: list[str] = []
    orf_complete_nullable: list[bool | object] = []
    stop_nullable: list[bool | object] = []

    for i in range(len(out)):
        if bool(stop_na.iloc[i]):
            stop_nullable.append(pd.NA)
        else:
            stop_nullable.append(bool(stop_pos.iloc[i]))

        if not bool(len_valid.iloc[i]) or bool(p_na.iloc[i]) or bool(stop_na.iloc[i]):
            if bool(p_incomplete.iloc[i]) or bool(stop_pos.iloc[i]):
                orf_status.append("incomplete_or_internal_stop")
                orf_complete_nullable.append(False)
            else:
                orf_status.append("unknown")
                orf_complete_nullable.append(pd.NA)
        elif bool(p_complete.iloc[i]) and bool(stop_zero.iloc[i]):
            orf_status.append("complete")
            orf_complete_nullable.append(True)
        else:
            orf_status.append("incomplete_or_internal_stop")
            orf_complete_nullable.append(False)

    out["orf_status"] = orf_status
    out["orf_is_complete"] = pd.Series(orf_complete_nullable, index=out.index, dtype="boolean")
    out["orf_has_internal_stop"] = pd.Series(stop_nullable, index=out.index, dtype="boolean")

    # 4. Contig-edge audit (missing rank/n_genes -> unknown, NOT interior)
    if contig_rank_col and contig_n_genes_col and contig_rank_col in out.columns and contig_n_genes_col in out.columns:
        rank = pd.to_numeric(out[contig_rank_col], errors="coerce")
        n_genes = pd.to_numeric(out[contig_n_genes_col], errors="coerce")
        valid_pos = rank.notna() & n_genes.notna() & (rank >= 1) & (n_genes >= rank)
        dist_start = rank - 1
        dist_end = n_genes - rank
        min_dist = np.where(valid_pos, np.minimum(dist_start, dist_end), np.nan)
        out["dist_from_contig_end_genes"] = min_dist

        edge_status: list[str] = []
        edge_nullable: list[bool | object] = []
        for i in range(len(out)):
            if not bool(valid_pos.iloc[i]):
                edge_status.append("unknown")
                edge_nullable.append(pd.NA)
            else:
                is_edge = bool((min_dist[i] < edge_gene_buffer) or (n_genes.iloc[i] <= (2 * edge_gene_buffer)))
                edge_status.append("edge" if is_edge else "interior")
                edge_nullable.append(is_edge)
        out["contig_edge_status"] = edge_status
        out["is_contig_edge"] = pd.Series(edge_nullable, index=out.index, dtype="boolean")
    else:
        out["dist_from_contig_end_genes"] = np.nan
        out["contig_edge_status"] = "unknown"
        out["is_contig_edge"] = pd.Series([pd.NA] * len(out), index=out.index, dtype="boolean")

    # 5. Candidate promotion eligibility (strictly requires known full sequence, validated coords, retained contig, complete ORF, and known edge status)
    if contig_in_cat_col and contig_in_cat_col in out.columns:
        contig_ok = out[contig_in_cat_col].astype(str).str.strip().isin({"True", "true", "1", "1.0"})
    else:
        contig_ok = pd.Series(False, index=out.index)

    if has_full_seq_col and has_full_seq_col in out.columns:
        seq_ok = out[has_full_seq_col].astype(str).str.strip().isin({"True", "true", "1", "1.0"})
    else:
        seq_ok = pd.Series(False, index=out.index)

    if coords_valid_col and coords_valid_col in out.columns:
        coords_ok = out[coords_valid_col].astype(str).str.strip().isin({"True", "true", "1", "1.0"})
    else:
        coords_ok = pd.Series(False, index=out.index)

    out["candidate_promotion_eligible"] = (
        (out["orf_status"] == "complete")
        & (out["contig_edge_status"].isin(["interior", "edge"]))
        & contig_ok
        & seq_ok
        & coords_ok
    )
    return out


def normalize_genus_groups(genomes: Sequence[str], genera: Sequence[str | None]) -> tuple[list[str], np.ndarray]:
    """Map unclassified/empty genus labels to unique singleton group IDs so W=0 identically.

    Returns (group_labels, is_unclassified_mask).
    """
    if len(genomes) != len(genera):
        raise ValueError("genomes and genera must have identical length")
    groups: list[str] = []
    unclassified: list[bool] = []
    for g, gen in zip(genomes, genera):
        g_str = "" if g is None or pd.isna(g) else str(g).strip()
        if not g_str:
            raise ValueError("Empty genome ID in normalize_genus_groups")
        s = "" if gen is None or pd.isna(gen) else str(gen).strip()
        if s in ("", "g__", "p__", "unclassified", "Unclassified", "nan", "None", "s__"):
            groups.append(f"unclassified_singleton::{g_str}")
            unclassified.append(True)
        else:
            groups.append(s)
            unclassified.append(False)
    return groups, np.asarray(unclassified, dtype=bool)


def decompose_trait_turnover(
    w_target: np.ndarray | Sequence[float],
    w_ref: np.ndarray | Sequence[float],
    trait: np.ndarray | Sequence[float],
    genomes: Sequence[str],
    genera: Sequence[str | None],
    *,
    normalize_to_mapped: bool = False,
) -> dict[str, float | bool]:
    """Exact symmetric 3-component decomposition of community-weighted trait shift (target - ref).

    Matches `proteome-integration-20261106/code/ecology.py`:
    Preserves the supplied weights by default (`normalize_to_mapped=False`) so all-read
    denominators retain differences in mapped fraction. Set `normalize_to_mapped=True`
    explicitly only when decomposing within the mapped-community fraction.

    For each taxonomic group k (with unclassified genomes assigned unique singleton groups):
      - If present in BOTH target and reference (m_a > 0 and m_b > 0):
          between_shared += (m_a - m_b) * (t_a + t_b) / 2
          within_shared  += (t_a - t_b) * (m_a + m_b) / 2
      - Else (exclusive to target or exclusive to reference):
          exclusive      += a_k @ t_k - b_k @ t_k
    Total difference (a - b) @ t == between_shared + within_shared + exclusive.
    """
    w_a = np.asarray(w_target, dtype="float64")
    w_b = np.asarray(w_ref, dtype="float64")
    t = np.asarray(trait, dtype="float64")
    if not (w_a.ndim == 1 and w_b.ndim == 1 and t.ndim == 1):
        raise ValueError("w_target, w_ref, and trait must be 1D arrays")
    if not (len(w_a) == len(w_b) == len(t) == len(genomes) == len(genera)):
        raise ValueError("Dimension mismatch in decompose_trait_turnover")
    if not (np.isfinite(w_a).all() and np.isfinite(w_b).all() and np.isfinite(t).all()):
        raise ValueError("w_target, w_ref, and trait must contain only finite values (no NaN/Inf)")
    if (w_a < 0.0).any() or (w_b < 0.0).any():
        raise ValueError("Abundance weights cannot be negative")
    sum_a = float(w_a.sum())
    sum_b = float(w_b.sum())
    if sum_a <= 0.0 or sum_b <= 0.0:
        raise ValueError("Target and reference abundance vectors must have positive sum")
    if normalize_to_mapped:
        a = w_a / sum_a
        b = w_b / sum_b
    else:
        a = w_a.copy()
        b = w_b.copy()

    groups, unclassified_mask = normalize_genus_groups(genomes, genera)
    unique_groups: dict[str, list[int]] = {}
    for idx, grp in enumerate(groups):
        unique_groups.setdefault(grp, []).append(idx)

    between_shared = 0.0
    within_shared = 0.0
    exclusive = 0.0
    for grp, idxs in unique_groups.items():
        ix = np.asarray(idxs, dtype=int)
        ma = float(a[ix].sum())
        mb = float(b[ix].sum())
        if ma > 0.0 and mb > 0.0:
            ta = float((a[ix] @ t[ix]) / ma)
            tb = float((b[ix] @ t[ix]) / mb)
            between_shared += (ma - mb) * (ta + tb) / 2.0
            within_shared += (ta - tb) * (ma + mb) / 2.0
        else:
            exclusive += float(a[ix] @ t[ix] - b[ix] @ t[ix])

    t_target = float(a @ t)
    t_ref = float(b @ t)
    total = float((a - b) @ t)
    recon = between_shared + within_shared + exclusive
    if abs(total - recon) > 1e-10:
        raise RuntimeError(f"Decomposition identity failed: total={total}, recon={recon}")

    unclass_mass = float((a[unclassified_mask].sum() + b[unclassified_mask].sum()) / 2.0)
    return {
        "normalize_to_mapped": bool(normalize_to_mapped),
        "target_weight_sum": float(a.sum()),
        "ref_weight_sum": float(b.sum()),
        "target_mean": t_target,
        "ref_mean": t_ref,
        "total_diff": total,
        "between_shared": between_shared,
        "within_shared": within_shared,
        "exclusive": exclusive,
        "unclassified_mass": unclass_mass,
        "identity_residual": abs(total - recon),
    }


def block_bootstrap_stat(
    blocks: Sequence[str],
    stat_fn: Callable[[np.ndarray], float | np.ndarray],
    *,
    n_boot: int = 1000,
    seed: int = 20261108,
) -> dict[str, np.ndarray | float | int]:
    """Cluster/block bootstrap resampling entire blocks (sites, genera, orders, or clades) with replacement.

    Explicitly rejects empty, null, or whitespace-only block identifiers, and never treats
    individual rows/genes within a block as independent observations.
    """
    if blocks is None or len(blocks) == 0:
        raise ValueError("Cannot bootstrap empty block array")
    clean_blocks: list[str] = []
    for b in blocks:
        if b is None or pd.isna(b):
            raise ValueError("Null/NaN block identifier found in block_bootstrap_stat")
        s = str(b).strip()
        if not s or s in ("nan", "None", "<NA>"):
            raise ValueError(f"Empty or invalid block identifier '{b}' in block_bootstrap_stat")
        clean_blocks.append(s)

    block_arr = np.asarray(clean_blocks, dtype=object)
    unique_blocks, inverse = np.unique(block_arr, return_inverse=True)
    n_blocks = len(unique_blocks)
    if n_blocks < 2:
        raise ValueError(f"Block bootstrap requires at least 2 distinct blocks, got {n_blocks}")

    block_to_indices = [np.where(inverse == i)[0] for i in range(n_blocks)]
    obs = np.asarray(stat_fn(np.arange(len(block_arr))), dtype="float64")
    if not np.isfinite(obs).all():
        raise ValueError("Observed statistic in block_bootstrap_stat is non-finite")

    rng = np.random.default_rng(seed)
    boot_draws = np.empty((n_boot,) + obs.shape, dtype="float64")
    for b_idx in range(n_boot):
        chosen = rng.integers(0, n_blocks, size=n_blocks)
        resampled_idx = np.concatenate([block_to_indices[c] for c in chosen])
        boot_draws[b_idx] = stat_fn(resampled_idx)

    ci_low_95 = np.nanpercentile(boot_draws, 2.5, axis=0)
    ci_high_95 = np.nanpercentile(boot_draws, 97.5, axis=0)
    ci_low_90 = np.nanpercentile(boot_draws, 5.0, axis=0)
    ci_high_90 = np.nanpercentile(boot_draws, 95.0, axis=0)
    se = np.nanstd(boot_draws, ddof=1, axis=0)

    return {
        "n_blocks": int(n_blocks),
        "n_obs": int(len(block_arr)),
        "estimate": obs.item() if obs.ndim == 0 else obs,
        "se": se.item() if se.ndim == 0 else se,
        "ci95_low": ci_low_95.item() if ci_low_95.ndim == 0 else ci_low_95,
        "ci95_high": ci_high_95.item() if ci_high_95.ndim == 0 else ci_high_95,
        "ci90_low": ci_low_90.item() if ci_low_90.ndim == 0 else ci_low_90,
        "ci90_high": ci_high_90.item() if ci_high_90.ndim == 0 else ci_high_90,
        "draws": boot_draws,
    }


def within_contig_circular_null(
    contig_ids: Sequence[str],
    target_mask: np.ndarray | Sequence[bool],
    flank_target_mask: np.ndarray | Sequence[bool],
    window: int = 5,
    *,
    n_perm: int = 1000,
    seed: int = 20261108,
) -> dict[str, float | int]:
    """Test whether target genes co-occur within +/- `window` genes of `flank_target_mask` on the same contig.

    Uses within-contig circular shifts so contig gene counts and local operon clumping of flanks are preserved.
    """
    if window < 1:
        raise ValueError("window must be >= 1")
    if len(contig_ids) == 0:
        raise ValueError("contig_ids cannot be empty")
    clean_contigs: list[str] = []
    for c in contig_ids:
        if c is None or pd.isna(c) or not str(c).strip():
            raise ValueError("Null or empty contig ID in within_contig_circular_null")
        clean_contigs.append(str(c).strip())
    contigs = np.asarray(clean_contigs, dtype=object)
    tgt = np.asarray(target_mask, dtype=bool)
    flk = np.asarray(flank_target_mask, dtype=bool)
    if not (len(contigs) == len(tgt) == len(flk)):
        raise ValueError("Length mismatch in within_contig_circular_null")

    unique_contigs, inv = np.unique(contigs, return_inverse=True)
    active_contigs = []
    for c_idx in range(len(unique_contigs)):
        idxs = np.where(inv == c_idx)[0]
        if tgt[idxs].any():
            active_contigs.append((tgt[idxs], flk[idxs]))

    if not active_contigs:
        return {
            "n_target_loci": 0,
            "obs_adjacent": 0,
            "null_mean": 0.0,
            "null_sd": 0.0,
            "enrichment_ratio": np.nan,
            "p_empirical_one_sided": np.nan,
        }

    def _has_flank_within_window(flk_vec: np.ndarray, w: int) -> np.ndarray:
        n = len(flk_vec)
        out = np.zeros(n, dtype=bool)
        flk_pos = np.where(flk_vec)[0]
        if len(flk_pos) == 0:
            return out
        for i in range(n):
            dists = np.abs(flk_pos - i)
            if np.any((dists >= 1) & (dists <= w)):
                out[i] = True
        return out

    obs_count = 0
    n_loci = 0
    contig_flank_hoods = []
    for t_vec, f_vec in active_contigs:
        hood = _has_flank_within_window(f_vec, window)
        n_t = int(t_vec.sum())
        n_loci += n_t
        obs_count += int((t_vec & hood).sum())
        contig_flank_hoods.append((t_vec, hood))

    rng = np.random.default_rng(seed)
    null_counts = np.zeros(n_perm, dtype=int)
    for p in range(n_perm):
        cnt = 0
        for t_vec, hood in contig_flank_hoods:
            n = len(t_vec)
            shift = int(rng.integers(0, n))
            shifted_t = np.roll(t_vec, shift)
            cnt += int((shifted_t & hood).sum())
        null_counts[p] = cnt

    null_mean = float(null_counts.mean())
    null_sd = float(null_counts.std(ddof=1)) if n_perm > 1 else 0.0
    p_emp = float((1 + np.sum(null_counts >= obs_count)) / (n_perm + 1))
    enrich = float(obs_count / null_mean) if null_mean > 0 else np.nan
    return {
        "n_target_loci": int(n_loci),
        "obs_adjacent": int(obs_count),
        "null_mean": null_mean,
        "null_sd": null_sd,
        "enrichment_ratio": enrich,
        "p_empirical_one_sided": p_emp,
    }


def classify_candidate_role(
    family_id: str,
    *,
    is_established_annotated_control: bool = False,
) -> str:
    """Ensure GrpE/DnaJ/HSP20-associated dark families are labelled as candidate examples, never positive controls."""
    fid = str(family_id).strip()
    if not fid:
        raise ValueError("Empty family_id in classify_candidate_role")
    if fid in PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES:
        if is_established_annotated_control:
            raise ValueError(
                f"{fid} is a prior unknown/dark heat-stress candidate example with unresolved "
                f"genus-aggregated effects, NOT an established annotated positive control."
            )
        return "prior_heat_stress_candidate_example"
    if is_established_annotated_control:
        return "annotated_positive_control"
    return "de_novo_candidate_family"


def evaluate_frequentist_tost(
    estimate: float,
    se: float,
    ci90_low: float,
    ci90_high: float,
    ci95_low: float,
    ci95_high: float,
    delta: float | None = None,
    *,
    margin_justified: bool = False,
    margin_justification: str = "",
) -> dict[str, float | str | bool]:
    """Evaluate frequentist TOST practical equivalence ONLY when an external margin is scientifically justified.

    - Defaults to `margin_justified=False`: reports estimate and CIs with status
      `'no_equivalence_claim_unjustified_margin'`.
    - If `margin_justified=True`, requires a non-empty `margin_justification` and positive `delta`.
    - Never claims equivalence from a null CI that overlaps +/-delta.
    - Must ONLY be applied to frequentist estimates/CIs (e.g. spatial block bootstrap), NEVER Bayesian INLA CIs.
    """
    if not margin_justified:
        return {
            "inference_framework": "frequentist_bootstrap",
            "estimate": float(estimate),
            "se": float(se),
            "ci90_low": float(ci90_low),
            "ci90_high": float(ci90_high),
            "ci95_low": float(ci95_low),
            "ci95_high": float(ci95_high),
            "delta": np.nan if delta is None else float(delta),
            "margin_justified": False,
            "margin_justification": margin_justification,
            "p_tost": np.nan,
            "within_90ci": False,
            "within_95ci": False,
            "status": "no_equivalence_claim_unjustified_margin",
        }

    if not margin_justification or not margin_justification.strip():
        raise ValueError("margin_justified=True requires a non-empty scientific margin_justification string")
    if delta is None or not np.isfinite(delta) or delta <= 0:
        raise ValueError(f"Justified equivalence margin delta must be finite and positive, got {delta}")
    if not (np.isfinite(estimate) and np.isfinite(se) and se > 0):
        return {
            "inference_framework": "frequentist_bootstrap",
            "estimate": float(estimate),
            "se": float(se),
            "ci90_low": float(ci90_low),
            "ci90_high": float(ci90_high),
            "ci95_low": float(ci95_low),
            "ci95_high": float(ci95_high),
            "delta": float(delta),
            "margin_justified": True,
            "margin_justification": margin_justification.strip(),
            "p_tost": np.nan,
            "within_90ci": False,
            "within_95ci": False,
            "status": "uncertain_degenerate_se",
        }

    z_lower = (estimate - (-delta)) / se
    z_upper = (delta - estimate) / se
    p_lower = float(stats.norm.sf(z_lower))
    p_upper = float(stats.norm.sf(z_upper))
    p_tost = max(p_lower, p_upper)

    within_90 = bool((ci90_low > -delta) and (ci90_high < delta))
    within_95 = bool((ci95_low > -delta) and (ci95_high < delta))
    excludes_zero_95 = bool((ci95_low > 0.0) or (ci95_high < 0.0))
    exceeds_delta_95 = bool((ci95_low >= delta) or (ci95_high <= -delta))

    if within_95 and p_tost < 0.05:
        status = "equivalent_within_bounds_strict95"
    elif within_90 and p_tost < 0.05:
        status = "equivalent_within_bounds_tost90"
    elif exceeds_delta_95:
        status = "exceeds_equivalence_bounds"
    elif excludes_zero_95:
        status = "non_zero_but_overlaps_equivalence_bound"
    else:
        status = "uncertain_inconclusive"

    return {
        "inference_framework": "frequentist_bootstrap",
        "estimate": float(estimate),
        "se": float(se),
        "ci90_low": float(ci90_low),
        "ci90_high": float(ci90_high),
        "ci95_low": float(ci95_low),
        "ci95_high": float(ci95_high),
        "delta": float(delta),
        "margin_justified": True,
        "margin_justification": margin_justification.strip(),
        "p_tost": float(p_tost),
        "within_90ci": within_90,
        "within_95ci": within_95,
        "status": status,
    }


def evaluate_bayesian_rope(
    posterior_mean: float,
    posterior_sd: float,
    cred95_low: float,
    cred95_high: float,
    rope_halfwidth: float | None = None,
    *,
    rope_justified: bool = False,
    rope_justification: str = "",
) -> dict[str, float | str | bool]:
    """Summarize Bayesian INLA posterior intervals without conflating them with frequentist TOST.

    - When `rope_justified=False` (default): reports posterior mean, SD, and 95% credible interval
      with status `'no_equivalence_claim_unjustified_rope'` (never claims equivalence from a credible
      interval crossing zero).
    - When `rope_justified=True`: requires a non-empty `rope_justification` and positive `rope_halfwidth`,
      and evaluates normal-approximated posterior mass inside `[-rope_halfwidth, +rope_halfwidth]` and
      95% credible interval inclusion.
    """
    excludes_zero = bool((cred95_low > 0.0) or (cred95_high < 0.0))
    if not rope_justified:
        return {
            "inference_framework": "bayesian_inla_posterior",
            "posterior_mean": float(posterior_mean),
            "posterior_sd": float(posterior_sd),
            "cred95_low": float(cred95_low),
            "cred95_high": float(cred95_high),
            "cred95_excludes_zero": excludes_zero,
            "rope_halfwidth": np.nan if rope_halfwidth is None else float(rope_halfwidth),
            "rope_justified": False,
            "rope_justification": rope_justification,
            "posterior_mass_in_rope": np.nan,
            "cred95_inside_rope": False,
            "status": "no_equivalence_claim_unjustified_rope",
        }

    if not rope_justification or not rope_justification.strip():
        raise ValueError("rope_justified=True requires a non-empty rope_justification string")
    if rope_halfwidth is None or not np.isfinite(rope_halfwidth) or rope_halfwidth <= 0:
        raise ValueError(f"ROPE halfwidth must be finite and positive, got {rope_halfwidth}")
    if not (np.isfinite(posterior_mean) and np.isfinite(posterior_sd) and posterior_sd > 0):
        return {
            "inference_framework": "bayesian_inla_posterior",
            "posterior_mean": float(posterior_mean),
            "posterior_sd": float(posterior_sd),
            "cred95_low": float(cred95_low),
            "cred95_high": float(cred95_high),
            "cred95_excludes_zero": excludes_zero,
            "rope_halfwidth": float(rope_halfwidth),
            "rope_justified": True,
            "rope_justification": rope_justification.strip(),
            "posterior_mass_in_rope": np.nan,
            "cred95_inside_rope": False,
            "status": "uncertain_degenerate_posterior_sd",
        }

    mass_in_rope = float(
        stats.norm.cdf(rope_halfwidth, loc=posterior_mean, scale=posterior_sd)
        - stats.norm.cdf(-rope_halfwidth, loc=posterior_mean, scale=posterior_sd)
    )
    cred95_inside = bool((cred95_low > -rope_halfwidth) and (cred95_high < rope_halfwidth))
    cred95_outside = bool((cred95_low >= rope_halfwidth) or (cred95_high <= -rope_halfwidth))

    if cred95_inside and mass_in_rope >= 0.95:
        status = "bayesian_cred95_inside_rope"
    elif cred95_outside:
        status = "bayesian_cred95_outside_rope"
    else:
        status = "bayesian_uncertain_overlaps_rope"

    return {
        "inference_framework": "bayesian_inla_posterior",
        "posterior_mean": float(posterior_mean),
        "posterior_sd": float(posterior_sd),
        "cred95_low": float(cred95_low),
        "cred95_high": float(cred95_high),
        "cred95_excludes_zero": excludes_zero,
        "rope_halfwidth": float(rope_halfwidth),
        "rope_justified": True,
        "rope_justification": rope_justification.strip(),
        "posterior_mass_in_rope": mass_in_rope,
        "cred95_inside_rope": cred95_inside,
        "status": status,
    }


def bh_fdr(pvals: np.ndarray | Sequence[float]) -> np.ndarray:
    """Exact Benjamini-Hochberg FDR q-values without requiring statsmodels."""
    p = np.asarray(pvals, dtype="float64")
    if p.ndim != 1:
        raise ValueError("pvals must be 1D")
    n = len(p)
    if n == 0:
        return np.array([], dtype="float64")
    if not np.isfinite(p).all() or (p < 0.0).any() or (p > 1.0).any():
        raise ValueError("pvals must be finite and in [0, 1]")
    order = np.argsort(p)
    ranked = p[order]
    q_sorted = ranked * n / np.arange(1, n + 1, dtype="float64")
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q = np.empty(n, dtype="float64")
    q[order] = q_sorted
    return q


def fit_clustered_binomial_glm(
    y_pos: np.ndarray | Sequence[float],
    n_total: np.ndarray | Sequence[float],
    X: np.ndarray,
    cluster_ids: Sequence[str],
    *,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> dict[str, np.ndarray | int]:
    """Fit grouped Binomial logistic GLM via IRLS with Liang-Zeger cluster-robust sandwich standard errors."""
    y = np.asarray(y_pos, dtype="float64")
    n = np.asarray(n_total, dtype="float64")
    X_mat = np.asarray(X, dtype="float64")
    if not (y.ndim == 1 and n.ndim == 1 and X_mat.ndim == 2):
        raise ValueError("Invalid dimensions in fit_clustered_binomial_glm")
    if not (len(y) == len(n) == X_mat.shape[0] == len(cluster_ids)):
        raise ValueError("Row length mismatch in fit_clustered_binomial_glm")
    if (n <= 0).any() or (y < 0).any() or (y > n).any():
        raise ValueError("Invalid binomial counts in fit_clustered_binomial_glm")

    k = X_mat.shape[1]
    beta = np.zeros(k, dtype="float64")
    for _ in range(max_iter):
        eta = np.clip(X_mat @ beta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(n * mu * (1.0 - mu), 1e-12)
        score = X_mat.T @ (y - n * mu)
        H = X_mat.T @ (X_mat * w[:, None])
        step = np.linalg.solve(H + 1e-10 * np.eye(k), score)
        beta = beta + step
        if float(np.max(np.abs(step))) < tol:
            break

    eta = np.clip(X_mat @ beta, -30.0, 30.0)
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = np.maximum(n * mu * (1.0 - mu), 1e-12)
    bread = np.linalg.inv(X_mat.T @ (X_mat * w[:, None]) + 1e-10 * np.eye(k))

    resid = y - n * mu
    unique_clusters, inv = np.unique(np.asarray(cluster_ids, dtype=object), return_inverse=True)
    G = len(unique_clusters)
    if G < 2:
        raise ValueError("Clustered sandwich SE requires at least 2 clusters")
    meat = np.zeros((k, k), dtype="float64")
    for c_idx in range(G):
        m = inv == c_idx
        u_c = X_mat[m].T @ resid[m]
        meat += np.outer(u_c, u_c)

    # Small-sample cluster adjustment matching statsmodels: G / (G - 1)
    cov_cluster = (G / float(G - 1)) * (bread @ meat @ bread)
    se = np.sqrt(np.maximum(np.diag(cov_cluster), 1e-16))
    z = beta / se
    pvals = 2.0 * stats.norm.sf(np.abs(z))
    return {
        "params": beta,
        "bse": se,
        "zvalues": z,
        "pvalues": pvals,
        "n_clusters": int(G),
    }


CODON_TABLE_11: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGC": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
TABLE_11_START_CODONS = frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"})
DNA_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement_dna(seq: str) -> str:
    """Return reverse-complement of a DNA string."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def _require_exact_int_coord(val: object) -> int | None:
    """Validate that a genomic coordinate is an exact integer (never truncating floats or booleans)."""
    if val is None or isinstance(val, bool) or pd.isna(val):
        return None
    if isinstance(val, (int, np.integer)):
        return int(val)
    if isinstance(val, (float, np.floating)):
        v = float(val)
        if not np.isfinite(v) or v != math.floor(v):
            return None
        return int(v)
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return None
        if s.startswith(("+", "-")):
            sign_part, digits_part = s[0], s[1:]
        else:
            sign_part, digits_part = "", s
        if not digits_part.isdigit():
            return None
        return int(sign_part + digits_part)
    return None


def clean_candidate_aa_sequence(raw_aa: object) -> dict[str, str | int | bool]:
    """Validate a candidate AA sequence string without deleting internal '*' stops.

    Permits stripping at most ONE terminal '*' stop character. Any remaining '*' inside
    the sequence body is flagged as an internal stop (`status='invalid_internal_stop_in_candidate_aa'`).
    Missing/empty sequences return `status='missing_candidate_sequence'`.
    """
    if raw_aa is None or isinstance(raw_aa, bool) or pd.isna(raw_aa):
        return {
            "status": "missing_candidate_sequence",
            "aa_body": "",
            "had_terminal_stop": False,
            "n_internal_stops_in_candidate": -1,
        }
    s = str(raw_aa).strip()
    if not s or s in {"missing", "nan", "None", "<NA>"}:
        return {
            "status": "missing_candidate_sequence",
            "aa_body": "",
            "had_terminal_stop": False,
            "n_internal_stops_in_candidate": -1,
        }
    if any(ch.isspace() for ch in s):
        return {
            "status": "invalid_whitespace_in_candidate_aa",
            "aa_body": "",
            "had_terminal_stop": False,
            "n_internal_stops_in_candidate": -1,
        }
    had_terminal_stop = s.endswith("*")
    body = s[:-1] if had_terminal_stop else s
    n_internal = int(body.count("*"))
    if len(body) == 0:
        return {
            "status": "missing_candidate_sequence",
            "aa_body": "",
            "had_terminal_stop": had_terminal_stop,
            "n_internal_stops_in_candidate": n_internal,
        }
    if n_internal > 0:
        return {
            "status": "invalid_internal_stop_in_candidate_aa",
            "aa_body": body,
            "had_terminal_stop": had_terminal_stop,
            "n_internal_stops_in_candidate": n_internal,
        }
    return {
        "status": "ok",
        "aa_body": body,
        "had_terminal_stop": had_terminal_stop,
        "n_internal_stops_in_candidate": 0,
    }


def translate_cds_table11(
    cds_dna: str,
    *,
    partial_flag: str = "00",
) -> dict[str, str | int | bool]:
    """Translate a CDS nucleotide string using NCBI Genetic Code 11 (Bacterial/Archaeal).

    Handles Prodigal start-codon semantics: when the 5' end is complete (`partial_flag` starts
    with '0', e.g. '00' or '01') and the first codon is a recognized Table 11 start codon
    (ATG, GTG, TTG, CTG, ATT, ATC, ATA), the initiating residue is translated as 'M'.
    When 5'-partial (`partial_flag` starts with '1', e.g. '10' or '11'), the first codon is
    translated literally. Strips a single terminal stop codon ('*') if present at the 3' end,
    records `has_admissible_start_codon`, and counts any internal stop codons separately.
    """
    dna = cds_dna.strip().upper()
    if len(dna) < 3 or (len(dna) % 3) != 0:
        return {
            "status": "non_triplet_or_empty_cds",
            "aa_sequence": "",
            "has_admissible_start_codon": False,
            "has_terminal_stop": False,
            "n_internal_stops": -1,
            "start_codon": dna[:3] if len(dna) >= 3 else "",
        }

    p_str = str(partial_flag).strip()
    if p_str in ("0", "False", "false", "0.0"):
        p_str = "00"
    five_prime_complete = p_str.startswith("0")

    codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
    has_admissible_start = bool(codons[0] in TABLE_11_START_CODONS)
    aa_chars: list[str] = []
    for idx, cod in enumerate(codons):
        if idx == 0 and five_prime_complete and has_admissible_start:
            aa_chars.append("M")
        else:
            aa_chars.append(CODON_TABLE_11.get(cod, "X"))

    has_terminal_stop = bool(aa_chars and aa_chars[-1] == "*")
    if has_terminal_stop:
        aa_body = aa_chars[:-1]
    else:
        aa_body = aa_chars

    n_internal_stops = int(sum(1 for c in aa_body if c == "*"))
    aa_seq = "".join(aa_body)
    return {
        "status": "ok",
        "aa_sequence": aa_seq,
        "has_admissible_start_codon": has_admissible_start,
        "has_terminal_stop": has_terminal_stop,
        "n_internal_stops": n_internal_stops,
        "start_codon": codons[0],
    }


def verify_locus_against_refined_contig(
    contigs_dict: dict[str, str],
    contig_id: str,
    start: object,
    end: object,
    strand: int | str,
    partial_flag: str,
    candidate_aa_seq: object,
) -> dict[str, str | int | bool]:
    """Independently verify contig retention, integer coordinate bounds, strand, start/stop codons, and CDS translation."""
    if not contig_id or contig_id not in contigs_dict:
        return {
            "refined_contig_present": False,
            "refined_contig_bp_len": 0,
            "locus_within_contig_bounds": False,
            "has_admissible_start_codon": False,
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": -1,
            "n_internal_stops_in_candidate": -1,
            "has_terminal_stop_codon": False,
            "verification_status": "contig_missing_in_refined_fasta",
            "final_promotion_gate_pass": False,
        }

    contig_dna = contigs_dict[contig_id]
    c_len = len(contig_dna)
    s_i = _require_exact_int_coord(start)
    e_i = _require_exact_int_coord(end)
    if s_i is None or e_i is None:
        return {
            "refined_contig_present": True,
            "refined_contig_bp_len": c_len,
            "locus_within_contig_bounds": False,
            "has_admissible_start_codon": False,
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": -1,
            "n_internal_stops_in_candidate": -1,
            "has_terminal_stop_codon": False,
            "verification_status": "invalid_noninteger_or_null_coordinates",
            "final_promotion_gate_pass": False,
        }

    nt_len = e_i - s_i + 1
    in_bounds = (s_i >= 1) and (e_i > s_i) and (e_i <= c_len) and ((nt_len % 3) == 0)
    if not in_bounds:
        return {
            "refined_contig_present": True,
            "refined_contig_bp_len": c_len,
            "locus_within_contig_bounds": False,
            "has_admissible_start_codon": False,
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": -1,
            "n_internal_stops_in_candidate": -1,
            "has_terminal_stop_codon": False,
            "verification_status": "out_of_bounds_or_non_triplet",
            "final_promotion_gate_pass": False,
        }

    strand_s = str(strand).strip()
    if strand_s in ("1", "+1", "+"):
        cds_slice = contig_dna[s_i - 1 : e_i]
    elif strand_s in ("-1", "-"):
        cds_slice = reverse_complement_dna(contig_dna[s_i - 1 : e_i])
    else:
        return {
            "refined_contig_present": True,
            "refined_contig_bp_len": c_len,
            "locus_within_contig_bounds": True,
            "has_admissible_start_codon": False,
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": -1,
            "n_internal_stops_in_candidate": -1,
            "has_terminal_stop_codon": False,
            "verification_status": "invalid_strand",
            "final_promotion_gate_pass": False,
        }

    p_str = str(partial_flag).strip()
    if p_str in ("0", "False", "false", "0.0"):
        p_str = "00"
    trans = translate_cds_table11(cds_slice, partial_flag=p_str)
    cand_info = clean_candidate_aa_sequence(candidate_aa_seq)

    if cand_info["status"] == "missing_candidate_sequence":
        return {
            "refined_contig_present": True,
            "refined_contig_bp_len": c_len,
            "locus_within_contig_bounds": True,
            "has_admissible_start_codon": bool(trans["has_admissible_start_codon"]),
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": int(trans["n_internal_stops"]),
            "n_internal_stops_in_candidate": -1,
            "has_terminal_stop_codon": bool(trans["has_terminal_stop"]),
            "verification_status": "unknown_missing_candidate_sequence",
            "final_promotion_gate_pass": False,
        }

    if cand_info["status"] != "ok":
        return {
            "refined_contig_present": True,
            "refined_contig_bp_len": c_len,
            "locus_within_contig_bounds": True,
            "has_admissible_start_codon": bool(trans["has_admissible_start_codon"]),
            "cds_translation_exact_match": False,
            "n_internal_stops_translated": int(trans["n_internal_stops"]),
            "n_internal_stops_in_candidate": int(cand_info["n_internal_stops_in_candidate"]),
            "has_terminal_stop_codon": bool(trans["has_terminal_stop"]),
            "verification_status": str(cand_info["status"]),
            "final_promotion_gate_pass": False,
        }

    cand_clean = str(cand_info["aa_body"])
    exact_match = bool(
        trans["status"] == "ok"
        and len(cand_clean) > 0
        and trans["aa_sequence"] == cand_clean
    )
    has_adm_start = bool(trans["has_admissible_start_codon"])
    has_term_stop = bool(trans["has_terminal_stop"])
    no_internal_stops = (int(trans["n_internal_stops"]) == 0) and (int(cand_info["n_internal_stops_in_candidate"]) == 0)
    is_complete_orf = (p_str == "00") and has_adm_start and has_term_stop and no_internal_stops

    if not exact_match:
        v_status = "mismatch_sequence_or_strand"
    elif not no_internal_stops:
        v_status = "internal_stop_in_translated_cds"
    elif p_str != "00":
        v_status = "verified_translation_but_partial_orf"
    elif not has_adm_start:
        v_status = "inadmissible_start_codon"
    elif not has_term_stop:
        v_status = "missing_terminal_stop_codon"
    else:
        v_status = "verified_exact_complete_cds"

    return {
        "refined_contig_present": True,
        "refined_contig_bp_len": c_len,
        "locus_within_contig_bounds": True,
        "has_admissible_start_codon": has_adm_start,
        "cds_translation_exact_match": exact_match,
        "n_internal_stops_translated": int(trans["n_internal_stops"]),
        "n_internal_stops_in_candidate": int(cand_info["n_internal_stops_in_candidate"]),
        "has_terminal_stop_codon": has_term_stop,
        "verification_status": v_status,
        "final_promotion_gate_pass": bool(exact_match and is_complete_orf),
    }


def parse_and_verify_refined_fasta_bytes(
    genome_id: str,
    raw_bytes: bytes,
    *,
    expected_sha256: str | None = None,
    expected_contigs: int | None = None,
    expected_bases: int | None = None,
) -> dict[str, object]:
    """Parse a refined genome FASTA byte buffer, rejecting duplicate contig IDs and verifying against `genome_manifest.tsv`."""
    computed_sha256 = hashlib.sha256(raw_bytes).hexdigest()
    text = raw_bytes.decode("ascii")
    contigs: dict[str, str] = {}
    cur_cid: str | None = None
    chunks: list[str] = []

    for raw_line in text.splitlines():
        if not raw_line:
            continue
        if raw_line.startswith(">"):
            if cur_cid is not None:
                seq_str = "".join(chunks)
                if not seq_str:
                    raise ValueError(f"{genome_id}: empty FASTA sequence for contig '{cur_cid}'")
                contigs[cur_cid] = seq_str
            fields = raw_line[1:].strip().split()
            if not fields:
                raise ValueError(f"{genome_id}: empty FASTA header line")
            hdr = fields[0]
            cid = hdr.split("|", 1)[1] if "|" in hdr else hdr
            if cid in contigs or cid == cur_cid:
                raise ValueError(f"{genome_id}: duplicate contig ID '{cid}' in FASTA")
            cur_cid = cid
            chunks = []
        else:
            if cur_cid is None:
                raise ValueError(f"{genome_id}: sequence data before first FASTA header")
            if any(ch.isspace() for ch in raw_line):
                raise ValueError(f"{genome_id}: whitespace inside sequence line on contig '{cur_cid}'")
            chunks.append(raw_line)

    if cur_cid is not None:
        seq_str = "".join(chunks)
        if not seq_str:
            raise ValueError(f"{genome_id}: empty FASTA sequence for contig '{cur_cid}'")
        contigs[cur_cid] = seq_str

    if not contigs:
        raise ValueError(f"{genome_id}: FASTA contains 0 contigs")

    n_contigs = len(contigs)
    total_bp = int(sum(len(s) for s in contigs.values()))

    if expected_sha256 is not None:
        sha_ok = (computed_sha256 == expected_sha256.strip().lower())
        contigs_ok = (expected_contigs is None) or (n_contigs == int(expected_contigs))
        bases_ok = (expected_bases is None) or (total_bp == int(expected_bases))
        if sha_ok and contigs_ok and bases_ok:
            prov_status = "verified_exact_match_genome_manifest"
        else:
            prov_status = "mismatch_against_genome_manifest"
    else:
        sha_ok = False
        contigs_ok = False
        bases_ok = False
        prov_status = "unverified_freshly_recorded_hash"

    return {
        "genome": genome_id,
        "computed_sha256": computed_sha256,
        "manifest_expected_sha256": expected_sha256 if expected_sha256 is not None else "missing_manifest",
        "sha256_matches_manifest": bool(sha_ok),
        "n_refined_contigs": n_contigs,
        "manifest_expected_contigs": int(expected_contigs) if expected_contigs is not None else -1,
        "contigs_match_manifest": bool(contigs_ok),
        "total_refined_bp": total_bp,
        "manifest_expected_bp": int(expected_bases) if expected_bases is not None else -1,
        "bp_matches_manifest": bool(bases_ok),
        "provenance_status": prov_status,
        "contigs_dict": contigs,
    }


def scale_plddt_explicit_source(val: object, source_scale: str) -> float:
    """Convert pLDDT to 0..1 using an explicit declared source scale ('0_to_1' or '0_to_100').

    Never guesses scale per row. Raises ValueError if `source_scale` is invalid or if `val`
    falls outside the declared bounded range.
    """
    if source_scale not in ("0_to_1", "0_to_100"):
        raise ValueError(f"Invalid pLDDT source_scale '{source_scale}'; must be '0_to_1' or '0_to_100'")
    if val is None or isinstance(val, bool) or pd.isna(val):
        return np.nan
    v = float(val)
    if not np.isfinite(v):
        return np.nan
    if source_scale == "0_to_1":
        if v < 0.0 or v > 1.0:
            raise ValueError(f"pLDDT value {v} out of bounds [0.0, 1.0] for source_scale='0_to_1'")
        return v
    if v < 0.0 or v > 100.0:
        raise ValueError(f"pLDDT value {v} out of bounds [0.0, 100.0] for source_scale='0_to_100'")
    return v / 100.0


def identify_mgnify_gtdb_duplicates(
    comp_df: pd.DataFrame,
    skani_df: pd.DataFrame,
    *,
    ani_threshold: float = 99.9,
) -> set[str]:
    """Identify MGnify comparator genomes that duplicate GTDB r232 species representatives.

    Combines exact assembly accession matches (`same_assembly_as_gtdb_rep == 'yes'`, n=41)
    and skani ANI >= 99.9% pairs between an MGnify genome and a GTDB representative (n=54 total).
    """
    validate_primary_key(comp_df, ["genome_id"], "v2_thermo_relatives_comparators")
    gtdb_ids = set(comp_df.loc[comp_df["source"] == "GTDB_r232_rep", "genome_id"].astype(str))
    mgn_ids = set(comp_df.loc[comp_df["source"] == "MGnify", "genome_id"].astype(str))

    same_asm = set(
        comp_df.loc[
            comp_df["same_assembly_as_gtdb_rep"].astype(str).str.strip().str.lower() == "yes",
            "genome_id",
        ].astype(str)
    )

    def _path_to_gid(p: object) -> str:
        name = Path(str(p)).name
        for suf in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
            if name.endswith(suf):
                return name[: -len(suf)]
        return name

    r_ids = skani_df["Ref_file"].map(_path_to_gid)
    q_ids = skani_df["Query_file"].map(_path_to_gid)
    anis = pd.to_numeric(skani_df["ANI"], errors="coerce")

    sk_dups: set[str] = set()
    for r, q, ani in zip(r_ids, q_ids, anis):
        if pd.notna(ani) and float(ani) >= ani_threshold and r != q:
            if r in mgn_ids and q in gtdb_ids:
                sk_dups.add(r)
            elif q in mgn_ids and r in gtdb_ids:
                sk_dups.add(q)

    return same_asm | sk_dups


def estimate_matched_joint_log_contrast(
    df: pd.DataFrame,
    *,
    tax_col: str,
    func_col: str,
    dune_col: str = "edgell_dune_type",
    ref_dune: str = "uruq",
    target_dune: str = "megabarchans",
    compartment_col: str = "compartment",
    ref_compartment: str = "Deep",
    trip_col: str = "trip",
    ref_trip: str = "trip1",
    log_reads_col: str = "log_reads",
    site_col: str = "site",
    block_col: str = "spatial_block",
    weighting: str = "sample_weighted",
    n_boot: int = 1000,
    seed: int = 20261108,
    custom_resample_indices: list[np.ndarray] | None = None,
) -> dict[str, object]:
    """Estimate the matched covariate-adjusted joint log contrast (`target_dune` vs `ref_dune`)
    for taxonomic Hill1 (`tax_col`) and functional Hill1 (`func_col`) on identical spatial-block resamples.

    Preserves the ecological estimand (`edgell_dune_type`: `megabarchans` vs `uruq`) and fixed-effect
    covariate design (`~ 1 + compartment + trip + log_reads + edgell_dune_type`) used in `fit.R`.
    Both endpoints are transformed to `log(Hill1)` (comparable multiplicative natural-log units),
    and the direct difference of the two log contrasts (`delta_log = beta_tax - beta_func`) is
    estimated jointly on each spatial block bootstrap replicate.
    Fixed base weights `1 / original_observations_per_site` are computed once on the original matched
    sample set and indexed by `idx` inside each resample so repeated block copies contribute repeated weight.
    Zero-cell or rank-deficient bootstrap replicates are recorded as invalid (`n_invalid_replicates`)
    and excluded rather than silently returning NaNs as completed inference.
    """
    if weighting not in ("sample_weighted", "site_balanced", "site_averaged"):
        raise ValueError(f"Unsupported weighting='{weighting}'")

    req_cols = [tax_col, func_col, dune_col, compartment_col, trip_col, log_reads_col, site_col, block_col]
    for c in req_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column '{c}' in estimate_matched_joint_log_contrast")

    work = df.copy()
    y_tax_raw = pd.to_numeric(work[tax_col], errors="coerce")
    y_func_raw = pd.to_numeric(work[func_col], errors="coerce")
    if np.isinf(y_tax_raw.to_numpy(dtype="float64")).any() or np.isinf(y_func_raw.to_numpy(dtype="float64")).any():
        raise ValueError("Non-finite (+inf/-inf) Hill1 values encountered; reject non-finite positive Hill values at input validation")
    l_reads = pd.to_numeric(work[log_reads_col], errors="coerce")
    valid_mask = (
        y_tax_raw.notna()
        & np.isfinite(y_tax_raw)
        & (y_tax_raw > 0)
        & y_func_raw.notna()
        & np.isfinite(y_func_raw)
        & (y_func_raw > 0)
        & l_reads.notna()
        & np.isfinite(l_reads)
        & work[dune_col].notna()
        & work[compartment_col].notna()
        & work[trip_col].notna()
        & work[site_col].notna()
        & work[block_col].notna()
    )
    work = work.loc[valid_mask].copy().reset_index(drop=True)
    if len(work) < 20:
        return {
            "status": "insufficient_matched_samples",
            "n_samples": len(work),
            "n_sites": int(work[site_col].nunique()) if not work.empty else 0,
            "n_blocks": 0,
            "n_boot_requested": n_boot,
            "n_valid_replicates": 0,
            "n_invalid_replicates": n_boot,
        }

    work["_log_tax"] = np.log(pd.to_numeric(work[tax_col], errors="coerce").to_numpy(dtype="float64"))
    work["_log_func"] = np.log(pd.to_numeric(work[func_col], errors="coerce").to_numpy(dtype="float64"))
    work["_log_reads"] = pd.to_numeric(work[log_reads_col], errors="coerce").to_numpy(dtype="float64")

    if weighting == "site_averaged":
        # Average within (site, compartment, edgell_dune_type, spatial_block) across repeat visits
        work = (
            work.groupby([site_col, compartment_col, dune_col, block_col], as_index=False, sort=False)
            .agg(
                _log_tax=("_log_tax", "mean"),
                _log_func=("_log_func", "mean"),
                _log_reads=("_log_reads", "mean"),
            )
            .reset_index(drop=True)
        )
        include_trip = False
    else:
        include_trip = True

    # Determine fixed design levels present in the full matched dataset
    comp_levels = [c for c in sorted(work[compartment_col].astype(str).unique()) if c != ref_compartment]
    trip_levels = (
        [t for t in sorted(work[trip_col].astype(str).unique()) if t != ref_trip]
        if include_trip
        else []
    )
    dune_levels = [d for d in sorted(work[dune_col].astype(str).unique()) if d != ref_dune]
    if target_dune not in dune_levels or (work[dune_col].astype(str) == ref_dune).sum() == 0:
        return {
            "status": "missing_target_or_reference_dune_in_cohort",
            "n_samples": len(work),
            "n_sites": int(work[site_col].nunique()),
            "n_blocks": int(work[block_col].nunique()),
            "n_boot_requested": n_boot,
            "n_valid_replicates": 0,
            "n_invalid_replicates": n_boot,
        }

    # Build full design matrix X on `work`
    col_names = ["intercept"]
    cols_list = [np.ones(len(work), dtype="float64")]
    for cl in comp_levels:
        col_names.append(f"comp_{cl}")
        cols_list.append((work[compartment_col].astype(str) == cl).to_numpy(dtype="float64"))
    for tl in trip_levels:
        col_names.append(f"trip_{tl}")
        cols_list.append((work[trip_col].astype(str) == tl).to_numpy(dtype="float64"))
    col_names.append("log_reads_z")
    lr_arr = work["_log_reads"].to_numpy(dtype="float64")
    lr_sd = float(np.std(lr_arr, ddof=1))
    cols_list.append((lr_arr - float(np.mean(lr_arr))) / (lr_sd if lr_sd > 1e-12 else 1.0))
    for dl in dune_levels:
        col_names.append(f"dune_{dl}")
        cols_list.append((work[dune_col].astype(str) == dl).to_numpy(dtype="float64"))

    X_full = np.column_stack(cols_list)
    target_col_idx = col_names.index(f"dune_{target_dune}")
    Y_full = np.column_stack([work["_log_tax"].to_numpy(dtype="float64"), work["_log_func"].to_numpy(dtype="float64")])
    sites_arr = work[site_col].astype(str).to_numpy()
    dunes_arr = work[dune_col].astype(str).to_numpy()

    # Precompute fixed base weights 1 / original_observations_per_site on the full sample set
    # so repeated block copies in a bootstrap resample retain their base weights and contribute repeated weight!
    _, orig_inv, orig_counts = np.unique(sites_arr, return_inverse=True, return_counts=True)
    fixed_base_weights = 1.0 / orig_counts[orig_inv].astype("float64")

    def _fit_slice(idx: np.ndarray) -> np.ndarray | None:
        # Check both ref_dune and target_dune are present in this sample/resample
        d_sub = dunes_arr[idx]
        if not ((d_sub == ref_dune).any() and (d_sub == target_dune).any()):
            return None
        X_s = X_full[idx]
        Y_s = Y_full[idx]
        p = X_s.shape[1]
        if X_s.shape[0] <= p:
            return None
        if weighting == "site_balanced":
            w = fixed_base_weights[idx]
            sw = np.sqrt(w)[:, None]
            Xw = X_s * sw
            Yw = Y_s * sw
        else:
            Xw = X_s
            Yw = Y_s

        if np.linalg.matrix_rank(Xw) < p:
            return None
        betas, _, rank, _ = np.linalg.lstsq(Xw, Yw, rcond=None)
        if rank < p:
            return None
        b_tax = float(betas[target_col_idx, 0])
        b_func = float(betas[target_col_idx, 1])
        if not (np.isfinite(b_tax) and np.isfinite(b_func)):
            return None
        return np.array([b_tax, b_func, b_tax - b_func], dtype="float64")

    full_est = _fit_slice(np.arange(len(work), dtype=np.int64))
    if full_est is None:
        return {
            "status": "rank_deficient_full_design",
            "n_samples": len(work),
            "n_sites": int(work[site_col].nunique()),
            "n_blocks": int(work[block_col].nunique()),
            "n_boot_requested": n_boot,
            "n_valid_replicates": 0,
            "n_invalid_replicates": n_boot,
        }

    clean_blocks: list[str] = []
    for b in work[block_col].tolist():
        if b is None or pd.isna(b):
            raise ValueError("Null/NaN block identifier found in estimate_matched_joint_log_contrast")
        s = str(b).strip()
        if not s or s in ("nan", "None", "<NA>"):
            raise ValueError(f"Empty or invalid block identifier '{b}' in estimate_matched_joint_log_contrast")
        clean_blocks.append(s)
    block_arr = np.asarray(clean_blocks, dtype=object)
    unique_blocks = sorted(set(clean_blocks))
    k_blocks = len(unique_blocks)
    if k_blocks < 2:
        raise ValueError(f"Spatial block bootstrap requires at least 2 distinct blocks, got {k_blocks}")
    block_to_indices = {b: np.flatnonzero(block_arr == b) for b in unique_blocks}

    valid_reps: list[np.ndarray] = []
    n_invalid = 0
    if custom_resample_indices is not None:
        n_boot_actual = len(custom_resample_indices)
        for boot_idx in custom_resample_indices:
            res = _fit_slice(np.asarray(boot_idx, dtype=np.int64))
            if res is None:
                n_invalid += 1
            else:
                valid_reps.append(res)
    else:
        n_boot_actual = n_boot
        rng = np.random.default_rng(seed)
        for _ in range(n_boot):
            sampled_b = rng.choice(unique_blocks, size=k_blocks, replace=True)
            boot_idx = np.concatenate([block_to_indices[b] for b in sampled_b])
            res = _fit_slice(boot_idx)
            if res is None:
                n_invalid += 1
            else:
                valid_reps.append(res)

    if len(valid_reps) < max(2 if custom_resample_indices is not None else 10, n_boot_actual // 5):
        return {
            "status": "insufficient_valid_bootstrap_replicates",
            "resampling_design": "spatial_block_bootstrap_conditional_on_target_dune_presence_not_replacement_for_INLA",
            "n_samples": len(work),
            "n_sites": int(work[site_col].nunique()),
            "n_blocks": k_blocks,
            "n_boot_requested": n_boot_actual,
            "n_valid_replicates": len(valid_reps),
            "n_invalid_replicates": n_invalid,
            "beta_log_tax": float(full_est[0]),
            "beta_log_func": float(full_est[1]),
            "delta_log_tax_minus_func": float(full_est[2]),
        }

    boot_mat = np.vstack(valid_reps)
    se_vec = np.std(boot_mat, axis=0, ddof=1)
    lo95_vec = np.quantile(boot_mat, 0.025, axis=0)
    hi95_vec = np.quantile(boot_mat, 0.975, axis=0)

    return {
        "status": "completed",
        "resampling_design": "spatial_block_bootstrap_conditional_on_target_dune_presence_not_replacement_for_INLA",
        "weighting": weighting,
        "n_samples": int(len(work)),
        "n_sites": int(work[site_col].nunique()),
        "n_megabarchan_samples": int((dunes_arr == target_dune).sum()),
        "n_uruq_samples": int((dunes_arr == ref_dune).sum()),
        "n_blocks": int(k_blocks),
        "n_boot_requested": int(n_boot_actual),
        "n_valid_replicates": int(len(valid_reps)),
        "n_invalid_replicates": int(n_invalid),
        "beta_log_tax": float(full_est[0]),
        "beta_log_tax_se": float(se_vec[0]),
        "beta_log_tax_ci95_low": float(lo95_vec[0]),
        "beta_log_tax_ci95_high": float(hi95_vec[0]),
        "beta_log_func": float(full_est[1]),
        "beta_log_func_se": float(se_vec[1]),
        "beta_log_func_ci95_low": float(lo95_vec[1]),
        "beta_log_func_ci95_high": float(hi95_vec[1]),
        "delta_log_tax_minus_func": float(full_est[2]),
        "delta_log_se": float(se_vec[2]),
        "delta_log_ci95_low": float(lo95_vec[2]),
        "delta_log_ci95_high": float(hi95_vec[2]),
        "fold_ratio_tax_over_func": float(math.exp(float(full_est[2]))),
        "fold_ratio_ci95_low": float(math.exp(float(lo95_vec[2]))),
        "fold_ratio_ci95_high": float(math.exp(float(hi95_vec[2]))),
    }


VALID_SIGNALP6_CLASSES: frozenset[str] = frozenset({"OTHER", "SP", "LIPO", "TAT", "TATLIPO", "PILIN"})

STAGE1_ROUTE_A1_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "candidate_gene",
    "pathway",
    "missing_reaction",
    "reaction_name",
    "ec",
    "review_priority",
    "inherited_catalogue_dark",
    "family_id",
    "contig",
    "start",
    "end",
    "strand",
    "partial_flag",
    "orf_status",
    "coord_aa_length",
    "has_full_seq",
    "coords_valid",
    "contig_in_catalogue",
    "contig_edge_status",
    "candidate_promotion_eligible",
    "n_flank_genes_pm5",
    "n_known_flank_genes_pm5",
    "flank_annotations_pm5",
    "context_score",
    "sequence_bitscore",
    "signalp_class",
    "sp_any",
    "tmbed_S",
    "tm_helices",
    "tm_beta",
    "checkm_completeness",
)

STAGE1_ROUTE_A2_REQUIRED_COLS: tuple[str, ...] = (
    "selection_route",
    "family_id",
    "in_route_a1_gapfill_priority",
    "primary_biological_category",
    "s_layer_and_enzyme_counterexample_audit",
    "n_member_genes_thermo237",
    "n_distinct_genomes",
    "n_distinct_orders",
    "n_independent_genus_blocks",
    "median_aa_length",
    "frac_candidate_promotion_eligible",
    "frac_contig_interior",
    "frac_sp_any",
    "frac_tmbed_S",
    "frac_sp_or_S",
    "mean_tm_helices",
    "mean_tm_beta",
    "rep_locus_gene",
    "rep_locus_n_flank_pm5",
    "rep_locus_n_known_flank_pm5",
    "rep_locus_flank_summary",
    "independent_selection_score",
)


def validate_stage1_to_stage1b_schema(
    a1_props_df: pd.DataFrame,
    a2_fams_df: pd.DataFrame,
) -> None:
    """Validate that Stage 1 output tables contain all required columns consumed by Stage 1b before running the locus audit."""
    missing_a1 = [c for c in STAGE1_ROUTE_A1_REQUIRED_COLS if c not in a1_props_df.columns]
    if missing_a1:
        raise ValueError(f"Stage 1 Route A1 proposals table is missing required columns for Stage 1b: {missing_a1}")
    missing_a2 = [c for c in STAGE1_ROUTE_A2_REQUIRED_COLS if c not in a2_fams_df.columns]
    if missing_a2:
        raise ValueError(f"Stage 1 Route A2 family table is missing required columns for Stage 1b: {missing_a2}")
    if a1_props_df.empty or a2_fams_df.empty:
        raise ValueError("Stage 1 Route A1 or Route A2 table is empty")


def validate_signalp_chunk_against_fasta(
    fasta_path: Path | str,
    pred_gz_path: Path | str,
    *,
    expected_fasta_sha256: str | None = None,
) -> list[tuple[str, str, str, bool, float]]:
    """Validate a SignalP 6.0h prediction output (`prediction_results.txt.gz`) against its input chunk FASTA.

    Enforces `external-signalp-final-checks.txt` Item 2:
    1. Verifies the input chunk FASTA SHA-256 against `expected_fasta_sha256` (when provided).
    2. Extracts all `(genome_id, gene_id)` headers from the input FASTA, rejecting empty or duplicate keys.
    3. Parses `pred_gz_path`, rejecting:
       - malformed rows (< 8 tab-separated columns),
       - missing/invalid headers without `genome_id|gene_id`,
       - duplicate `(genome_id, gene_id)` records (so duplicated replacement rows with unchanged counts fail!),
       - unknown SignalP prediction classes outside `{'OTHER', 'SP', 'LIPO', 'TAT', 'TATLIPO', 'PILIN'}`,
       - non-finite or out-of-range `[0.0, 1.0]` probabilities or rows whose 6 class probabilities do not sum to ~1.0.
    4. Asserts exact set equality `scored_keys == fasta_keys` and exact count equality.
    Returns a list of validated `(genome_id, gene_id, pred_class, has_sp, non_other_prob)` tuples.
    """
    import gzip as _gzip

    f_path = Path(fasta_path)
    p_path = Path(pred_gz_path)
    if not f_path.is_file():
        raise FileNotFoundError(f"Missing SignalP chunk FASTA: {f_path}")
    if not p_path.is_file():
        raise FileNotFoundError(f"Missing SignalP chunk prediction file: {p_path}")

    raw_fasta = f_path.read_bytes()
    actual_sha = hashlib.sha256(raw_fasta).hexdigest()
    if expected_fasta_sha256 is not None and actual_sha != expected_fasta_sha256.strip().lower():
        raise ValueError(
            f"SignalP chunk FASTA SHA-256 mismatch for {f_path}: expected {expected_fasta_sha256}, got {actual_sha}"
        )

    fasta_keys: set[tuple[str, str]] = set()
    for raw_line in raw_fasta.decode("utf-8", errors="replace").splitlines():
        if raw_line.startswith(">"):
            tok = raw_line[1:].strip().split()[0] if raw_line[1:].strip() else ""
            if "|" not in tok:
                raise ValueError(f"Invalid chunk FASTA header (missing 'genome|gene'): {raw_line}")
            gid, gene_id = tok.split("|", 1)
            gid, gene_id = gid.strip(), gene_id.strip()
            if not gid or not gene_id:
                raise ValueError(f"Empty genome_id or gene_id in chunk FASTA header: {raw_line}")
            key = (gid, gene_id)
            if key in fasta_keys:
                raise ValueError(f"Duplicate (genome_id, gene_id) key {key} in chunk FASTA {f_path}")
            fasta_keys.add(key)

    if not fasta_keys:
        raise ValueError(f"Chunk FASTA {f_path} contains 0 sequences")

    scored_keys: set[tuple[str, str]] = set()
    validated_rows: list[tuple[str, str, str, bool, float]] = []

    with _gzip.open(p_path, "rt", encoding="utf-8") as fh:
        for line_no, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                raise ValueError(
                    f"Malformed SignalP prediction row at {p_path}:{line_no} (expected >=8 columns, got {len(parts)})"
                )
            hdr = parts[0].strip()
            pred = parts[1].strip()
            if pred not in VALID_SIGNALP6_CLASSES:
                raise ValueError(
                    f"Invalid SignalP prediction class '{pred}' at {p_path}:{line_no}; expected one of {sorted(VALID_SIGNALP6_CLASSES)}"
                )
            tok = hdr.split()[0] if hdr else ""
            if "|" not in tok:
                raise ValueError(f"Missing 'genome|gene' prefix in SignalP ID '{hdr}' at {p_path}:{line_no}")
            gid, gene_id = tok.split("|", 1)
            gid, gene_id = gid.strip(), gene_id.strip()
            if not gid or not gene_id:
                raise ValueError(f"Empty genome_id or gene_id in SignalP ID '{hdr}' at {p_path}:{line_no}")
            key = (gid, gene_id)
            if key in scored_keys:
                raise ValueError(f"Duplicate SignalP prediction for key {key} at {p_path}:{line_no}")

            probs: list[float] = []
            for p_str in parts[2:8]:
                try:
                    pv = float(p_str)
                except ValueError as exc:
                    raise ValueError(f"Non-numeric SignalP probability '{p_str}' at {p_path}:{line_no}") from exc
                if not math.isfinite(pv) or pv < -1e-6 or pv > 1.000001:
                    raise ValueError(f"Non-finite or out-of-range SignalP probability {pv} at {p_path}:{line_no}")
                probs.append(pv)
            prob_sum = sum(probs)
            if abs(prob_sum - 1.0) > 0.02:
                raise ValueError(
                    f"SignalP class probabilities sum to {prob_sum:.4f} (expected ~1.0) at {p_path}:{line_no}"
                )

            scored_keys.add(key)
            p_other = max(0.0, min(1.0, probs[0]))
            validated_rows.append((gid, gene_id, pred, pred != "OTHER", 1.0 - p_other))

    if scored_keys != fasta_keys:
        missing_in_pred = fasta_keys - scored_keys
        unexpected_in_pred = scored_keys - fasta_keys
        raise ValueError(
            f"SignalP ID coverage mismatch for {p_path} vs {f_path}: "
            f"n_fasta={len(fasta_keys)}, n_scored={len(scored_keys)}, "
            f"n_missing={len(missing_in_pred)} (e.g. {list(sorted(missing_in_pred))[:3]}), "
            f"n_unexpected={len(unexpected_in_pred)} (e.g. {list(sorted(unexpected_in_pred))[:3]})"
        )

    return validated_rows


CANONICAL_REFINED_REPRESENTATIVES_SHA256: str = (
    "dbe9db11362e77c4274543431ad40ee92ed6b181c2c0f7ddc055c4f5a32971d6"
)


def parse_gtdb_taxonomy(tax_str: object) -> dict[str, str]:
    """Parse a semicolon-delimited GTDB taxonomy string into rank keys (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`)."""
    out = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }
    if tax_str is None or pd.isna(tax_str):
        return out
    prefix_map = {
        "d__": "domain",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }
    for tok in str(tax_str).strip().split(";"):
        t = tok.strip()
        for pref, rname in prefix_map.items():
            if t.startswith(pref):
                out[rname] = t
                break
    return out


def join_and_validate_canonical_eq_taxonomy(
    eq_df: pd.DataFrame,
    canonical_tsv_path: Path | str,
    *,
    genome_col: str = "genome",
    legacy_taxonomy_col: str = "gtdb_classification",
    expected_sha256: str | None = CANONICAL_REFINED_REPRESENTATIVES_SHA256,
    expected_n_eq: int | None = 237,
    require_canonical_diff_count: int | None = 42,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Join the 237 EQ Thermoplasmatota genomes to canonical `refined_representatives.tsv` (`mag`) with strict provenance and uniqueness gates.

    Enforces `canonical-taxonomy-sensitivity-design.txt` (Stage 2b section):
    1. Verifies the SHA-256 of `canonical_tsv_path` against `expected_sha256` (when provided).
    2. Verifies primary-key uniqueness and non-nullness on both `eq_df[genome_col]` and `canonical_df['mag']`.
    3. Verifies exact expected EQ genome count (`expected_n_eq`, default 237) and 1:1 membership in `canonical_df['mag']`
       (rejecting same-count wrong-genome substitutions).
    4. Verifies every joined EQ genome has canonical phylum `p__Thermoplasmatota` in `gtdb_r232_classification`
       (rejecting bacterial/other-phylum MAG substitutions even if present in `refined_representatives.tsv`).
    5. Verifies exact completeness and contamination agreement (`atol=1e-3`) when present in `eq_df`.
    6. Extracts both canonical GTDB R232 rank labels (`canonical_class_label`, `canonical_order_label`,
       `canonical_family_label`, `canonical_genus_label`, `canonical_genus_or_singleton_block`) and legacy analytical
       rank labels (`legacy_class_label`, `legacy_order_label`, `legacy_family_label`, `legacy_genus_label`,
       `legacy_genus_or_singleton_block`), explicitly labelling genus-or-singleton blocks as resampling blocks rather
       than independent phylogenetic lineages.
    7. Verifies `n_taxonomy_strings_differ == require_canonical_diff_count` (when provided) so stale taxonomy
       label substitution cannot pass silently as canonical.
    """
    c_path = Path(canonical_tsv_path)
    if not c_path.is_file():
        raise FileNotFoundError(f"Missing canonical refined_representatives.tsv: {c_path}")

    actual_sha = hashlib.sha256(c_path.read_bytes()).hexdigest()
    if expected_sha256 is not None and actual_sha != expected_sha256.strip().lower():
        raise ValueError(
            f"Canonical refined_representatives.tsv SHA-256 mismatch for {c_path}: "
            f"expected {expected_sha256}, got {actual_sha}"
        )

    can_df = pd.read_csv(c_path, sep="\t", low_memory=False)
    for req_c in ("mag", "Completeness", "Contamination", "gtdb_r232_classification"):
        if req_c not in can_df.columns:
            raise ValueError(f"Canonical refined_representatives.tsv missing required column '{req_c}'")

    validate_primary_key(can_df, ["mag"], "canonical_refined_representatives")
    if genome_col not in eq_df.columns:
        raise ValueError(f"eq_df missing genome column '{genome_col}'")
    actual_leg_col = legacy_taxonomy_col
    if actual_leg_col not in eq_df.columns and "taxonomy" in eq_df.columns:
        actual_leg_col = "taxonomy"
    if actual_leg_col not in eq_df.columns:
        raise ValueError(f"eq_df missing legacy taxonomy column '{legacy_taxonomy_col}'")

    validate_primary_key(eq_df, [genome_col], "eq_df")
    if expected_n_eq is not None and len(eq_df) != int(expected_n_eq):
        raise ValueError(f"Expected {expected_n_eq} EQ genomes, got {len(eq_df)}")

    can_sub = can_df[["mag", "Completeness", "Contamination", "gtdb_r232_classification"]].rename(
        columns={
            "mag": genome_col,
            "Completeness": "canonical_completeness",
            "Contamination": "canonical_contamination",
            "gtdb_r232_classification": "canonical_gtdb_r232_classification",
        }
    )
    joined = validated_merge(
        eq_df.copy(),
        can_sub,
        on=[genome_col],
        how="left",
        validate="one_to_one",
        left_name="eq_df",
        right_name="canonical_refined_representatives",
        require_all_left=True,
    )

    comp_col = "completeness" if "completeness" in joined.columns else ("checkm_completeness" if "checkm_completeness" in joined.columns else None)
    if comp_col is not None:
        c_diff = np.abs(
            pd.to_numeric(joined[comp_col], errors="raise").to_numpy(dtype=float)
            - pd.to_numeric(joined["canonical_completeness"], errors="raise").to_numpy(dtype=float)
        )
        if np.any(c_diff > 1e-3):
            raise ValueError(f"Completeness mismatch between eq_df and canonical refined_representatives (max diff={c_diff.max()})")
    if "contamination" in joined.columns:
        ct_diff = np.abs(
            pd.to_numeric(joined["contamination"], errors="raise").to_numpy(dtype=float)
            - pd.to_numeric(joined["canonical_contamination"], errors="raise").to_numpy(dtype=float)
        )
        if np.any(ct_diff > 1e-3):
            raise ValueError(f"Contamination mismatch between eq_df and canonical refined_representatives (max diff={ct_diff.max()})")

    joined["legacy_analytical_taxonomy"] = joined[actual_leg_col].astype(str).str.strip()
    joined["canonical_gtdb_r232_classification"] = joined["canonical_gtdb_r232_classification"].astype(str).str.strip()

    can_phyla = []
    can_classes, can_orders, can_families, can_genera, can_blocks = [], [], [], [], []
    leg_classes, leg_orders, leg_families, leg_genera, leg_blocks = [], [], [], [], []

    for row in joined.itertuples(index=False):
        gid = str(getattr(row, genome_col)).strip()
        c_tax = str(row.canonical_gtdb_r232_classification)
        l_tax = str(row.legacy_analytical_taxonomy)
        c_ranks = parse_gtdb_taxonomy(c_tax)
        l_ranks = parse_gtdb_taxonomy(l_tax)

        if c_ranks["phylum"] != "p__Thermoplasmatota":
            raise ValueError(
                f"Joined genome '{gid}' has non-Thermoplasmatota canonical phylum '{c_ranks['phylum']}' "
                f"(possible wrong-genome substitution)"
            )
        can_phyla.append(c_ranks["phylum"])
        can_classes.append(c_ranks["class"])
        can_orders.append(c_ranks["order"])
        can_families.append(c_ranks["family"])
        can_genera.append(c_ranks["genus"])
        can_blocks.append(
            f"UNCLASSIFIED_GENUS_SINGLETON::{gid}" if c_ranks["genus"] == "g__" else f"GENUS::{c_ranks['genus']}"
        )

        leg_classes.append(l_ranks["class"])
        leg_orders.append(l_ranks["order"])
        leg_families.append(l_ranks["family"])
        leg_genera.append(l_ranks["genus"])
        leg_blocks.append(
            f"UNCLASSIFIED_GENUS_SINGLETON::{gid}" if l_ranks["genus"] == "g__" else f"GENUS::{l_ranks['genus']}"
        )

    joined["canonical_class_label"] = can_classes
    joined["canonical_order_label"] = can_orders
    joined["canonical_family_label"] = can_families
    joined["canonical_genus_label"] = can_genera
    joined["canonical_genus_or_singleton_block"] = can_blocks

    joined["legacy_class_label"] = leg_classes
    joined["legacy_order_label"] = leg_orders
    joined["legacy_family_label"] = leg_families
    joined["legacy_genus_label"] = leg_genera
    joined["legacy_genus_or_singleton_block"] = leg_blocks

    joined["taxonomy_string_differs"] = joined["canonical_gtdb_r232_classification"] != joined["legacy_analytical_taxonomy"]
    joined["class_label_differs"] = joined["canonical_class_label"] != joined["legacy_class_label"]
    joined["order_label_differs"] = joined["canonical_order_label"] != joined["legacy_order_label"]
    joined["family_label_differs"] = joined["canonical_family_label"] != joined["legacy_family_label"]
    joined["genus_label_differs"] = joined["canonical_genus_label"] != joined["legacy_genus_label"]

    n_tax_diff = int(joined["taxonomy_string_differs"].sum())
    if require_canonical_diff_count is not None and n_tax_diff != int(require_canonical_diff_count):
        raise ValueError(
            f"Expected {require_canonical_diff_count} EQ taxonomy strings to differ between canonical and legacy, "
            f"got {n_tax_diff} (possible stale taxonomy label substitution)"
        )

    summary = {
        "canonical_tsv_path": str(c_path),
        "canonical_tsv_sha256": actual_sha,
        "n_eq_genomes": int(len(joined)),
        "n_taxonomy_strings_differ": n_tax_diff,
        "n_class_labels_differ": int(joined["class_label_differs"].sum()),
        "n_order_labels_differ": int(joined["order_label_differs"].sum()),
        "n_family_labels_differ": int(joined["family_label_differs"].sum()),
        "n_genus_labels_differ": int(joined["genus_label_differs"].sum()),
        "canonical_n_named_orders": int(len({o for o in can_orders if o != "o__"})),
        "canonical_n_unassigned_order_genomes": int(sum(1 for o in can_orders if o == "o__")),
        "canonical_n_named_genera": int(len({g for g in can_genera if g != "g__"})),
        "canonical_n_unassigned_genus_genomes": int(sum(1 for g in can_genera if g == "g__")),
        "canonical_n_genus_or_singleton_resampling_blocks": int(len(set(can_blocks))),
        "legacy_n_named_orders": int(len({o for o in leg_orders if o != "o__"})),
        "legacy_n_unassigned_order_genomes": int(sum(1 for o in leg_orders if o == "o__")),
        "legacy_n_named_genera": int(len({g for g in leg_genera if g != "g__"})),
        "legacy_n_unassigned_genus_genomes": int(sum(1 for g in leg_genera if g == "g__")),
        "legacy_n_genus_or_singleton_resampling_blocks": int(len(set(leg_blocks))),
        "block_semantics_note": "genus_or_unclassified_singleton_resampling_blocks_not_independent_phylogenetic_lineages",
    }
    return joined, summary


def compute_thermo_proteomes_tar_sha256_map(
    tar_gz_path: Path | str,
    target_genome_ids: set[str] | None = None,
) -> dict[str, str]:
    """Compute exact SHA-256 digests of `.faa` proteome members inside `thermo_proteomes.tar.gz`."""
    import tarfile as _tarfile

    t_path = require_file(tar_gz_path, "thermo_proteomes.tar.gz")
    sha_map: dict[str, str] = {}
    with _tarfile.open(t_path, "r:gz") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".faa"):
                continue
            gid = Path(member.name).name[:-4]
            if target_genome_ids is not None and gid not in target_genome_ids:
                continue
            if gid in sha_map:
                raise ValueError(f"Duplicate .faa member for genome '{gid}' in {t_path}")
            fobj = tf.extractfile(member)
            if fobj is None:
                raise ValueError(f"Could not extract .faa member for genome '{gid}' in {t_path}")
            sha_map[gid] = hashlib.sha256(fobj.read()).hexdigest()
    if target_genome_ids is not None and set(sha_map.keys()) != set(target_genome_ids):
        missing = sorted(set(target_genome_ids) - set(sha_map.keys()))
        raise ValueError(
            f"Missing {len(missing)} target genome .faa members in {t_path} (e.g. {missing[:5]})"
        )
    return sha_map


def validate_full_stage2b_signalp_and_uniref50_provenance(
    df_484: pd.DataFrame,
    ext_sig_df: pd.DataFrame,
    ext_manifest_df: pd.DataFrame,
    eq_sp_raw: dict[tuple[str, str], tuple[str, float]],
    *,
    expected_ext_genome_faa_sha256: dict[str, str],
) -> pd.DataFrame:
    """Validate and join full Stage 2b SignalP 6.0h + direct UniRef50 inputs (`external-join-final-launch.txt` Item 3).

    Enforces:
    1. Primary-key uniqueness on `df_484` `(genome, gene)`, `ext_sig_df` `(genome_id, gene_id)`, and
       `ext_manifest_df` `(genome_id)` (rejecting duplicate SignalP keys or silent dict overwrites).
    2. Valid SignalP 6.0h prediction classes (`VALID_SIGNALP6_CLASSES`) across all external and EQ rows.
    3. Exact 247-genome set equality and per-genome protein count equality across `df_484` external rows,
       `ext_sig_df`, and `ext_manifest_df`.
    4. Exact per-genome `.faa` SHA-256 verification of `ext_manifest_df['proteome_faa_sha256']` against
       `expected_ext_genome_faa_sha256` (from the direct-search source `thermo_proteomes.tar.gz`), failing
       immediately if a file exists with wrong provenance.
    5. Exact `(genome, gene)` key-set equality between `df_484` external rows and `ext_sig_df`, and exact
       `gene` coverage + conflict-free prediction lookup for all `df_484` EQ rows.
    """
    for req_c in ("genome", "gene", "cohort"):
        if req_c not in df_484.columns:
            raise ValueError(f"df_484 missing required column '{req_c}'")
    validate_primary_key(df_484, ["genome", "gene"], "df_484_direct_uniref50")

    for req_c in ("genome_id", "gene_id", "signalp_prediction"):
        if req_c not in ext_sig_df.columns:
            raise ValueError(f"ext_sig_df missing required column '{req_c}'")
    validate_primary_key(ext_sig_df, ["genome_id", "gene_id"], "external247_signalp6_predictions")

    for req_c in ("genome_id", "n_proteins", "proteome_faa_sha256"):
        if req_c not in ext_manifest_df.columns:
            raise ValueError(f"ext_manifest_df missing required column '{req_c}'")
    validate_primary_key(ext_manifest_df, ["genome_id"], "external_signalp_manifest")

    df_eq = df_484[df_484["cohort"] == "EQ_Thermoplasmatota_237"]
    df_ext = df_484[df_484["cohort"] == "External_Comparators_Deduplicated_247"]
    if len(df_eq) + len(df_ext) != len(df_484) or len(df_eq) == 0 or len(df_ext) == 0:
        raise ValueError("Unexpected or empty cohort partition in df_484")

    ext_genomes_484 = set(df_ext["genome"].astype(str).unique())
    ext_genomes_man = set(ext_manifest_df["genome_id"].astype(str).unique())
    if ext_genomes_484 != ext_genomes_man:
        raise ValueError(
            f"External genome set mismatch between df_484 ({len(ext_genomes_484)}) and "
            f"external_signalp_manifest ({len(ext_genomes_man)})"
        )
    if set(expected_ext_genome_faa_sha256.keys()) != ext_genomes_man:
        raise ValueError(
            "expected_ext_genome_faa_sha256 keys do not match external_signalp_manifest genome_ids"
        )

    # Verify per-genome SHA-256 and protein counts in external_signalp_manifest against direct-search source
    man_counts: dict[str, int] = {}
    for row in ext_manifest_df.itertuples(index=False):
        gid = str(row.genome_id).strip()
        man_sha = str(row.proteome_faa_sha256).strip().lower()
        exp_sha = str(expected_ext_genome_faa_sha256[gid]).strip().lower()
        if not man_sha or len(man_sha) != 64 or man_sha != exp_sha:
            raise ValueError(
                f"Provenance SHA-256 mismatch for external genome '{gid}': "
                f"manifest={man_sha}, direct_search_source={exp_sha}"
            )
        man_counts[gid] = int(row.n_proteins)

    ext_484_counts = df_ext.groupby("genome", sort=False).size().to_dict()
    for gid, n_exp in man_counts.items():
        if int(ext_484_counts.get(gid, -1)) != n_exp:
            raise ValueError(
                f"Per-genome protein count mismatch for '{gid}': manifest={n_exp}, df_484={ext_484_counts.get(gid)}"
            )

    # Build ext_sp_map without silent overwrite and verify valid SignalP classes
    ext_sp_map: dict[tuple[str, str], tuple[str, bool]] = {}
    for row in ext_sig_df.itertuples(index=False):
        key = (str(row.genome_id).strip(), str(row.gene_id).strip())
        if key in ext_sp_map:
            raise ValueError(f"Duplicate SignalP key {key} in ext_sig_df")
        pred = str(row.signalp_prediction).strip()
        if pred not in VALID_SIGNALP6_CLASSES:
            raise ValueError(f"Invalid SignalP prediction class '{pred}' for {key}")
        ext_sp_map[key] = (pred, pred in {"SP", "LIPO", "TAT", "TATLIPO"})

    ext_484_keys = set(zip(df_ext["genome"].astype(str), df_ext["gene"].astype(str)))
    if set(ext_sp_map.keys()) != ext_484_keys:
        missing_k = sorted(ext_484_keys - set(ext_sp_map.keys()))
        extra_k = sorted(set(ext_sp_map.keys()) - ext_484_keys)
        raise ValueError(
            f"External SignalP (genome, gene) key set mismatch vs df_484: "
            f"n_missing={len(missing_k)} (e.g. {missing_k[:3]}), n_extra={len(extra_k)} (e.g. {extra_k[:3]})"
        )

    # Build eq_sp_map checking for conflicting predictions across source tags
    eq_sp_map: dict[str, tuple[str, bool]] = {}
    for (_, gene_id), (pred, _) in eq_sp_raw.items():
        g_key = str(gene_id).strip()
        p_clean = str(pred).strip()
        if p_clean not in VALID_SIGNALP6_CLASSES:
            raise ValueError(f"Invalid EQ SignalP prediction class '{p_clean}' for gene '{g_key}'")
        if g_key in eq_sp_map and eq_sp_map[g_key][0] != p_clean:
            raise ValueError(
                f"Conflicting EQ SignalP predictions for gene '{g_key}': '{eq_sp_map[g_key][0]}' vs '{p_clean}'"
            )
        eq_sp_map[g_key] = (p_clean, p_clean in {"SP", "LIPO", "TAT", "TATLIPO"})

    eq_484_genes = set(df_eq["gene"].astype(str))
    missing_eq_genes = sorted(eq_484_genes - set(eq_sp_map.keys()))
    if missing_eq_genes:
        raise ValueError(
            f"Missing EQ SignalP predictions for {len(missing_eq_genes)} genes in df_484 (e.g. {missing_eq_genes[:5]})"
        )

    sp_pred_col: list[str] = []
    sp_any_col: list[bool] = []
    for g_id, gene_id, coh in zip(df_484["genome"].astype(str), df_484["gene"].astype(str), df_484["cohort"].astype(str)):
        if coh == "EQ_Thermoplasmatota_237":
            pred, is_sp = eq_sp_map[gene_id]
        else:
            pred, is_sp = ext_sp_map[(g_id, gene_id)]
        sp_pred_col.append(pred)
        sp_any_col.append(bool(is_sp))

    out = df_484.copy()
    out["has_signalp6"] = True
    out["signalp6_prediction"] = sp_pred_col
    out["sp_any"] = sp_any_col
    return out



