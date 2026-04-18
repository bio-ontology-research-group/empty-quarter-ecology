"""Diversity and composition helpers that operate on the cached
feature-table / taxonomy / metadata produced by :mod:`eq.loader`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def relative_abundance(ft: pd.DataFrame) -> pd.DataFrame:
    """Column-wise relative abundance (ASVs x samples)."""
    col_sum = ft.sum(axis=0)
    col_sum = col_sum.replace(0, np.nan)
    return ft.div(col_sum, axis=1).fillna(0.0)


def aggregate_to_rank(
    ft: pd.DataFrame, tax: pd.DataFrame, rank: str = "phylum"
) -> pd.DataFrame:
    """Sum ASVs per taxonomic rank. Returns (taxa x samples) DataFrame."""
    if rank not in tax.columns:
        raise ValueError(f"rank {rank!r} not in taxonomy columns: {list(tax.columns)}")
    asv_to_rank = tax[rank].fillna("unclassified").replace("", "unclassified")
    joined = ft.join(asv_to_rank.rename("rank_label"))
    grouped = joined.groupby("rank_label").sum(numeric_only=True)
    return grouped


def top_n_groups(
    grouped: pd.DataFrame, n: int = 10, group_col: str = "Other"
) -> pd.DataFrame:
    """Keep the ``n`` most abundant rank labels, collapse the rest into ``group_col``."""
    totals = grouped.sum(axis=1).sort_values(ascending=False)
    top = totals.head(n).index.tolist()
    other = grouped.loc[~grouped.index.isin(top)].sum(axis=0)
    out = grouped.loc[top].copy()
    out.loc[group_col] = other
    return out


def shannon(ft: pd.DataFrame) -> pd.Series:
    """Shannon diversity (base e) per sample (column)."""
    rel = relative_abundance(ft)
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(rel.where(rel > 0))
    h = -(rel * logs.fillna(0)).sum(axis=0)
    return pd.Series(h.values, index=ft.columns, name="shannon")


def richness(ft: pd.DataFrame) -> pd.Series:
    """Observed ASV count (>0 reads) per sample."""
    return (ft > 0).sum(axis=0).rename("observed_asvs")


def simpson(ft: pd.DataFrame) -> pd.Series:
    """1 - sum(p_i^2) Simpson index per sample."""
    rel = relative_abundance(ft)
    d = 1 - (rel ** 2).sum(axis=0)
    return d.rename("simpson")
