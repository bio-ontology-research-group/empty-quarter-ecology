"""Core / transient / opportunist classification following the
Shade & Handelsman (2012, *Environ Microbiol*) occupancy-abundance
framework and the temporal categories of Dini-Andreote et al. (2015,
*ISME J*). Replaces the naïve "present in ≥4 of 5 trips" rule which
conflates detection with ecology and collapses to the trivial result
that abundant taxa are always core.

Definitions used here:

- **presence** in a sample  := relative abundance ≥ ``presence_ra``
  (default 0.001, i.e. 0.1 %); conservative vs. naïve "≥1 read"
- **global prevalence**     := fraction of all samples with presence
- **per-trip prevalence**   := fraction of that trip's samples with presence
- **temporal core**         := per-trip prevalence ≥ ``core_prev`` (0.5) in
  ``min_trips`` (≥ 4) of the 5 trips
- **persistent opportunist**:= per-trip prevalence ≥ ``opp_prev`` (0.3) in
  1–3 trips and < ``opp_prev`` in the rest
- **rare transient**        := global prevalence < ``rare_prev`` (0.05)
- **conditionally rare**    := global prevalence < 0.2 but max relative
  abundance across samples ≥ 0.01 at least once (Shade 2014)

Numbers with these thresholds are far more ecologically interpretable
than the ≥1-read-per-trip rule — most of the ~1,600 genera fall into
meaningful intermediate categories instead of collapsing into core or
rare.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def classify_genera(
    gen: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    presence_ra: float = 0.001,
    core_prev: float = 0.5,
    opp_prev: float = 0.3,
    rare_prev: float = 0.05,
    min_trips_for_core: int = 4,
    cr_max_ra: float = 0.01,
    cr_prev: float = 0.2,
) -> pd.DataFrame:
    """Classify genera (rows of ``gen``) by temporal + occupancy pattern.

    Parameters
    ----------
    gen : DataFrame of shape (n_genera, n_samples), **raw read counts**
        (we convert to relative abundance internally).
    meta : per-sample metadata with a ``trip`` column indexed by sample.
    presence_ra : minimum relative abundance to count as "present".

    Returns
    -------
    DataFrame indexed by genus with columns:
        * ``global_prev``      fraction of all samples present
        * ``mean_ra``          mean relative abundance across samples
        * ``max_ra``           per-sample maximum relative abundance
        * ``trips_with_core_prev`` count of trips with per-trip prev ≥ core_prev
        * ``category``         one of
          {"core", "opportunist", "conditionally_rare", "rare_transient", "other"}
    """
    assert (gen >= 0).all().all(), "expected non-negative counts"
    col_sum = gen.sum(axis=0).replace(0, np.nan)
    rel = gen.div(col_sum, axis=1).fillna(0.0)
    present = rel >= presence_ra

    # Global prevalence
    global_prev = present.mean(axis=1)

    # Per-trip prevalence
    trip_col = meta.loc[gen.columns, "trip"].astype(int)
    trip_prev = pd.DataFrame(
        {t: present.loc[:, trip_col == t].mean(axis=1).fillna(0.0)
         for t in sorted(trip_col.unique())}
    )
    trips_with_core_prev = (trip_prev >= core_prev).sum(axis=1)
    trips_with_opp_prev = (trip_prev >= opp_prev).sum(axis=1)

    mean_ra = rel.mean(axis=1)
    max_ra = rel.max(axis=1)

    def _cat(i: int) -> str:
        gp = global_prev.iloc[i]
        twc = trips_with_core_prev.iloc[i]
        two = trips_with_opp_prev.iloc[i]
        if twc >= min_trips_for_core:
            return "core"
        if two >= 1 and twc < min_trips_for_core:
            return "opportunist"
        if max_ra.iloc[i] >= cr_max_ra and gp < cr_prev:
            return "conditionally_rare"
        if gp < rare_prev:
            return "rare_transient"
        return "other"

    category = pd.Series(
        [_cat(i) for i in range(gen.shape[0])],
        index=gen.index, name="category",
    )

    out = pd.DataFrame({
        "global_prev": global_prev,
        "mean_ra": mean_ra,
        "max_ra": max_ra,
        "trips_with_core_prev": trips_with_core_prev,
        "trips_with_opp_prev": trips_with_opp_prev,
        "category": category,
    })
    return out
