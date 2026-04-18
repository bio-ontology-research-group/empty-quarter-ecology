#!/usr/bin/env python
"""Run β-NTI for all three compartments and save outputs under cache/bnti/."""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from eq import CACHE_DIR  # noqa: E402
from eq.bnti import beta_nti, cophenetic_matrix_ape  # noqa: E402

N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 999
N_TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 500
N_CORES = int(sys.argv[3]) if len(sys.argv) > 3 else 8

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("bnti")

OUT = CACHE_DIR / "bnti"
OUT.mkdir(exist_ok=True, parents=True)

log.info("loading caches")
ft = pd.read_parquet(CACHE_DIR / "feature_table.parquet")
meta = pd.read_parquet(CACHE_DIR / "metadata.parquet").dropna(subset=["compartment"])
ft = ft[meta.index]
tree_nwk = str(CACHE_DIR / "trees" / "asv_tree_fasttree.nwk")

for comp in ("surface", "deep", "rhizosphere"):
    samps = meta.index[meta["compartment"] == comp]
    log.info("=== %s: %d samples ===", comp, len(samps))
    ft_c = ft[samps]
    prev = (ft_c > 0).sum(axis=1)
    top = prev.sort_values(ascending=False).head(N_TOP).index.tolist()
    ft_c = ft_c.loc[top]

    t0 = time.time()
    phy_d = cophenetic_matrix_ape(tree_nwk, top)
    log.info("  cophenetic matrix %s in %.1fs", phy_d.shape, time.time() - t0)

    rel = (ft_c / ft_c.sum(axis=0)).T.to_numpy()

    t0 = time.time()
    obs, null_mean, bnti = beta_nti(rel, phy_d, n_null=N_NULL, n_cores=N_CORES)
    log.info("  β-NTI in %.1fs", time.time() - t0)

    sample_ids = samps.tolist()
    pd.DataFrame(obs, index=sample_ids, columns=sample_ids).to_parquet(
        OUT / f"obs_bMNTD_{comp}.parquet"
    )
    pd.DataFrame(bnti, index=sample_ids, columns=sample_ids).to_parquet(
        OUT / f"bNTI_{comp}.parquet"
    )
    pd.DataFrame(null_mean, index=sample_ids, columns=sample_ids).to_parquet(
        OUT / f"null_mean_{comp}.parquet"
    )

    v = bnti[np.triu_indices_from(bnti, k=1)]
    v = v[np.isfinite(v)]
    summary = {
        "compartment": comp,
        "n_samples": int(len(sample_ids)),
        "n_asvs": int(len(top)),
        "n_pairs": int(len(v)),
        "mean_bnti": float(np.mean(v)),
        "sd_bnti": float(np.std(v)),
        "frac_abs_gt2": float(np.mean(np.abs(v) > 2)),
        "frac_variable_selection": float(np.mean(v > 2)),
        "frac_homogeneous_selection": float(np.mean(v < -2)),
        "frac_stochastic": float(np.mean(np.abs(v) < 2)),
    }
    log.info("  %s", summary)
    pd.Series(summary).to_json(OUT / f"summary_{comp}.json")

log.info("done")
