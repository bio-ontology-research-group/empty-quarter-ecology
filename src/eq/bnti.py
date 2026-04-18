"""Vectorised β-NTI / Raup-Crick (Stegen et al. 2012, 2013).

The bottleneck in picante::comdistnt is the abundance-weighted mean
nearest-neighbour phylogenetic distance computation, which is naturally
expressed as a matrix product once we precompute per-sample nearest-
distance vectors. This implementation runs the observed step in seconds
and the nulls in parallel via multiprocessing.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd
from ete3 import Tree

log = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Tree / cophenetic helpers
# ----------------------------------------------------------------------
def cophenetic_matrix(tree_nwk: str, asv_order: list[str]) -> np.ndarray:
    """Cophenetic (patristic) distance matrix for the given tip set.

    Prunes the tree to ``asv_order``, recomputes root-to-tip distances,
    and fills a symmetric matrix. O(n^2) memory; feasible for ≤ 5k tips.
    """
    t = Tree(tree_nwk, format=1)
    tips_in_tree = set(t.get_leaf_names())
    keep = [a for a in asv_order if a in tips_in_tree]
    t.prune(keep, preserve_branch_length=True)
    name2leaf = {leaf.name: leaf for leaf in t.get_leaves()}
    n = len(asv_order)
    D = np.full((n, n), np.nan, dtype=np.float64)
    # We compute pairwise via the LCA of each pair; efficient enough
    for i, ai in enumerate(asv_order):
        li = name2leaf.get(ai)
        if li is None:
            continue
        for j in range(i, n):
            aj = asv_order[j]
            lj = name2leaf.get(aj)
            if lj is None:
                continue
            if i == j:
                D[i, j] = 0.0
            else:
                D[i, j] = D[j, i] = t.get_distance(li, lj)
    return D


def cophenetic_matrix_ape(tree_nwk: str, asv_order: list[str]) -> np.ndarray:
    """Fast cophenetic matrix via a simple DFS accumulator.

    Uses only the dendropy-free ete3 tree traversal. For ≤ 5k tips the
    resulting (n × n) matrix is ≤ 100 MB.
    """
    t = Tree(tree_nwk, format=1)
    tips_in_tree = {leaf.name for leaf in t.get_leaves()}
    wanted = [a for a in asv_order if a in tips_in_tree]
    t.prune(wanted, preserve_branch_length=True)

    # Index each tip to its row in D
    idx = {a: i for i, a in enumerate(asv_order)}
    n = len(asv_order)
    D = np.full((n, n), np.nan, dtype=np.float64)

    # For each internal node, the tips descending from its left and right
    # branches are separated through that node; distance = sum of branch
    # lengths to that LCA from both sides. We accumulate by pre-computing
    # for each internal node the list of (tip, distance-to-this-node)
    # pairs for each child subtree, then combining across pairs of children.
    for node in t.traverse("postorder"):
        if node.is_leaf():
            node._subtree_tips = [(node.name, 0.0)]
        else:
            # Combine children: pairwise tips across different children
            child_sets = [ch._subtree_tips for ch in node.children]
            for a_idx, ca in enumerate(child_sets):
                for cb in child_sets[a_idx + 1:]:
                    for (ta, da) in ca:
                        if ta not in idx: continue
                        i = idx[ta]
                        for (tb, db) in cb:
                            if tb not in idx: continue
                            j = idx[tb]
                            d = da + db + (0.0 if ta == tb else 0.0)
                            D[i, j] = D[j, i] = d
            # Merge subtree lists with branch length to parent
            combined = []
            for ch in node.children:
                bl = float(ch.dist)
                combined.extend((t, d + bl) for (t, d) in ch._subtree_tips)
                del ch._subtree_tips
            node._subtree_tips = combined
    np.fill_diagonal(D, 0.0)
    # Any nan cells correspond to ASVs not in the tree: leave nan, caller handles
    return D


# ----------------------------------------------------------------------
# β-MNTD, vectorised
# ----------------------------------------------------------------------
def beta_mntd_weighted(
    rel: np.ndarray, phy_d: np.ndarray, large_value: float = 1e12
) -> np.ndarray:
    """Abundance-weighted pairwise β-MNTD.

    Parameters
    ----------
    rel : (n_samples, n_asvs), non-negative, rows sum to ~1
    phy_d : (n_asvs, n_asvs), symmetric, zero diagonal

    Returns
    -------
    (n_samples, n_samples) symmetric β-MNTD.
    """
    n_s, n_a = rel.shape
    # For each sample j, min_to_j[k] = min over y ∈ present(j) of phy_d[k, y]
    # Use a large_value sentinel for absent y so np.min works row-wise.
    min_to_j = np.empty((n_s, n_a), dtype=np.float64)
    for j in range(n_s):
        mask = rel[j] > 0
        if not mask.any():
            min_to_j[j] = np.nan
            continue
        min_to_j[j] = phy_d[:, mask].min(axis=1)

    # MNTD(i→j) = sum_k rel[i,k] * min_to_j[j, k]
    # Symmetrised: (MNTD(i→j) + MNTD(j→i)) / 2
    mntd_ij = rel @ min_to_j.T
    bmntd = 0.5 * (mntd_ij + mntd_ij.T)
    np.fill_diagonal(bmntd, 0.0)
    return bmntd


# Workers use module-level globals set at pool init time to avoid
# re-pickling large arrays on every call.
_WORKER_REL: np.ndarray | None = None
_WORKER_PHY: np.ndarray | None = None


def _pool_init(rel: np.ndarray, phy_d: np.ndarray) -> None:
    global _WORKER_REL, _WORKER_PHY
    _WORKER_REL = rel
    _WORKER_PHY = phy_d


def _null_one(seed: int) -> np.ndarray:
    assert _WORKER_REL is not None and _WORKER_PHY is not None
    rng = np.random.default_rng(int(seed))
    perm = rng.permutation(_WORKER_PHY.shape[0])
    phy_d_shuf = _WORKER_PHY[np.ix_(perm, perm)]
    return beta_mntd_weighted(_WORKER_REL, phy_d_shuf)


def beta_nti(
    rel: np.ndarray, phy_d: np.ndarray,
    n_null: int = 999, seed: int = 1, n_cores: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute β-NTI using the Stegen framework.

    Returns (observed β-MNTD, null mean, β-NTI).
    """
    n_cores = n_cores or min(os.cpu_count() or 1, 16)
    log.info("observed β-MNTD …")
    obs = beta_mntd_weighted(rel, phy_d)
    log.info("null β-MNTD (%d permutations, %d cores) …", n_null, n_cores)

    seeds = np.arange(seed, seed + n_null)
    n_s = rel.shape[0]
    sm = np.zeros((n_s, n_s), dtype=np.float64)
    sq = np.zeros((n_s, n_s), dtype=np.float64)
    count = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(n_cores, initializer=_pool_init, initargs=(rel, phy_d)) as pool:
        for k, mat in enumerate(
            pool.imap_unordered(
                _null_one, seeds, chunksize=max(1, n_null // (n_cores * 4))
            ), 1
        ):
            sm += mat
            sq += mat * mat
            count += 1
            if k % 100 == 0:
                log.info("  null %d / %d", k, n_null)
    null_mean = sm / count
    null_var = np.maximum(0.0, sq / count - null_mean ** 2)
    null_sd = np.sqrt(null_var)
    bnti = (obs - null_mean) / np.where(null_sd > 0, null_sd, np.nan)
    return obs, null_mean, bnti
