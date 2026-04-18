"""Compositional co-occurrence networks for the Empty Quarter dataset.

Uses CLR-transformed genus abundances + Spearman rank correlation
(the pragmatic proxy for SparCC / SPIEC-EASI when full-compositional
inference is too expensive). Multiple-testing controlled via
Benjamini-Hochberg. Community detection with Louvain modularity.
Keystones identified by a composite centrality score (degree × BC
× closeness, normalised).

For ≤ 500 genera and ≤ 500 samples per compartment this runs in a few
seconds and is reproducible from the cached feature table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def compositional_correlation(
    gen_count: pd.DataFrame,
    *,
    min_prevalence: float = 0.20,
    presence_ra: float = 0.001,
    pseudo: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """CLR-transform a (genera × samples) count table and return the
    Spearman correlation matrix plus pairwise p-values.

    Genera not passing ``min_prevalence`` (fraction of samples where
    relative abundance ≥ ``presence_ra``) are dropped before correlation.
    """
    rel = gen_count.div(gen_count.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    keep = rel.index[(rel >= presence_ra).mean(axis=1) >= min_prevalence]
    sub = gen_count.loc[keep]
    # CLR with pseudo-count
    X = np.log(sub.to_numpy() + pseudo)
    X = X - X.mean(axis=0, keepdims=True)
    # Spearman on CLR values
    rho, p = spearmanr(X, axis=1)
    rho = pd.DataFrame(rho, index=keep, columns=keep)
    p = pd.DataFrame(p, index=keep, columns=keep)
    return rho, p


def build_network(
    rho: pd.DataFrame,
    p: pd.DataFrame,
    *,
    rho_threshold: float = 0.4,
    q_threshold: float = 0.01,
) -> nx.Graph:
    """Threshold the correlation matrix and build an undirected network.

    Only the upper triangle is used; self-loops are excluded. Edge weight
    is |rho|; sign of correlation is stored as ``sign``.
    """
    iu = np.triu_indices(rho.shape[0], k=1)
    r = rho.to_numpy()[iu]
    pv = p.to_numpy()[iu]
    # BH across all pairs
    _, q, _, _ = multipletests(pv, method="fdr_bh")
    names = rho.index.tolist()
    idx_i, idx_j = iu
    G = nx.Graph()
    G.add_nodes_from(names)
    for (i, j, rij, qij) in zip(idx_i, idx_j, r, q):
        if abs(rij) >= rho_threshold and qij <= q_threshold:
            G.add_edge(names[i], names[j],
                       weight=abs(rij), rho=float(rij),
                       sign=int(np.sign(rij)), q=float(qij))
    return G


def louvain_modules(G: nx.Graph, seed: int = 1) -> dict[str, int]:
    """Louvain modularity — only considers positive-sign edges."""
    G_pos = nx.Graph()
    G_pos.add_nodes_from(G.nodes)
    for u, v, d in G.edges(data=True):
        if d["sign"] > 0:
            G_pos.add_edge(u, v, weight=d["weight"])
    if G_pos.number_of_edges() == 0:
        return {n: 0 for n in G.nodes}
    try:
        partition = nx.community.louvain_communities(G_pos, seed=seed)
    except AttributeError:
        partition = nx.algorithms.community.louvain_communities(G_pos, seed=seed)
    mod = {}
    for i, nodes in enumerate(partition):
        for n in nodes:
            mod[n] = i
    return mod


def keystone_score(G: nx.Graph) -> pd.DataFrame:
    """Composite keystone score per node.

    ``keystone = 0.5·degree_norm + 0.3·betweenness_norm + 0.2·closeness_norm``
    where each component is scaled to the [0, 1] range across nodes.
    """
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["node", "degree", "betweenness", "closeness",
                                     "keystone"])

    deg = dict(G.degree())
    if G.number_of_edges() == 0:
        return pd.DataFrame({
            "node": list(G.nodes), "degree": [0]*G.number_of_nodes(),
            "betweenness": [0]*G.number_of_nodes(),
            "closeness": [0]*G.number_of_nodes(),
            "keystone": [0.0]*G.number_of_nodes(),
        })

    bc = nx.betweenness_centrality(G, normalized=True)
    cc = nx.closeness_centrality(G)

    def _rel(d: dict) -> dict:
        vals = np.array(list(d.values()), dtype=float)
        mx = vals.max() if vals.max() > 0 else 1.0
        return {k: v / mx for k, v in d.items()}
    deg_r = _rel(deg); bc_r = _rel(bc); cc_r = _rel(cc)
    ks = {n: 0.5 * deg_r[n] + 0.3 * bc_r[n] + 0.2 * cc_r[n] for n in G.nodes}
    out = pd.DataFrame({
        "node": list(G.nodes),
        "degree": [deg[n] for n in G.nodes],
        "betweenness": [bc[n] for n in G.nodes],
        "closeness": [cc[n] for n in G.nodes],
        "keystone": [ks[n] for n in G.nodes],
    }).sort_values("keystone", ascending=False)
    return out
