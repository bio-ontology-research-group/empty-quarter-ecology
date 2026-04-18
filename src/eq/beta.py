"""Beta-diversity utilities: compositional transforms, distances, PCoA."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA


def clr(ft: pd.DataFrame, pseudo: float = 0.5) -> pd.DataFrame:
    """Centered Log-Ratio on counts. Rows are ASVs, columns are samples.

    Returns (ASV × sample) DataFrame where each column is CLR-transformed.
    """
    X = ft.to_numpy(dtype=float) + pseudo
    logX = np.log(X)
    out = logX - logX.mean(axis=0, keepdims=True)
    return pd.DataFrame(out, index=ft.index, columns=ft.columns)


def aitchison_distance(clr_df: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Euclidean distance between CLR-transformed samples.

    Equivalent to Aitchison distance on the original relative abundances.
    """
    from scipy.spatial.distance import pdist

    X = clr_df.T.to_numpy()  # samples x features
    D = pdist(X, metric="euclidean")
    D = squareform(D)
    return pd.DataFrame(D, index=clr_df.columns, columns=clr_df.columns)


def bray_curtis(ft: pd.DataFrame) -> pd.DataFrame:
    """Bray-Curtis on relative abundance (ASV × sample)."""
    from scipy.spatial.distance import pdist

    rel = ft.div(ft.sum(axis=0), axis=1).fillna(0).to_numpy().T
    D = pdist(rel, metric="braycurtis")
    return pd.DataFrame(squareform(D), index=ft.columns, columns=ft.columns)


def pcoa(distance: pd.DataFrame, n_components: int = 4) -> tuple[pd.DataFrame, np.ndarray]:
    """Principal coordinates analysis on a square distance matrix.

    Returns (samples × components) DataFrame and the eigenvalue fractions.
    """
    D2 = distance.to_numpy() ** 2
    n = D2.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals, evecs = evals[idx], evecs[:, idx]
    positive = evals > 0
    k = min(n_components, positive.sum())
    coords = evecs[:, :k] * np.sqrt(evals[:k])
    frac = evals[:k] / evals[positive].sum()
    out = pd.DataFrame(
        coords,
        index=distance.index,
        columns=[f"PCo{i+1}" for i in range(k)],
    )
    return out, frac


def permanova_by_term(
    distance: pd.DataFrame, meta: pd.DataFrame, terms: list[str], n_perm: int = 999,
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Per-term PERMANOVA F and R² on a distance matrix.

    This is a simple per-term implementation; for marginal (type-III)
    partitioning use R's vegan::adonis2.
    """
    from skbio.stats.distance import DistanceMatrix, permanova

    arr = np.ascontiguousarray(distance.to_numpy(), dtype=float)
    # Symmetrise rounding-error asymmetry
    arr = (arr + arr.T) / 2
    np.fill_diagonal(arr, 0.0)
    dm = DistanceMatrix(arr, ids=list(distance.index))
    rows = []
    for term in terms:
        grp = meta.loc[list(distance.index), term].astype(str)
        try:
            res = permanova(dm, grouping=grp, permutations=n_perm)
            rows.append({
                "term": term,
                "F": float(res["test statistic"]),
                "p": float(res["p-value"]),
                "n_groups": int(res["number of groups"]),
                "n": int(res["sample size"]),
            })
        except Exception as exc:  # pragma: no cover
            rows.append({"term": term, "F": np.nan, "p": np.nan,
                         "n_groups": np.nan, "n": np.nan, "error": str(exc)})
    return pd.DataFrame(rows)
