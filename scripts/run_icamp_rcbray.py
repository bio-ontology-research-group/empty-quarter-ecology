#!/usr/bin/env python3
"""Tier-1 #1: Complete iCAMP -- compute Raup-Crick on Bray-Curtis (RCbray)
to break the stochastic class into homogenizing dispersal vs dispersal
limitation vs drift. Uses pre-computed per-compartment betaNTI matrices
in cache/bnti/.

Inputs:
  cache/bnti/bNTI_{rhizosphere,surface,deep}.parquet  -- 474x474 etc square
                                                          matrices (NaN diag)
  cache/feature_table.parquet  -- counts (ASV x sample)

Outputs:
  cache/icamp/RCbray_{compartment}.parquet
  cache/icamp/process_attribution_{compartment}.tsv
  cache/icamp/process_summary_all.tsv
  cache/icamp/process_summary_all.txt
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "cache"
OUT = CACHE / "icamp"
OUT.mkdir(parents=True, exist_ok=True)

N_PERM = 999
N_ASV_TOP = 500  # match the existing bNTI feature set
RNG = np.random.default_rng(20260509)


def select_top_asvs(ft_sub: pd.DataFrame, n: int) -> pd.DataFrame:
    """Keep top-n ASVs by total abundance within the compartment subset."""
    totals = ft_sub.sum(axis=1)
    keep = totals.sort_values(ascending=False).head(n).index
    return ft_sub.loc[keep]


def rc_bray_matrix(A: np.ndarray, obs_bc: np.ndarray) -> np.ndarray:
    """Vectorized RCbray: for each pair, count perms with BC > / == / < obs.

    Args:
      A:       (n_samples, n_asvs) RELATIVE abundance matrix (rows sum to 1)
      obs_bc:  (n_samples, n_samples) observed Bray-Curtis distances.

    Returns:
      rc:      (n_samples, n_samples) RCbray in [-1, 1]
    """
    n = A.shape[0]
    m = A.shape[1]
    # Marginal counts: total reads per sample (we'll preserve), and per-ASV
    # occupancy weight for sampling.
    sample_totals = A.sum(axis=1)
    occ = (A > 0).mean(axis=0)
    occ = np.where(occ > 0, occ, 1e-9)
    occ = occ / occ.sum()
    abund_pool = A.mean(axis=0)
    abund_pool = np.where(abund_pool > 0, abund_pool, 1e-9)
    abund_pool = abund_pool / abund_pool.sum()

    sample_richness = (A > 0).sum(axis=1).astype(int)

    above = np.zeros((n, n), dtype=np.int32)
    below = np.zeros((n, n), dtype=np.int32)

    for p in range(N_PERM):
        Aperm = np.zeros_like(A)
        for s in range(n):
            r = sample_richness[s]
            if r == 0:
                continue
            chosen = RNG.choice(m, size=int(r), replace=False, p=occ)
            w = abund_pool[chosen]
            w = w / w.sum() * sample_totals[s]
            Aperm[s, chosen] = w

        # Pairwise BC: pdist returns condensed; squareform expands
        d = squareform(pdist(Aperm, metric="braycurtis"))
        # Compare to obs only on upper triangle (incl diagonal NaN-safe)
        gt = d > obs_bc
        lt = d < obs_bc
        above[gt] += 1
        below[lt] += 1
        if (p + 1) % 100 == 0:
            print(f"  perm {p+1}/{N_PERM}", flush=True)

    ties = N_PERM - above - below
    rc = 2.0 * (above + 0.5 * ties) / N_PERM - 1.0
    np.fill_diagonal(rc, np.nan)
    return rc


def parse_compartment(sid: str) -> str:
    import re
    s = sid.split("_")[-1]
    m = re.match(r"^[A-Z]*[0-9]+([A-Z]+)r[0-9]+$", s)
    if not m:
        return "?"
    return {"PR": "rhizosphere", "S": "surface", "D": "deep"}.get(m.group(1), "?")


def attribute(comp: str, bnti_pq: Path, rc_mat: pd.DataFrame, samples: list) -> dict:
    bnti = pd.read_parquet(bnti_pq)
    bnti = bnti.reindex(index=samples, columns=samples)
    rc = rc_mat.reindex(index=samples, columns=samples)
    rec = []
    n = len(samples)
    for i in range(n):
        for j in range(i + 1, n):
            bn = bnti.iat[i, j]
            r = rc.iat[i, j]
            if pd.isna(bn) or pd.isna(r):
                continue
            if bn < -2:
                proc = "homogeneous_selection"
            elif bn > 2:
                proc = "variable_selection"
            elif r < -0.95:
                proc = "homogenizing_dispersal"
            elif r > 0.95:
                proc = "dispersal_limitation"
            else:
                proc = "drift_or_weak"
            rec.append({"sample_i": samples[i], "sample_j": samples[j],
                        "bNTI": float(bn), "RCbray": float(r), "process": proc})
    df = pd.DataFrame(rec)
    df.to_csv(OUT / f"process_attribution_{comp}.tsv", sep="\t", index=False)
    counts = df["process"].value_counts()
    fracs = counts / counts.sum()
    summary = {"compartment": comp, "n_pairs": int(len(df))}
    for k in ["homogeneous_selection", "variable_selection",
              "homogenizing_dispersal", "dispersal_limitation", "drift_or_weak"]:
        summary[f"frac_{k}"] = float(fracs.get(k, 0.0))
        summary[f"n_{k}"] = int(counts.get(k, 0))
    return summary


def main():
    ft = pd.read_parquet(CACHE / "feature_table.parquet")
    print(f"feature_table: {ft.shape}", flush=True)

    summaries = []
    for comp in ["rhizosphere", "surface", "deep"]:
        bnti_pq = CACHE / "bnti" / f"bNTI_{comp}.parquet"
        if not bnti_pq.exists():
            print(f"[{comp}] no bNTI parquet -- skipping")
            continue

        # The bNTI sample set
        bnti_idx = pd.read_parquet(bnti_pq).index.tolist()
        # keep only samples present in feature_table
        samples = [s for s in bnti_idx if s in ft.columns]
        print(f"\n[{comp}] samples in bNTI ∩ feature_table: {len(samples)}", flush=True)
        sub = ft[samples]
        sub = select_top_asvs(sub, N_ASV_TOP)
        # Convert to relative
        rel = sub.div(sub.sum(axis=0).replace(0, 1), axis=1).T.values  # samples x ASVs
        # Observed BC
        obs_bc = squareform(pdist(rel, metric="braycurtis"))

        rc_path = OUT / f"RCbray_{comp}.parquet"
        if rc_path.exists():
            print(f"[{comp}] RCbray parquet exists -- loading")
            rc_mat = pd.read_parquet(rc_path)
        else:
            print(f"[{comp}] computing RCbray ({N_PERM} perms, "
                  f"n={len(samples)}, m={rel.shape[1]})", flush=True)
            rc = rc_bray_matrix(rel, obs_bc)
            rc_mat = pd.DataFrame(rc, index=samples, columns=samples)
            rc_mat.to_parquet(rc_path)
            print(f"[{comp}] -> {rc_path.name}", flush=True)

        summaries.append(attribute(comp, bnti_pq, rc_mat, samples))

    if summaries:
        sdf = pd.DataFrame(summaries)
        sdf.to_csv(OUT / "process_summary_all.tsv", sep="\t", index=False)
        with open(OUT / "process_summary_all.txt", "w") as fh:
            fh.write("iCAMP process attribution by compartment (Tier-1 #1)\n")
            fh.write("=" * 70 + "\n\n")
            fh.write(sdf.to_string(index=False))
            fh.write("\n\n")
            fh.write("Process classes:\n")
            fh.write("  variable_selection      |bNTI|>2, bNTI>2  (heterogeneous selection)\n")
            fh.write("  homogeneous_selection   |bNTI|>2, bNTI<-2 (uniform selection)\n")
            fh.write("  homogenizing_dispersal  |bNTI|<2, RCbray<-0.95\n")
            fh.write("  dispersal_limitation    |bNTI|<2, RCbray>+0.95\n")
            fh.write("  drift_or_weak           |bNTI|<2, |RCbray|<0.95\n")
        print("\n=== Process summary ===")
        print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
