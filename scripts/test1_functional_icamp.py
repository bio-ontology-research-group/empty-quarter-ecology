#!/usr/bin/env python3
"""TEST 1: Functional iCAMP decoupling.

Compare process attribution at the PICRUSt2-functional level vs the
taxonomic level. iCAMP needs a phylogeny for betaNTI; for functions we
substitute a "functional bMNTD" using KEGG pathway BC (Stegen 2013's
Bray-Curtis-based RCbray captures the stochastic-vs-deterministic axis
without a tree). We then fit a 3-class (variable-selection / homogeneous-
selection / stochastic) classifier using:
  - Functional NTI proxy: Mantel-derived "z-score" against the abundance-
    randomized null
  - RCbray for stochastic split

Simpler approach used here:
  - Per (compartment) compute pairwise-functional-BC and pairwise-functional
    RCbray (999 perm null preserving sample richness, abundance-pool weights).
  - Compare RCbray distribution to taxonomic RCbray (same 5-class iCAMP we
    already ran).
  - Mantel: functional-BC vs taxonomic-BC (does taxonomy predict function?)
  - Mantel partial: functional-BC vs environment | distance — does
    environment select functions more than taxa?

Outputs:
  cache/test1_functional_icamp/RCbray_functional_{compartment}.parquet
  cache/test1_functional_icamp/process_attrib_functional_{compartment}.tsv
  cache/test1_functional_icamp/taxonomy_vs_function_mantel.tsv
  cache/test1_functional_icamp/summary.txt
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _sample_parse import parse_samples_to_df

CACHE = REPO / "cache"
DATA = REPO / "data"
OUT = CACHE / "test1_functional_icamp"
OUT.mkdir(parents=True, exist_ok=True)

N_PERM = 999
N_TOP_PATHWAYS = 200   # match the order-of-magnitude of taxonomic top-500 ASVs
RNG = np.random.default_rng(20260509)


def parse_compartment(sid: str) -> str:
    import re
    s = sid.split("_")[-1]
    m = re.match(r"^[A-Z]*[0-9]+([A-Z]+)r[0-9]+", s)
    if not m:
        return "?"
    return {"PR": "rhizosphere", "S": "surface", "D": "deep"}.get(m.group(1), "?")


def select_top(table: pd.DataFrame, n: int) -> pd.DataFrame:
    totals = table.sum(axis=1)
    keep = totals.sort_values(ascending=False).head(n).index
    return table.loc[keep]


def rc_bray_matrix(A: np.ndarray, obs_bc: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    m = A.shape[1]
    sample_totals = A.sum(axis=1)
    occ = (A > 0).mean(axis=0)
    occ = np.where(occ > 0, occ, 1e-9); occ = occ / occ.sum()
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
            if r == 0: continue
            chosen = RNG.choice(m, size=int(r), replace=False, p=occ)
            w = abund_pool[chosen]; w = w / w.sum() * sample_totals[s]
            Aperm[s, chosen] = w
        d = squareform(pdist(Aperm, metric="braycurtis"))
        above[d > obs_bc] += 1
        below[d < obs_bc] += 1
        if (p + 1) % 100 == 0:
            print(f"    perm {p+1}/{N_PERM}", flush=True)
    ties = N_PERM - above - below
    rc = 2.0 * (above + 0.5 * ties) / N_PERM - 1.0
    np.fill_diagonal(rc, np.nan)
    return rc


def main():
    # Load PICRUSt2 pathways
    path = pd.read_csv(DATA / "functional" / "picrust2" / "path_abun_unstrat.tsv",
                       sep="\t", index_col=0)
    print(f"path table: {path.shape} (pathways x samples)", flush=True)

    # Build sample compartment map for functional table
    samp_comp = {s: parse_compartment(s) for s in path.columns}
    samp_comp_df = pd.Series(samp_comp).to_frame("compartment")
    print("functional samples per compartment:")
    print(samp_comp_df["compartment"].value_counts().to_string())

    summaries = []
    for comp in ["rhizosphere", "surface", "deep"]:
        cols = [s for s, c in samp_comp.items() if c == comp]
        sub = path[cols].copy()
        sub = select_top(sub, N_TOP_PATHWAYS)
        rel = sub.div(sub.sum(axis=0).replace(0, 1), axis=1).T.values
        obs_bc = squareform(pdist(rel, metric="braycurtis"))
        print(f"\n[{comp}] computing functional RCbray ({N_PERM} perms, "
              f"n={len(cols)}, m={rel.shape[1]})", flush=True)
        rc = rc_bray_matrix(rel, obs_bc)
        rc_df = pd.DataFrame(rc, index=cols, columns=cols)
        rc_df.to_parquet(OUT / f"RCbray_functional_{comp}.parquet")

        # Classify pairs by RCbray alone (no betaNTI for functions —
        # use deterministic-vs-stochastic split):
        #   |RCbray| < 0.95  -> stochastic (drift)
        #   RCbray > +0.95   -> dispersal limitation / divergence
        #   RCbray < -0.95   -> homogenizing dispersal
        iu = np.triu_indices(len(cols), k=1)
        rc_vals = rc[iu]
        valid = np.isfinite(rc_vals)
        rcv = rc_vals[valid]
        attribs = {
            "homogenizing_dispersal": float(np.mean(rcv < -0.95)),
            "dispersal_limitation":  float(np.mean(rcv > 0.95)),
            "drift_or_weak":         float(np.mean(np.abs(rcv) <= 0.95)),
            "n_pairs": int(valid.sum()),
        }
        summaries.append({"compartment": comp, **attribs})
        with open(OUT / f"process_attrib_functional_{comp}.tsv", "w") as fh:
            fh.write("class\tfraction\n")
            for k, v in attribs.items():
                fh.write(f"{k}\t{v}\n")

    sdf = pd.DataFrame(summaries)
    sdf.to_csv(OUT / "process_attrib_functional_summary.tsv", sep="\t", index=False)
    print("\n=== FUNCTIONAL RCbray attribution ===")
    print(sdf.to_string(index=False))

    # Compare to taxonomic iCAMP (existing)
    tax_summary = pd.read_csv(CACHE / "icamp" / "process_summary_all.tsv",
                               sep="\t")
    print("\n=== TAXONOMIC iCAMP (existing) ===")
    print(tax_summary.to_string(index=False))

    # Compute Mantel: functional BC vs taxonomic BC per compartment
    # (using site-aggregated BC matrices we already have)
    bc_func_path = OUT / "bc_functional_per_compartment.parquet"
    func_bc_records = []
    # Aggregate per (site, trip, compartment) the functional vector
    smeta = parse_samples_to_df(path.columns)
    smeta["site"] = smeta["site"].astype(int)
    rel_path = path.div(path.sum(axis=0).replace(0, 1), axis=1)
    func_comm = {}
    for (s, t, c), grp in smeta.groupby(["site", "trip", "compartment"]):
        cols = grp["sample"].tolist()
        if len(cols) == 0: continue
        v = rel_path[cols].mean(axis=1).values
        func_comm[(int(s), int(t), c)] = v

    # Pull existing taxonomic BC
    bc_tax = pd.read_parquet(CACHE / "wind_dispersal" / "bc_per_compartment_trip.parquet")

    def mantel(x, y, n_perm=999):
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() < 10:
            return (np.nan, np.nan, int(valid.sum()))
        x = x[valid]; y = y[valid]
        r_obs = np.corrcoef(x, y)[0, 1]
        cnt = 0
        for _ in range(n_perm):
            yp = RNG.permutation(y)
            if abs(np.corrcoef(x, yp)[0, 1]) >= abs(r_obs):
                cnt += 1
        return (float(r_obs), (cnt + 1) / (n_perm + 1), int(valid.sum()))

    mantel_rows = []
    for comp in ["rhizosphere", "surface", "deep"]:
        for trip in range(1, 6):
            tax_sub = bc_tax[(bc_tax["compartment"] == comp) & (bc_tax["trip"] == trip)]
            # Build matched functional BC for the same site pairs
            sites = sorted(set(tax_sub["site_i"]).union(tax_sub["site_j"]))
            n_s = len(sites)
            if n_s < 5: continue
            sidx = {s: i for i, s in enumerate(sites)}
            bc_func = np.full((n_s, n_s), np.nan)
            for i in range(n_s):
                vi = func_comm.get((sites[i], trip, comp))
                if vi is None: continue
                for j in range(i + 1, n_s):
                    vj = func_comm.get((sites[j], trip, comp))
                    if vj is None: continue
                    num = np.abs(vi - vj).sum()
                    den = vi.sum() + vj.sum()
                    bc_func[i, j] = bc_func[j, i] = num / den if den > 0 else np.nan
            np.fill_diagonal(bc_func, 0.0)
            # Build matched taxonomic BC
            bc_taxm = np.full((n_s, n_s), np.nan)
            for _, r in tax_sub.iterrows():
                if int(r["site_i"]) in sidx and int(r["site_j"]) in sidx:
                    i = sidx[int(r["site_i"])]; j = sidx[int(r["site_j"])]
                    bc_taxm[i, j] = bc_taxm[j, i] = float(r["bc"])
            iu = np.triu_indices(n_s, k=1)
            r, p, n = mantel(bc_func[iu], bc_taxm[iu])
            mantel_rows.append({"compartment": comp, "trip": trip,
                                 "n_pairs": n, "r_func_vs_tax_BC": r,
                                 "p": p})
    mdf = pd.DataFrame(mantel_rows)
    mdf.to_csv(OUT / "taxonomy_vs_function_mantel.tsv", sep="\t", index=False)
    print("\n=== Mantel: functional BC vs taxonomic BC per (comp, trip) ===")
    print(mdf.to_string(index=False))

    # Summary
    with open(OUT / "summary.txt", "w") as fh:
        fh.write("Test 1: Functional iCAMP decoupling\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"PICRUSt2 pathway table: {path.shape}\n")
        fh.write(f"Top {N_TOP_PATHWAYS} pathways used (by total abundance)\n")
        fh.write(f"Permutations per RCbray: {N_PERM}\n\n")

        fh.write("FUNCTIONAL RCbray attribution (% of pairs):\n")
        fh.write(sdf.to_string(index=False))

        fh.write("\n\nTAXONOMIC iCAMP attribution (existing, for comparison):\n")
        fh.write(tax_summary[["compartment", "n_pairs",
                               "frac_homogenizing_dispersal",
                               "frac_dispersal_limitation",
                               "frac_drift_or_weak",
                               "frac_homogeneous_selection",
                               "frac_variable_selection"]].to_string(index=False))

        fh.write("\n\nDECOUPLING TEST:\n")
        for c in ["rhizosphere", "surface", "deep"]:
            fmech = sdf[sdf["compartment"] == c].iloc[0]
            tmech = tax_summary[tax_summary["compartment"] == c].iloc[0]
            f_homog = fmech["homogenizing_dispersal"]
            t_homog = tmech["frac_homogenizing_dispersal"]
            t_select = (tmech["frac_homogeneous_selection"]
                         + tmech["frac_variable_selection"])
            f_select_proxy = 1.0 - fmech["homogenizing_dispersal"] \
                            - fmech["dispersal_limitation"] - fmech["drift_or_weak"]
            # functional only has 3 classes — can't measure selection same way
            fh.write(f"  {c:>11s}: TAXONOMIC homog-disp={t_homog:.1%} | "
                     f"FUNCTIONAL homog-disp={f_homog:.1%} | "
                     f"diff={(f_homog - t_homog):+.1%}\n")

        fh.write("\nMantel functional vs taxonomic BC (per (comp, trip)):\n")
        fh.write(mdf.to_string(index=False))
        fh.write("\n\nINTERPRETATION KEY:\n")
        fh.write("  if functional homog-disp >> taxonomic: functions move with bodies\n")
        fh.write("  if functional homog-disp << taxonomic: functions are MORE selected\n"
                 "    than taxa (the surprising 'bodies dispersed, genes selected' case)\n")
        fh.write("  Mantel r close to 1: function tracks taxonomy (no decoupling)\n")
        fh.write("  Mantel r << 1:        function decoupled from taxonomy\n")
    print(f"\nWrote {OUT}/summary.txt")


if __name__ == "__main__":
    main()
