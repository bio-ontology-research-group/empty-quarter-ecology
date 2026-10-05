import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description="Replay the author's within-plate position permutation test.")
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
WELLS_FILE = args.output/"wells.tsv"
BC_FILE = args.output/"bc_matrix.tsv"
META_FILE = args.metadata
STRATA_OUT = args.output/"strata_summary.tsv"
NULL_OUT = args.output/"null_distribution.tsv"
EXCESS_OUT = args.output/"well_excess_similarity.tsv"
SUMMARY_OUT = args.output/"summary.txt"
result_records = []
N_PERM = 999
NEIGHBOUR_DEFS = {"ortho+diag": 1.5, "ortho_only": 1.0}  # radius: d<=sqrt(2) vs d=1 only
MIN_WELLS = 10
SEED = 42

# --- Load ---
wells_all = pd.read_csv(WELLS_FILE, sep="\t")
wells_all = wells_all.merge(pd.read_csv(META_FILE, sep="\t")[["profile", "trip"]], left_on="sid", right_on="profile").drop(columns="profile")
bc = pd.read_csv(BC_FILE, sep="\t", index_col=0)

def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))

# --- Analysis for one trip ---
def analyse(wells, label, radius, ndef, rng):
    wells = wells.reset_index(drop=True)
    n = len(wells)
    bc_w = bc.loc[wells.sid, wells.sid].to_numpy()

    # within-plate, different-site pairs; strata = geo-decile x same-compartment, deciles within this trip
    i_idx, j_idx = np.triu_indices(n, 1)
    keep = (wells.Plate.values[i_idx] == wells.Plate.values[j_idx]) & (wells.site.values[i_idx] != wells.site.values[j_idx])
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    pairs = pd.DataFrame({
        "i": i_idx, "j": j_idx, "bc": bc_w[i_idx, j_idx],
        "geo_km": haversine(wells.lat.values[i_idx], wells.lon.values[i_idx], wells.lat.values[j_idx], wells.lon.values[j_idx]),
        "same_comp": wells.compartment.values[i_idx] == wells.compartment.values[j_idx],
    })
    pairs["geo_decile"] = pd.qcut(pairs.geo_km, 10, labels=False, duplicates="drop")
    pairs["stratum"] = pairs.geo_decile.astype(str) + "_" + np.where(pairs.same_comp, "same", "diff")
    strat_codes, strat_labels = pd.factorize(pairs.stratum)
    n_strata = len(strat_labels)
    rows, cols = wells.row.values.astype(float), wells.col.values.astype(float)
    bc_v = pairs.bc.values

    def adjacency(perm):
        d = np.hypot(rows[perm[i_idx]] - rows[perm[j_idx]], cols[perm[i_idx]] - cols[perm[j_idx]])
        return d <= radius

    def statistic(adj):
        sum_adj = np.bincount(strat_codes, weights=bc_v * adj, minlength=n_strata)
        n_adj = np.bincount(strat_codes, weights=adj, minlength=n_strata)
        sum_non = np.bincount(strat_codes, weights=bc_v * ~adj, minlength=n_strata)
        n_non = np.bincount(strat_codes, weights=~adj, minlength=n_strata)
        ok = (n_adj > 0) & (n_non > 0)
        delta = sum_adj[ok] / n_adj[ok] - sum_non[ok] / n_non[ok]
        return np.average(delta, weights=n_adj[ok]), n_adj, n_non

    # observed + null (shuffle well positions within each plate)
    identity = np.arange(n)
    adj_obs = adjacency(identity)
    T_obs, n_adj_obs, n_non_obs = statistic(adj_obs)
    pairs["adjacent"] = adj_obs
    plate_groups = [np.flatnonzero(wells.Plate.values == p) for p in sorted(wells.Plate.unique())]
    T_null = np.empty(N_PERM)
    for k in range(N_PERM):
        perm = identity.copy()
        for g in plate_groups:
            perm[g] = rng.permutation(g)
        T_null[k] = statistic(adjacency(perm))[0]
    p_one = (1 + np.sum(T_null <= T_obs)) / (1 + N_PERM)
    z = (T_obs - T_null.mean()) / T_null.std(ddof=1)

    strata = pd.DataFrame({"trip": label, "neighbours": ndef, "stratum": strat_labels, "n_adjacent": n_adj_obs.astype(int), "n_nonadjacent": n_non_obs.astype(int)})
    strata["mean_bc_adjacent"] = strata.stratum.map(pairs[pairs.adjacent].groupby("stratum").bc.mean())
    strata["mean_bc_nonadjacent"] = strata.stratum.map(pairs[~pairs.adjacent].groupby("stratum").bc.mean())
    strata["delta"] = strata.mean_bc_adjacent - strata.mean_bc_nonadjacent
    strata["geo_km_range"] = strata.stratum.map(pairs.groupby("stratum").geo_km.agg(lambda x: f"{x.min():.0f}-{x.max():.0f}"))

    # per-well excess similarity: mean BC(matched non-neighbours) - mean BC(neighbours)
    excess = np.full(n, np.nan)
    n_neigh = np.zeros(n, dtype=int)
    n_fallback = 0
    adj_pairs, non_pairs = pairs[pairs.adjacent], pairs[~pairs.adjacent]
    for w in range(n):
        mine = adj_pairs[(adj_pairs.i == w) | (adj_pairs.j == w)]
        if mine.empty:
            continue
        n_neigh[w] = len(mine)
        nonadj = non_pairs[(non_pairs.i == w) | (non_pairs.j == w)]
        diffs = []
        for _, r in mine.iterrows():
            matched = nonadj[nonadj.stratum == r.stratum]
            if matched.empty:
                matched = nonadj[nonadj.same_comp == r.same_comp]  # fallback: drop the geo-decile constraint
                n_fallback += 1
            if not matched.empty:
                diffs.append(matched.bc.mean() - r.bc)
        if diffs:
            excess[w] = np.mean(diffs)
    wells = wells.assign(n_neighbors=n_neigh, excess_similarity=excess)
    pseudo = wells.dna_conc[wells.dna_conc > 0].min() / 2
    wells["log10_dna"] = np.log10(wells.dna_conc + pseudo)
    d = wells.dropna(subset=["excess_similarity", "dna_conc"])
    rho, p_rho = spearmanr(d.excess_similarity, d.dna_conc)
    r, p_r = pearsonr(d.excess_similarity, d.log10_dna)
    d_pos = d[d.dna_conc > 0]
    r_pos, p_r_pos = pearsonr(d_pos.excess_similarity, np.log10(d_pos.dna_conc))

    lines = [
        f"=== TRIP {label} | neighbours: {ndef} (d<={radius}) ===",
        f"Wells: {n} ({wells.sid.nunique()} unique samples, plates {sorted(wells.Plate.unique().tolist())})",
        f"Within-plate different-site pairs: {len(pairs)} (adjacent: {adj_obs.sum()}, non-adjacent: {(~adj_obs).sum()}); strata: {n_strata}",
        "Global permutation test",
        f"  Observed T (weighted mean BC_adj - BC_nonadj): {T_obs:+.4f}",
        f"  Null ({N_PERM} within-plate position permutations): mean {T_null.mean():+.4f}, sd {T_null.std(ddof=1):.4f}, z = {z:.2f}",
        f"  One-sided p (adjacent less dissimilar): {p_one:.3f}",
        "Low-yield vulnerability diagnostic",
        f"  Wells with excess score and DNA conc: {len(d)} ({(d.dna_conc == 0).sum()} with conc = 0); {n_fallback}/{adj_obs.sum()} neighbour pairs used compartment-only matching",
        f"  Spearman rho (excess vs dna_conc): {rho:+.3f}, p = {p_rho:.3f}",
        f"  Pearson r (excess vs log10(conc + {pseudo:.1e})): {r:+.3f}, p = {p_r:.3f}",
        f"  Pearson r, conc > 0 only (n={len(d_pos)}): {r_pos:+.3f}, p = {p_r_pos:.3f}",
        "",
    ]
    result_records.append({
        "trip": int(label), "neighbours": ndef, "radius": radius,
        "samples": n, "plates": sorted(wells.Plate.unique().tolist()),
        "pairs": len(pairs), "adjacent_pairs": int(adj_obs.sum()),
        "strata": n_strata, "delta": float(T_obs), "p_one_sided": float(p_one),
        "permutations": N_PERM,
    })
    return lines, strata, pd.DataFrame({"trip": label, "neighbours": ndef, "T_null": T_null}), wells.assign(neighbours=ndef)

# --- Run per trip ---
rng = np.random.default_rng(SEED)
out_lines, out_strata, out_null, out_wells = [], [], [], []
for trip, grp in wells_all.groupby("trip"):
    if len(grp) < MIN_WELLS:
        out_lines.append(f"=== TRIP {trip} === skipped ({len(grp)} wells)\n")
        continue
    for ndef, radius in NEIGHBOUR_DEFS.items():
        lines, strata, null, w = analyse(grp, trip, radius, ndef, rng)
        out_lines += lines
        out_strata.append(strata); out_null.append(null); out_wells.append(w)

summary = "\n".join(out_lines)
print(summary)
with open(SUMMARY_OUT, "w") as f:
    f.write(summary + "\n")
pd.concat(out_strata).to_csv(STRATA_OUT, sep="\t", index=False, float_format="%.4f")
pd.concat(out_null).to_csv(NULL_OUT, sep="\t", index=False, float_format="%.6f")
pd.concat(out_wells).to_csv(EXCESS_OUT, sep="\t", index=False, float_format="%.5f")

(args.output/"results.json").write_text(json.dumps({
    "seed": SEED, "permutations": N_PERM,
    "matching": "within-trip geographic distance decile x same/different compartment",
    "permutation": "sample-to-well assignment within each plate",
    "results": result_records,
}, indent=2)+"\n")
