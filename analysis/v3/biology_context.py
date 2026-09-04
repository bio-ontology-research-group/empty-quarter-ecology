#!/usr/bin/env python3
"""Descriptive biology behind the statistical results of the ecology manuscript.

Second round of the co-author review of 2 Sep 2026 (R. Gruenberg): the paper
reported statistics without organisms or substrates in several places.  This
module adds the missing descriptions from the canonical inputs and the
committed canonical result bundles.  Every analysis is descriptive or a
declared multiple-testing family; none changes a tracked conclusion.

Analyses
--------
1. Landform.  Each core site carries one field landform label (sand dune,
   saline pan, desert oasis, gravel, dune slack, aeolian lake, oil spill).
   Site-level alpha diversity by landform; composition of the 44 dune versus
   the 9 saline-pan sites in the geographic-model coordinates (omnibus
   pseudo-F with 9,999 label permutations, per-genus Mann-Whitney tests with
   Benjamini-Hochberg over 200 genera), with and without adjustment for route
   position; and the dune-only sensitivity of the route model, the route
   genus correlations and the climate-diversity correlations.
2. Compartment genus family.  Per-genus paired contrasts between compartments
   in the block-centred site tensor of the canonical compartment analysis
   (200 genera x 3 contrasts = 600 sign-flip tests, one BH family), with a
   site bootstrap interval.  The mean differences must reproduce the
   committed displacement loadings.
3. Compartment pathways.  Leading supported PICRUSt2 pathway contrasts per
   compartment pair and direction (from the committed pathway family).
4. Laboratory-XRF axis.  Element loadings of the committed elemental axis and
   the genera whose site-level CLR values track the axis (Spearman, BH over
   200).
5. pH.  Genera whose site-level CLR values track site-mean archived-soil pH
   (Spearman, BH over 200).
6. Core genera per compartment.  After pooling campaigns within site and
   compartment and subsampling to 12,865 reads, the genera detected at
   >= 90 % of sites in each compartment, as a descriptive explanation of the
   compartment difference in distance decay.

Relative abundances are shares of reads; CLR values are centred log ratios
over the 200 primary genera; SILVA 138.2 names.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compartment_composition_rescue as ccr  # noqa: E402
import taxon_context as tc  # noqa: E402

SCHEMA_VERSION = "1.0"
POSITIONS = ("Surface", "Deep", "Rhizosphere")
LABEL = {"Surface": "surface", "Deep": "shallow_subsurface", "Rhizosphere": "root_adjacent"}
CONTRASTS = (("Deep", "Surface"), ("Rhizosphere", "Surface"), ("Rhizosphere", "Deep"))
CORE_DEPTH = 12865
CORE_OCCUPANCY = 0.90


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.10g", lineterminator="\n")


def bh(p: np.ndarray) -> np.ndarray:
    return tc.benjamini_hochberg(np.asarray(p, dtype=float))


def route_r2(site_means: pd.DataFrame, position: np.ndarray) -> float:
    """Share of the summed site-level CLR variance explained by a linear plus
    quadratic route model (the pre-specified geographic model)."""
    design = np.column_stack([np.ones_like(position), position, position**2])
    y = site_means.to_numpy(dtype=float)
    y = y - y.mean(axis=0)
    beta, *_ = np.linalg.lstsq(design - design.mean(axis=0), y, rcond=None)
    fitted = (design - design.mean(axis=0)) @ beta
    return float(1 - np.square(y - fitted).sum() / np.square(y).sum())


def residualise(site_means: pd.DataFrame, position: np.ndarray) -> pd.DataFrame:
    design = np.column_stack([np.ones_like(position), position, position**2])
    y = site_means.to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return pd.DataFrame(y - design @ beta, index=site_means.index, columns=site_means.columns)


def group_pseudo_f(matrix: np.ndarray, labels: np.ndarray) -> float:
    grand = matrix.mean(axis=0)
    ss_total = float(np.square(matrix - grand).sum())
    ss_within = 0.0
    for label in np.unique(labels):
        block = matrix[labels == label]
        ss_within += float(np.square(block - block.mean(axis=0)).sum())
    ss_between = ss_total - ss_within
    k = len(np.unique(labels))
    n = matrix.shape[0]
    return (ss_between / (k - 1)) / (ss_within / (n - k))


def landform_omnibus(matrix: np.ndarray, labels: np.ndarray, permutations: int, rng: np.random.Generator) -> tuple[float, float]:
    observed = group_pseudo_f(matrix, labels)
    exceed = 0
    for _ in range(permutations):
        exceed += group_pseudo_f(matrix, rng.permutation(labels)) >= observed
    return observed, (exceed + 1) / (permutations + 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=HERE.parents[1])
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--core-draws", type=int, default=20)
    args = parser.parse_args()

    root = args.project_root.resolve()
    output = (args.output_dir or root / "analysis/v3/biology_context").resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    inputs = {
        "genus_counts": root / "analysis/v2/review/cache/genus_counts.tsv",
        "sample_metadata": root / "data/processed/functional/picrust2/merged/sample_metadata.tsv",
        "site_coordinates": root / "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv",
        "climate_site_summary": root / "analysis/v3/environment_associations/climate_site_summary.tsv",
        "climate_alpha_correlations": root / "analysis/v3/environment_associations/climate_alpha_correlations.tsv",
        "transect_replacement": root / "analysis/v3/taxon_context/transect_replacement.tsv",
        "genus_composition": root / "analysis/v3/taxon_context/genus_composition.tsv",
        "displacement_loadings": root / "analysis/v3/compartment_composition/paired_displacement_loadings.tsv",
        "pathway_position_effects": root / "analysis/v3/picrust2_ecology/pathway_position_effects.tsv",
        "xrf_axis": root / "analysis/v3/xrf_community_rescue/laboratory_xrf_axis.tsv",
        "xrf_loadings": root / "analysis/v3/xrf_community_rescue/elemental_pc1_loadings.tsv",
        "ph_sample_profile_join": root / "analysis/v3/ph_shared_v1/ecology/ph_sample_profile_join.tsv",
    }
    missing = [str(p) for p in inputs.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing inputs:\n" + "\n".join(missing))

    # ---- profiles, sites, landforms -------------------------------------
    genus = pd.read_csv(inputs["genus_counts"], sep="\t", index_col=0)
    genus = genus.loc[~genus.index.isna()]
    genus.index = genus.index.astype(str)
    meta = pd.DataFrame([m for sid in genus.columns if (m := tc.sample_metadata(sid)) is not None])
    core = meta[meta["site"].between(1, 60) & meta["position"].isin(POSITIONS)].copy()
    coordinates = pd.read_csv(inputs["site_coordinates"], sep="\t")
    coordinates = coordinates[coordinates["site"].between(1, 60)].sort_values("transect_km").reset_index(drop=True)
    ordered = coordinates["site"].astype(int).tolist()
    thirds = {s: tc.THIRDS[min(i // 20, 2)] for i, s in enumerate(ordered)}
    position_km = coordinates.set_index("site")["transect_km"]

    sample_meta = pd.read_csv(inputs["sample_metadata"], sep="\t")
    sample_meta = sample_meta[sample_meta["site"].between(1, 60)]
    landform_by_site = sample_meta.groupby("site")["feature"].agg(lambda v: sorted(set(v.dropna())))
    if (landform_by_site.map(len) != 1).any():
        raise ValueError("A core site carries more than one landform label")
    landform = landform_by_site.map(lambda v: v[0]).rename("landform")
    landform.index = landform.index.astype(int)
    site_table = pd.DataFrame(
        {
            "site": ordered,
            "transect_km": [float(position_km[s]) for s in ordered],
            "third": [thirds[s] for s in ordered],
            "landform": [landform[s] for s in ordered],
        }
    )
    write_tsv(site_table, output / "site_landforms.tsv")

    # ---- grouped counts, primary genera, site-level CLR means -----------
    grouped = tc.grouped_counts(genus, core, 2000)
    ranked = tc.rank_taxa(grouped, 0.20)
    primary = ranked[:200]
    replacement_ref = pd.read_csv(inputs["transect_replacement"], sep="\t")
    if set(replacement_ref["genus"]) != set(primary):
        raise ValueError("Primary genus set differs from the taxon-context module")
    _, site_means = tc.transect_replacement(grouped, primary, coordinates, thirds, 0.5)
    site_means = site_means.reindex(ordered)
    km = np.asarray([position_km[s] for s in site_means.index], dtype=float)

    # ---- 1. landform -------------------------------------------------------
    climate = pd.read_csv(inputs["climate_site_summary"], sep="\t").set_index("site")
    alpha_rows = []
    for name, block in site_table.groupby("landform"):
        sites = block["site"].tolist()
        alpha_rows.append(
            {
                "landform": name,
                "n_sites": len(sites),
                "sites": ",".join(map(str, sites)),
                "thirds": ",".join(f"{t}:{n}" for t, n in block["third"].value_counts().sort_index().items()),
                "mean_shannon": float(climate.loc[sites, "shannon"].mean()),
                "mean_expected_richness_25k": float(climate.loc[sites, "expected_richness_25k"].mean()),
                "mean_normalized_evenness": float(climate.loc[sites, "normalized_evenness"].mean()),
                "mean_air_temperature_c": float(climate.loc[sites, "mean_air_temperature_c"].mean()),
            }
        )
    alpha_by_landform = pd.DataFrame(alpha_rows).sort_values("n_sites", ascending=False, kind="mergesort")
    dune_sites = site_table.loc[site_table["landform"] == "sand dune", "site"].tolist()
    pan_sites = site_table.loc[site_table["landform"] == "saline pan", "site"].tolist()
    for column in ("shannon", "expected_richness_25k", "normalized_evenness"):
        u, p = stats.mannwhitneyu(climate.loc[dune_sites, column], climate.loc[pan_sites, column], alternative="two-sided")
        alpha_by_landform[f"dune_vs_pan_mannwhitney_p_{column}"] = p
    write_tsv(alpha_by_landform, output / "alpha_by_landform.tsv")

    two = site_means.loc[dune_sites + pan_sites]
    labels = np.asarray(["sand dune"] * len(dune_sites) + ["saline pan"] * len(pan_sites))
    km_two = np.asarray([position_km[s] for s in two.index], dtype=float)
    omnibus_rows = []
    per_genus_rows = []
    for adjustment, matrix in (
        ("none", two),
        ("route_linear_quadratic", residualise(two, km_two)),
    ):
        f_obs, p_obs = landform_omnibus(matrix.to_numpy(dtype=float), labels, args.permutations, rng)
        omnibus_rows.append(
            {
                "adjustment": adjustment,
                "n_dune_sites": len(dune_sites),
                "n_saline_pan_sites": len(pan_sites),
                "pseudo_f": f_obs,
                "permutation_p": p_obs,
                "permutations": args.permutations,
            }
        )
        pvals = []
        for g in primary:
            a = matrix.loc[dune_sites, g]
            b = matrix.loc[pan_sites, g]
            _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            pvals.append(p)
            per_genus_rows.append(
                {
                    "adjustment": adjustment,
                    "genus": g,
                    "mean_clr_dune": float(a.mean()),
                    "mean_clr_saline_pan": float(b.mean()),
                    "pan_minus_dune_clr": float(b.mean() - a.mean()),
                    "mannwhitney_p": p,
                }
            )
    per_genus = pd.DataFrame(per_genus_rows)
    for adjustment in per_genus["adjustment"].unique():
        mask = per_genus["adjustment"] == adjustment
        per_genus.loc[mask, "q_bh_200"] = bh(per_genus.loc[mask, "mannwhitney_p"].to_numpy())
    per_genus["supported_q_lt_0_05"] = per_genus["q_bh_200"] < 0.05
    composition = pd.read_csv(inputs["genus_composition"], sep="\t").set_index("genus")
    per_genus["phylum"] = per_genus["genus"].map(composition["phylum"]).fillna("")
    per_genus = per_genus.sort_values(["adjustment", "pan_minus_dune_clr"], kind="mergesort")
    write_tsv(pd.DataFrame(omnibus_rows), output / "landform_composition_omnibus.tsv")
    write_tsv(per_genus, output / "landform_genus_contrasts.tsv")

    # Dune-only sensitivity of the route results.
    dune_means = site_means.loc[dune_sites]
    km_dune = np.asarray([position_km[s] for s in dune_sites], dtype=float)
    route_rows = [
        {"cohort": "all_60_sites", "n_sites": 60, "route_r2": route_r2(site_means, km)},
    ]
    ref = replacement_ref.set_index("genus")
    dune_genus = []
    for g in primary:
        rho, p = stats.spearmanr(km_dune, dune_means[g])
        dune_genus.append({"genus": g, "phylum": composition["phylum"].get(g, ""), "spearman_rho_dune_sites": rho, "p_value": p,
                           "spearman_rho_all_sites": float(ref.loc[g, "spearman_rho_route_position"]),
                           "supported_all_sites": bool(ref.loc[g, "supported_q_lt_0_05"])})
    dune_genus = pd.DataFrame(dune_genus)
    dune_genus["q_bh_200"] = bh(dune_genus["p_value"].to_numpy())
    dune_genus["supported_dune_sites"] = dune_genus["q_bh_200"] < 0.05
    dune_genus["same_sign"] = np.sign(dune_genus["spearman_rho_dune_sites"]) == np.sign(dune_genus["spearman_rho_all_sites"])
    dune_genus = dune_genus.sort_values("spearman_rho_dune_sites", kind="mergesort")
    write_tsv(dune_genus, output / "route_genus_correlations_dune_sites.tsv")
    route_rows.append(
        {
            "cohort": "sand_dune_sites",
            "n_sites": len(dune_sites),
            "route_r2": route_r2(dune_means, km_dune),
            "supported_route_genera_all_sites": int(ref["supported_q_lt_0_05"].sum()),
            "supported_route_genera_dune_sites": int(dune_genus["supported_dune_sites"].sum()),
            "supported_in_both": int((dune_genus["supported_dune_sites"] & dune_genus["supported_all_sites"]).sum()),
            "supported_in_both_same_sign": int((dune_genus["supported_dune_sites"] & dune_genus["supported_all_sites"] & dune_genus["same_sign"]).sum()),
        }
    )
    route_table = pd.DataFrame(route_rows)
    write_tsv(route_table, output / "route_model_dune_sensitivity.tsv")

    climate_ref = pd.read_csv(inputs["climate_alpha_correlations"], sep="\t")
    clim_rows = []
    for variable in ("mean_air_temperature_c", "mean_monthly_rain_mm", "mean_relative_humidity_pct"):
        for measure in ("shannon", "expected_richness_25k", "normalized_evenness"):
            rho_all, p_all = stats.spearmanr(climate.loc[ordered, variable], climate.loc[ordered, measure])
            rho_d, p_d = stats.spearmanr(climate.loc[dune_sites, variable], climate.loc[dune_sites, measure])
            clim_rows.append(
                {
                    "climate_variable": variable,
                    "diversity_measure": measure,
                    "spearman_rho_all_60": rho_all,
                    "p_all_60": p_all,
                    "spearman_rho_dune_44": rho_d,
                    "p_dune_44": p_d,
                }
            )
    clim = pd.DataFrame(clim_rows)
    clim["q_bh_9_all_60"] = bh(clim["p_all_60"].to_numpy())
    clim["q_bh_9_dune_44"] = bh(clim["p_dune_44"].to_numpy())
    write_tsv(clim, output / "climate_diversity_dune_sensitivity.tsv")

    # ---- 2. compartment genus family --------------------------------------
    loadings = pd.read_csv(inputs["displacement_loadings"], sep="\t")
    family_rows = []
    for first, second in CONTRASTS:
        sites, tensor, n_blocks = ccr.block_centred_site_tensor(
            grouped, primary, (first, second), "pseudocount_0.5", None
        )
        differences = tensor[:, 0, :] - tensor[:, 1, :]
        mean_difference = differences.mean(axis=0)
        n_sites = differences.shape[0]
        signs = rng.choice((-1.0, 1.0), size=(args.permutations, n_sites))
        null = signs @ differences / n_sites
        p_values = (1 + (np.abs(null) >= np.abs(mean_difference)).sum(axis=0)) / (args.permutations + 1)
        idx = rng.integers(0, n_sites, size=(args.bootstrap, n_sites))
        boot = differences[idx].mean(axis=1)
        lo = np.quantile(boot, 0.025, axis=0)
        hi = np.quantile(boot, 0.975, axis=0)
        contrast = f"{first}-{second}"
        reference = loadings[loadings["contrast"] == contrast].set_index("genus")["mean_clr_difference"]
        for j, g in enumerate(primary):
            family_rows.append(
                {
                    "contrast": contrast,
                    "first_compartment": LABEL[first],
                    "second_compartment": LABEL[second],
                    "genus": g,
                    "phylum": composition["phylum"].get(g, ""),
                    "n_sites": n_sites,
                    "n_blocks": n_blocks,
                    "mean_clr_difference": float(mean_difference[j]),
                    "ci_low": float(lo[j]),
                    "ci_high": float(hi[j]),
                    "sign_flip_p": float(p_values[j]),
                    "committed_loading": float(reference.get(g, np.nan)),
                }
            )
    family = pd.DataFrame(family_rows)
    shared = family["committed_loading"].notna()
    max_dev = float((family.loc[shared, "mean_clr_difference"] - family.loc[shared, "committed_loading"]).abs().max())
    if max_dev > 1e-6:
        raise ValueError(f"Per-genus mean differences do not reproduce the committed loadings ({max_dev})")
    family["q_bh_600"] = bh(family["sign_flip_p"].to_numpy())
    family["supported_q_lt_0_05"] = family["q_bh_600"] < 0.05
    family["higher_in"] = np.where(family["mean_clr_difference"] > 0, family["first_compartment"], family["second_compartment"])
    family = family.sort_values(["contrast", "mean_clr_difference"], kind="mergesort")
    write_tsv(family, output / "compartment_genus_family.tsv")

    # ---- 3. compartment pathways ------------------------------------------
    pathways = pd.read_csv(inputs["pathway_position_effects"], sep="\t")
    supported = pathways[pathways["supported_q_lt_0_05"]].copy()
    supported["higher_in"] = np.where(
        supported["mean_clr_difference"] > 0,
        supported["contrast"].str.split("-").str[0],
        supported["contrast"].str.split("-").str[1],
    )
    top_rows = []
    for (contrast, higher), block in supported.groupby(["contrast", "higher_in"]):
        block = block.assign(magnitude=block["mean_clr_difference"].abs()).sort_values("magnitude", ascending=False, kind="mergesort").head(8)
        for rank, row in enumerate(block.itertuples(), start=1):
            top_rows.append(
                {
                    "contrast": contrast,
                    "higher_in": LABEL[higher],
                    "rank": rank,
                    "pathway": row.pathway,
                    "description": row.description,
                    "mean_clr_difference": row.mean_clr_difference,
                    "q_global_600": row.q_global_600,
                }
            )
    write_tsv(pd.DataFrame(top_rows), output / "compartment_pathway_leaders.tsv")

    # ---- 4. XRF axis -------------------------------------------------------
    axis = pd.read_csv(inputs["xrf_axis"], sep="\t")
    axis = axis[axis["Site"].between(1, 60)]
    axis_site = axis.groupby("Site")["elemental_pc1"].mean()
    axis_sites = [s for s in ordered if s in axis_site.index]
    xrf_rows = []
    for g in primary:
        rho, p = stats.spearmanr(axis_site.loc[axis_sites], site_means.loc[axis_sites, g])
        xrf_rows.append({"genus": g, "phylum": composition["phylum"].get(g, ""), "n_sites": len(axis_sites), "spearman_rho_elemental_axis": rho, "p_value": p})
    xrf_genus = pd.DataFrame(xrf_rows)
    xrf_genus["q_bh_200"] = bh(xrf_genus["p_value"].to_numpy())
    xrf_genus["supported_q_lt_0_05"] = xrf_genus["q_bh_200"] < 0.05
    xrf_genus = xrf_genus.sort_values("spearman_rho_elemental_axis", kind="mergesort")
    write_tsv(xrf_genus, output / "xrf_axis_genus_correlations.tsv")
    rho_axis_route, p_axis_route = stats.spearmanr([position_km[s] for s in axis_sites], axis_site.loc[axis_sites])
    xrf_loadings = pd.read_csv(inputs["xrf_loadings"], sep="\t")
    xrf_loadings["sign"] = np.where(xrf_loadings["pc1_loading"] > 0, "positive", "negative")
    write_tsv(xrf_loadings, output / "xrf_axis_loadings.tsv")

    # ---- 5. pH -----------------------------------------------------------------
    ph = pd.read_csv(inputs["ph_sample_profile_join"], sep="\t")
    ph = ph[ph["disposition"] == "ADMITTED_MEASUREMENT"]
    ph_site = ph.groupby("site")["ph_value"].mean()
    ph_sites = [s for s in ordered if s in ph_site.index]
    ph_rows = []
    for g in primary:
        rho, p = stats.spearmanr(ph_site.loc[ph_sites], site_means.loc[ph_sites, g])
        ph_rows.append({"genus": g, "phylum": composition["phylum"].get(g, ""), "n_sites": len(ph_sites), "spearman_rho_site_ph": rho, "p_value": p})
    ph_genus = pd.DataFrame(ph_rows)
    ph_genus["q_bh_200"] = bh(ph_genus["p_value"].to_numpy())
    ph_genus["supported_q_lt_0_05"] = ph_genus["q_bh_200"] < 0.05
    ph_genus["route_supported"] = ph_genus["genus"].map(ref["supported_q_lt_0_05"])
    ph_genus["route_rho"] = ph_genus["genus"].map(ref["spearman_rho_route_position"])
    ph_genus = ph_genus.sort_values("spearman_rho_site_ph", kind="mergesort")
    write_tsv(ph_genus, output / "ph_genus_correlations.tsv")

    # ---- 6. core genera per compartment ------------------------------------
    named = genus.loc[[g for g in genus.index if g != "unclassified_genus"]]
    core_rows = []
    occupancy_tables = {}
    for position in POSITIONS:
        ids = core.loc[core["position"] == position, ["sample_id", "site"]]
        pooled = named[ids["sample_id"].tolist()].T.groupby(ids["site"].to_numpy()).sum().T
        pooled = pooled.loc[:, pooled.sum(axis=0) >= CORE_DEPTH]
        matrix = pooled.to_numpy(dtype=np.int64)
        detected = np.zeros(matrix.shape, dtype=float)
        for _ in range(args.core_draws):
            for c in range(matrix.shape[1]):
                detected[:, c] += rng.multivariate_hypergeometric(matrix[:, c], CORE_DEPTH) > 0
        present = detected / args.core_draws >= 0.5
        occupancy = pd.Series(present.mean(axis=1), index=pooled.index)
        occupancy_tables[position] = occupancy
        core_set = occupancy[occupancy >= CORE_OCCUPANCY].sort_values(ascending=False, kind="mergesort")
        core_rows.append(
            {
                "compartment": LABEL[position],
                "n_sites_pooled": int(pooled.shape[1]),
                "n_genera_detected_any_site": int((occupancy > 0).sum()),
                "n_core_genera_ge_90pct_sites": int(len(core_set)),
                "mean_site_occupancy_of_detected_genera": float(occupancy[occupancy > 0].mean()),
                "median_rarefied_genera_per_site": float(np.median(present.sum(axis=0))),
                "core_genera": ",".join(core_set.index),
            }
        )
    core_table = pd.DataFrame(core_rows)
    occupancy_frame = pd.DataFrame({LABEL[p]: occupancy_tables[p] for p in POSITIONS}).fillna(0.0)
    occupancy_frame.index.name = "genus"
    occupancy_frame = occupancy_frame.reset_index()
    occupancy_frame["phylum"] = occupancy_frame["genus"].map(composition["phylum"]).fillna("")
    write_tsv(core_table, output / "core_genera_by_compartment.tsv")
    write_tsv(occupancy_frame.sort_values("root_adjacent", ascending=False, kind="mergesort"), output / "genus_site_occupancy_by_compartment.tsv")
    core_sets = {row.compartment: set(row.core_genera.split(",")) if row.core_genera else set() for row in core_table.itertuples()}
    shared_all = set.intersection(*core_sets.values())

    # ---- summary -------------------------------------------------------------
    supported_family = family[family["supported_q_lt_0_05"]]
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": "descriptive_biology_context_complete",
        "scope": (
            "Descriptive taxon-level context for tracked results; declared BH "
            "families (600 compartment genus tests; 200 tests each for landform, "
            "dune-only route, XRF axis and pH); no tracked conclusion changed."
        ),
        "landform": {
            "sites_per_landform": alpha_by_landform.set_index("landform")["n_sites"].to_dict(),
            "saline_pan_sites_by_third": site_table[site_table["landform"] == "saline pan"]["third"].value_counts().to_dict(),
            "omnibus": omnibus_rows,
            "supported_genera_unadjusted": int(per_genus[(per_genus["adjustment"] == "none") & per_genus["supported_q_lt_0_05"]].shape[0]),
            "supported_genera_route_adjusted": int(per_genus[(per_genus["adjustment"] == "route_linear_quadratic") & per_genus["supported_q_lt_0_05"]].shape[0]),
            "route_model": route_table.to_dict(orient="records"),
            "climate_diversity_dune_all_nine_supported": bool((clim["q_bh_9_dune_44"] < 0.05).all()),
        },
        "compartment_family": {
            "n_tests": int(len(family)),
            "n_supported": int(len(supported_family)),
            "supported_by_contrast_and_direction": {
                f"{c}:{h}": int(n) for (c, h), n in supported_family.groupby(["contrast", "higher_in"]).size().items()
            },
            "max_abs_deviation_from_committed_loadings": max_dev,
        },
        "xrf_axis": {
            "positive_loading_elements": xrf_loadings.loc[xrf_loadings["pc1_loading"] > 0.2, "element"].tolist(),
            "negative_loading_elements": xrf_loadings.loc[xrf_loadings["pc1_loading"] < -0.2, "element"].tolist(),
            "n_sites": len(axis_sites),
            "spearman_rho_axis_vs_route": float(rho_axis_route),
            "p_axis_vs_route": float(p_axis_route),
            "supported_genera": int(xrf_genus["supported_q_lt_0_05"].sum()),
        },
        "ph": {"n_sites": len(ph_sites), "supported_genera": int(ph_genus["supported_q_lt_0_05"].sum()),
               "supported_also_route_supported": int((ph_genus["supported_q_lt_0_05"] & ph_genus["route_supported"].fillna(False)).sum())},
        "core": {
            "depth": CORE_DEPTH,
            "occupancy_threshold": CORE_OCCUPANCY,
            "core_genera_per_compartment": {row.compartment: int(row.n_core_genera_ge_90pct_sites) for row in core_table.itertuples()},
            "core_shared_by_all_three": sorted(shared_all),
        },
        "parameters": {"seed": args.seed, "permutations": args.permutations, "bootstrap": args.bootstrap, "core_draws": args.core_draws},
        "inputs": {name: {"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha256(p)} for name, p in sorted(inputs.items())},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    (output / "run_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    def top(frame: pd.DataFrame, col: str, n: int, ascending: bool) -> str:
        block = frame[frame["supported_q_lt_0_05"]].sort_values(col, ascending=ascending, kind="mergesort").head(n)
        return ", ".join(f"{r.genus} ({getattr(r, col):+.2f})" for r in block.itertuples())

    readme = [
        "# Biology behind the statistics (descriptive context)",
        "",
        "- Status: `descriptive_biology_context_complete`",
        "",
        "## Landform",
        "",
        f"- Sites per landform: {summary['landform']['sites_per_landform']}",
        f"- Saline-pan sites by transect third: {summary['landform']['saline_pan_sites_by_third']}",
    ]
    for row in omnibus_rows:
        readme.append(f"- Dune vs saline pan composition, adjustment {row['adjustment']}: pseudo-F {row['pseudo_f']:.2f}, p {row['permutation_p']:.4f}")
    readme.append(f"- Genera differing (BH q<0.05): {summary['landform']['supported_genera_unadjusted']} unadjusted, {summary['landform']['supported_genera_route_adjusted']} route-adjusted")
    pan = per_genus[per_genus["adjustment"] == "route_linear_quadratic"]
    readme.append("- Higher in saline pans (route-adjusted): " + top(pan, "pan_minus_dune_clr", 10, False))
    readme.append("- Higher in dunes (route-adjusted): " + top(pan, "pan_minus_dune_clr", 10, True))
    for row in route_table.itertuples():
        readme.append(f"- Route model R2, {row.cohort}: {row.route_r2:.3f}")
    last = route_table.iloc[-1]
    readme.append(f"- Route genera supported: {int(last['supported_route_genera_all_sites'])} (all sites) vs {int(last['supported_route_genera_dune_sites'])} (dune sites); {int(last['supported_in_both_same_sign'])} supported in both with the same sign")
    readme.append(f"- Climate-diversity correlations on dune sites only: all nine supported = {summary['landform']['climate_diversity_dune_all_nine_supported']}")
    readme += ["", "## Compartment genus family (600 tests)", "", f"- Supported: {summary['compartment_family']['n_supported']} of 600; by contrast and direction {summary['compartment_family']['supported_by_contrast_and_direction']}"]
    for first, second in CONTRASTS:
        block = family[family["contrast"] == f"{first}-{second}"]
        readme.append(f"- {LABEL[first]} vs {LABEL[second]}: higher in {LABEL[first]}: " + top(block, "mean_clr_difference", 10, False))
        readme.append(f"  higher in {LABEL[second]}: " + top(block, "mean_clr_difference", 10, True))
    readme += ["", "## XRF elemental axis", "", f"- Positive loadings: {summary['xrf_axis']['positive_loading_elements']}; negative: {summary['xrf_axis']['negative_loading_elements']}",
               f"- Axis vs route position: rho {rho_axis_route:+.2f} (p {p_axis_route:.2g}); genera tracking the axis: {summary['xrf_axis']['supported_genera']}",
               "- Positive (evaporite/carbonate side): " + top(xrf_genus, "spearman_rho_elemental_axis", 10, False),
               "- Negative (quartz side): " + top(xrf_genus, "spearman_rho_elemental_axis", 10, True)]
    readme += ["", "## pH", "", f"- Sites {len(ph_sites)}; genera tracking site pH: {summary['ph']['supported_genera']} ({summary['ph']['supported_also_route_supported']} also route-supported)",
               "- Higher at higher pH: " + top(ph_genus, "spearman_rho_site_ph", 10, False),
               "- Higher at lower pH: " + top(ph_genus, "spearman_rho_site_ph", 10, True)]
    readme += ["", f"## Core genera per compartment (>= {CORE_OCCUPANCY:.0%} of sites after subsampling to {CORE_DEPTH} reads)", ""]
    for row in core_table.itertuples():
        readme.append(f"- {row.compartment}: {row.n_core_genera_ge_90pct_sites} core genera of {row.n_genera_detected_any_site} detected; mean occupancy {row.mean_site_occupancy_of_detected_genera:.3f}; median {row.median_rarefied_genera_per_site:.0f} genera per site")
    readme.append(f"- Core shared by all three compartments: {len(shared_all)}")
    readme += ["", "## Permitted wording", "",
               "- Landform contrasts are site-level and descriptive; saline pans sit at both ends of the route, so the route-adjusted contrast is the informative one.",
               "- Compartment genus contrasts are paired within site and campaign; they describe which genera carry the compartment difference and do not identify a mechanism.",
               "- The XRF axis is an elemental axis (Ca, Mg, Na, S, Cl, Fe, Ti positive; Si negative); do not call it salinity.",
               "- Core counts describe occupancy after subsampling and do not test a hypothesis.", ""]
    (output / "README.md").write_text("\n".join(readme), encoding="utf-8")
    lines = [f"{sha256(p)}  {p.name}" for p in sorted(output.iterdir()) if p.is_file() and p.name != "SHA256SUMS"]
    (output / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(readme))


if __name__ == "__main__":
    main()
