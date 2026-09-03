#!/usr/bin/env python3
"""Descriptive taxon and pathway context for the ecology manuscript.

Co-author review (R. Gruenberg, 2 Sep 2026) asked which taxa the Empty
Quarter communities contain, which taxa replace one another along the
transect, how the within-transect differences compare with the difference to
another desert, and which predicted metabolism dominates.  This module
answers those questions descriptively from the same canonical inputs that the
inferential analyses use.  It adds no new hypothesis test to the claim
ledger; the one family of correlation tests it reports (200 genera against
route position) is corrected with Benjamini--Hochberg and is descriptive of
the already supported geographic structure.

Analyses
--------
1. Phylum, class and genus composition of the 1,227 core-site profiles:
   mean relative abundance of total reads, prevalence, and the same
   quantities by compartment and by transect third.
2. Transect replacement: for the 200 primary genera in the same CLR
   coordinates as the geographic model (campaign-by-compartment means
   removed, site averages), the Spearman correlation with route position and
   the east-minus-west mean difference with a delete-one-site jackknife
   interval.
3. Genus-set overlap after rarefaction to a common depth: between the western,
   central and eastern thirds of the transect, and between the Empty Quarter
   and the public Atacama pit profiles (PRJEB39249) that the cross-desert
   module already uses.  Only presence at a fixed prevalence is compared;
   feature tables are not merged and different primers and DNA fractions
   remain a stated limit.
4. Predicted pathway dominance: the leading PICRUSt2 MetaCyc pathways by mean
   relative abundance, a keyword-based grouping into biosynthesis,
   degradation/utilization, energy/central metabolism and other, the rank of
   the autotrophic carbon-fixation and photosynthesis pathways, and the class
   composition of the pathways already supported in the route and compartment
   tests.

The taxonomy is SILVA 138.2 nomenclature as delivered by the canonical
workflow.  Relative abundances are shares of reads and are not cell counts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
from scipy import stats

SCHEMA_VERSION = "1.0"
POSITIONS = ("Surface", "Deep", "Rhizosphere")
POSITION_LABELS = {
    "Surface": "surface",
    "Deep": "shallow_subsurface",
    "Rhizosphere": "root_adjacent",
}
THIRDS = ("west", "central", "east")
RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
SAMPLE_RE = re.compile(
    r"^(?:e\d+_)?(?P<prefix>[TFSV])?(?P<site>\d+)"
    r"(?P<compartment>PR|P|D|S)r(?P<replicate>\d+)"
)
PREFIX_CAMPAIGN = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
CODE_COMPARTMENT = {"D": "Deep", "S": "Surface", "P": "Rhizosphere", "PR": "Rhizosphere"}

PATHWAY_CLASS_RULES = (
    (
        "energy_central_metabolism",
        re.compile(
            r"respiration|fermentation|glycolysis|gluconeogenesis|TCA cycle|"
            r"Calvin|photosynthesis|pentose phosphate|Entner|glyoxylate|"
            r"methanogenesis|oxidation|reductive|electron|acetyl-CoA|"
            r"superpathway of glucose|Bifidobacterium shunt|heterolactic|"
            r"incomplete reductive|methylaspartate cycle|glycolate",
            re.IGNORECASE,
        ),
    ),
    (
        "degradation_utilization",
        re.compile(
            r"degradation|utilization|assimilation|catabolism|cleavage|"
            r"Kdo|hydrolysis|conversion|interconversion",
            re.IGNORECASE,
        ),
    ),
    ("biosynthesis", re.compile(r"biosynthesis|salvage|elongation", re.IGNORECASE)),
)
AUTOTROPHY_PATHWAYS = {
    "CALVIN-PWY": "Calvin-Benson-Bassham cycle",
    "PWY-101": "photosynthesis light reactions",
    "P23-PWY": "reductive TCA cycle I",
    "PWY-5392": "reductive TCA cycle II",
    "CODH-PWY": "reductive acetyl coenzyme A pathway",
    "PWY-5743": "3-hydroxypropanoate cycle",
    "PWY-5744": "glyoxylate assimilation",
    "PWY-5789": "3-hydroxypropanoate/4-hydroxybutanate cycle",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_metadata(sample_id: str) -> dict[str, Any] | None:
    match = SAMPLE_RE.match(str(sample_id).replace(" ", ""))
    if match is None:
        return None
    prefix = match.group("prefix") or ""
    return {
        "sample_id": sample_id,
        "campaign": PREFIX_CAMPAIGN[prefix],
        "site": int(match.group("site")),
        "position": CODE_COMPARTMENT[match.group("compartment")],
        "replicate": int(match.group("replicate")),
    }


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    array = np.asarray(p_values, dtype=float)
    order = np.argsort(array)
    ranked = array[order] * array.size / (np.arange(array.size) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(array.size)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def parse_ranks(taxon: pd.Series) -> pd.DataFrame:
    """Split ``Kingdom;Phylum;...;Species;confidence`` strings into ranks."""
    splits = taxon.astype(str).str.split(";", expand=True)
    frame = pd.DataFrame(index=taxon.index)
    for index, rank in enumerate(RANKS):
        column = splits[index] if index < splits.shape[1] else pd.Series("", index=taxon.index)
        frame[rank] = column.fillna("").str.strip().replace({"NA": "", "nan": ""})
    return frame


def rank_taxa(counts: pd.DataFrame, prevalence_threshold: float) -> list[str]:
    """Rank prevalent taxa by mean relative abundance across groups."""
    prevalence = (counts > 0).mean(axis=0)
    relative = counts.div(counts.sum(axis=1), axis=0)
    eligible = prevalence[prevalence >= prevalence_threshold].index
    ranking = relative[eligible].mean(axis=0).sort_values(ascending=False, kind="mergesort")
    return ranking.index.tolist()


def clr(counts: pd.DataFrame, pseudocount: float) -> pd.DataFrame:
    logged = np.log(counts.to_numpy(dtype=float) + pseudocount)
    centred = logged - logged.mean(axis=1, keepdims=True)
    return pd.DataFrame(centred, index=counts.index, columns=counts.columns)


# --------------------------------------------------------------------------
# 1. Composition from the canonical feature table
# --------------------------------------------------------------------------
def accumulate_feature_table(
    feature_table: Path,
    taxonomy: Path,
    profiles: list[str],
    chunk_size: int,
) -> dict[str, pd.DataFrame]:
    """Sum ASV counts per phylum, class and genus for the requested profiles."""
    tax = pd.read_csv(taxonomy, sep="\t", index_col=0)
    ranks = parse_ranks(tax["Taxon"])
    level_maps = {
        "phylum": ranks["Phylum"],
        "class": ranks["Class"],
        "genus": ranks["Genus"],
    }
    genus_phylum = (
        ranks.loc[ranks["Genus"] != "", ["Genus", "Phylum", "Class"]]
        .groupby("Genus")
        .agg(lambda values: values.value_counts().index[0])
    )
    accumulators: dict[str, dict[str, np.ndarray]] = {key: {} for key in level_maps}
    depth = np.zeros(len(profiles))
    asv_counts: dict[str, dict[str, int]] = {key: {} for key in level_maps}
    n_asvs = 0
    for chunk in pd.read_csv(
        feature_table, sep="\t", index_col=0, skiprows=[0], chunksize=chunk_size
    ):
        if "Taxon" in chunk.columns:
            chunk = chunk.drop(columns=["Taxon"])
        chunk = chunk[profiles]
        values = chunk.to_numpy(dtype=np.float64)
        depth += values.sum(axis=0)
        n_asvs += len(chunk)
        for level, mapping in level_maps.items():
            labels = mapping.reindex(chunk.index).fillna("").to_numpy()
            labels = np.where(labels == "", f"unclassified_{level}", labels)
            frame = pd.DataFrame(values, index=labels)
            summed = frame.groupby(level=0, sort=False).sum()
            for label, row in zip(summed.index, summed.to_numpy()):
                acc = accumulators[level]
                if label in acc:
                    acc[label] += row
                else:
                    acc[label] = row.copy()
            for label, count in pd.Series(labels).value_counts().items():
                asv_counts[level][label] = asv_counts[level].get(label, 0) + int(count)
    tables = {
        level: pd.DataFrame(acc, index=profiles).T for level, acc in accumulators.items()
    }
    tables["depth"] = pd.DataFrame({"reads": depth}, index=profiles)
    tables["genus_phylum"] = genus_phylum
    for level in level_maps:
        tables[f"{level}_asv_counts"] = pd.Series(asv_counts[level], name="n_asvs")
    tables["n_asvs_total"] = pd.DataFrame({"n": [n_asvs]})
    return tables


def composition_table(
    counts: pd.DataFrame,
    depth: pd.Series,
    metadata: pd.DataFrame,
    asv_counts: pd.Series,
    level: str,
) -> pd.DataFrame:
    """Mean share of total reads and prevalence, overall and by stratum."""
    relative = counts.div(depth.reindex(counts.columns), axis=1)
    rows = []
    for taxon in relative.index:
        row = {
            level: taxon,
            "n_asvs": int(asv_counts.get(taxon, 0)),
            "mean_relative_abundance": float(relative.loc[taxon].mean()),
            "median_relative_abundance": float(relative.loc[taxon].median()),
            "prevalence": float((counts.loc[taxon] > 0).mean()),
        }
        for position in POSITIONS:
            ids = metadata.loc[metadata["position"] == position, "sample_id"]
            row[f"mean_{POSITION_LABELS[position]}"] = float(relative.loc[taxon, ids].mean())
        for third in THIRDS:
            ids = metadata.loc[metadata["third"] == third, "sample_id"]
            row[f"mean_{third}_third"] = float(relative.loc[taxon, ids].mean())
            row[f"prevalence_{third}_third"] = float((counts.loc[taxon, ids] > 0).mean())
        rows.append(row)
    frame = pd.DataFrame(rows).sort_values(
        "mean_relative_abundance", ascending=False, kind="mergesort"
    )
    frame.insert(0, "rank", np.arange(1, len(frame) + 1))
    return frame.reset_index(drop=True)


# --------------------------------------------------------------------------
# 2. Transect replacement in the geographic-model coordinates
# --------------------------------------------------------------------------
def grouped_counts(
    genus: pd.DataFrame, metadata: pd.DataFrame, minimum_group_reads: int
) -> pd.DataFrame:
    values = genus[metadata["sample_id"].tolist()].T
    values["campaign"] = metadata["campaign"].to_numpy()
    values["site"] = metadata["site"].to_numpy()
    values["position"] = metadata["position"].to_numpy()
    grouped = values.groupby(["campaign", "site", "position"], sort=True).sum(numeric_only=True)
    return grouped.loc[grouped.sum(axis=1) >= minimum_group_reads]


def transect_replacement(
    grouped: pd.DataFrame,
    taxa: list[str],
    coordinates: pd.DataFrame,
    thirds: dict[int, str],
    pseudocount: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Site-level CLR profiles with campaign-by-position means removed."""
    transformed = clr(grouped[taxa], pseudocount)
    frame = transformed.reset_index()
    centred = frame.copy()
    centred[taxa] = frame[taxa] - frame.groupby(["campaign", "position"])[taxa].transform("mean")
    site_means = centred.groupby("site")[taxa].mean()
    site_means = site_means.reindex(sorted(site_means.index))
    position = coordinates.set_index("site").loc[site_means.index, "transect_km"].to_numpy()
    third_labels = np.asarray([thirds[int(site)] for site in site_means.index])
    west = third_labels == "west"
    east = third_labels == "east"
    rows = []
    for genus_name in taxa:
        y = site_means[genus_name].to_numpy()
        rho, p_value = stats.spearmanr(position, y)
        difference = float(y[east].mean() - y[west].mean())
        # Delete-one-site jackknife over the 40 sites that enter the contrast.
        pseudo = []
        for drop in np.flatnonzero(west | east):
            keep = np.ones(y.size, dtype=bool)
            keep[drop] = False
            pseudo.append(y[east & keep].mean() - y[west & keep].mean())
        pseudo = np.asarray(pseudo)
        n_units = pseudo.size
        standard_error = float(np.sqrt((n_units - 1) / n_units * np.sum((pseudo - pseudo.mean()) ** 2)))
        t_value = stats.t.ppf(0.975, df=n_units - 1)
        rows.append(
            {
                "genus": genus_name,
                "n_sites": int(y.size),
                "spearman_rho_route_position": float(rho),
                "p_value": float(p_value),
                "east_minus_west_mean_clr": difference,
                "east_minus_west_ci_low": difference - t_value * standard_error,
                "east_minus_west_ci_high": difference + t_value * standard_error,
                "mean_clr_west_third": float(y[west].mean()),
                "mean_clr_central_third": float(y[third_labels == "central"].mean()),
                "mean_clr_east_third": float(y[east].mean()),
            }
        )
    result = pd.DataFrame(rows)
    result["q_bh_200"] = benjamini_hochberg(result["p_value"].to_numpy())
    result["supported_q_lt_0_05"] = result["q_bh_200"] < 0.05
    result["direction"] = np.where(
        result["spearman_rho_route_position"] > 0, "increases_eastward", "decreases_eastward"
    )
    result = result.sort_values("spearman_rho_route_position", kind="mergesort").reset_index(drop=True)
    return result, site_means


# --------------------------------------------------------------------------
# 3. Genus-set overlap after rarefaction
# --------------------------------------------------------------------------
def rarefied_prevalence(
    counts: pd.DataFrame, depth: int, draws: int, rng: np.random.Generator
) -> tuple[pd.Series, int]:
    """Mean fraction of profiles in which each taxon is detected after
    subsampling every profile (columns) to ``depth`` reads, over ``draws``
    independent subsamples.  Profiles below ``depth`` are excluded."""
    matrix = counts.to_numpy(dtype=np.int64)
    totals = matrix.sum(axis=0)
    keep = totals >= depth
    matrix = matrix[:, keep]
    if matrix.shape[1] == 0:
        raise ValueError("No profile reaches the rarefaction depth")
    detected = np.zeros(matrix.shape[0])
    for _ in range(draws):
        for column in range(matrix.shape[1]):
            sample = rng.multivariate_hypergeometric(matrix[:, column], depth)
            detected += sample > 0
    prevalence = detected / (draws * matrix.shape[1])
    return pd.Series(prevalence, index=counts.index), int(matrix.shape[1])


def overlap_row(
    label_a: str,
    set_a: set[str],
    label_b: str,
    set_b: set[str],
    top_a: list[str],
    top_b: list[str],
    n_profiles_a: int,
    n_profiles_b: int,
    note: str,
) -> dict[str, Any]:
    shared = set_a & set_b
    union = set_a | set_b
    return {
        "set_a": label_a,
        "set_b": label_b,
        "n_profiles_a": n_profiles_a,
        "n_profiles_b": n_profiles_b,
        "n_genera_a": len(set_a),
        "n_genera_b": len(set_b),
        "n_shared": len(shared),
        "jaccard": len(shared) / len(union) if union else float("nan"),
        "fraction_of_a_in_b": len(shared) / len(set_a) if set_a else float("nan"),
        "fraction_of_b_in_a": len(shared) / len(set_b) if set_b else float("nan"),
        "top50_a_detected_in_b": int(sum(1 for g in top_a[:50] if g in set_b)),
        "top50_b_detected_in_a": int(sum(1 for g in top_b[:50] if g in set_a)),
        "note": note,
    }


# --------------------------------------------------------------------------
# 4. Predicted pathway dominance
# --------------------------------------------------------------------------
def classify_pathway(description: str) -> str:
    for label, pattern in PATHWAY_CLASS_RULES:
        if pattern.search(description):
            return label
    return "other"


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.10g", lineterminator="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=20260903)
    parser.add_argument("--rarefaction-depth", type=int, default=8000)
    parser.add_argument("--rarefaction-draws", type=int, default=20)
    parser.add_argument("--detection-prevalence", type=float, default=0.10)
    parser.add_argument("--minimum-group-reads", type=int, default=2000)
    parser.add_argument("--group-prevalence", type=float, default=0.20)
    parser.add_argument("--top-genera", type=int, default=200)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--chunk-size", type=int, default=10000)
    args = parser.parse_args()

    root = args.project_root.resolve()
    output = (args.output_dir or root / "analysis/v3/taxon_context").resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    inputs = {
        "genus_counts": root / "analysis/v2/review/cache/genus_counts.tsv",
        "feature_table": root / "data/metadata/taxonomy/feature-table-trips1-5.tsv",
        "taxonomy": root / "data/metadata/taxonomy/taxonomy-trips1-5.tsv",
        "site_coordinates": root / "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv",
        "picrust_pathways": root / "data/processed/functional/picrust2/merged/path_abun_unstrat.tsv",
        "picrust_descriptions": root / "data/processed/functional/picrust2/path_abun_unstrat_descriptions.tsv",
        "picrust_position_effects": root / "analysis/v3/picrust2_ecology/pathway_position_effects.tsv",
        "picrust_geographic_correlations": root / "analysis/v3/picrust2_ecology/pathway_geographic_correlations.tsv",
        "atacama_pit_asv_table": root / "data/metadata/comparators/atacama/pit/ASV_table.tsv",
        "atacama_pit_taxonomy": root / "data/metadata/comparators/atacama/pit/ASV_tax.silva_138_2.tsv",
        "atacama_pit_depth_map": root / "data/metadata/comparators/atacama/pit/sample_depth_map.tsv",
        "climate_site_summary": root / "analysis/v3/environment_associations/climate_site_summary.tsv",
        "ph_sample_profile_join": root / "analysis/v3/ph_shared_v1/ecology/ph_sample_profile_join.tsv",
    }
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing inputs:\n" + "\n".join(missing))

    # ---- profiles and strata -------------------------------------------
    genus_cache = pd.read_csv(inputs["genus_counts"], sep="\t", index_col=0)
    genus_cache = genus_cache.loc[~genus_cache.index.isna()]
    genus_cache.index = genus_cache.index.astype(str)
    parsed = [m for sid in genus_cache.columns if (m := sample_metadata(sid)) is not None]
    metadata = pd.DataFrame(parsed)
    if len(metadata) != genus_cache.shape[1]:
        raise ValueError("Some quality-controlled profile identifiers did not parse")
    core = metadata[metadata["site"].between(1, 60) & metadata["position"].isin(POSITIONS)].copy()
    coordinates = pd.read_csv(inputs["site_coordinates"], sep="\t")
    coordinates = coordinates[coordinates["site"].between(1, 60)].sort_values("transect_km")
    if len(coordinates) != 60:
        raise ValueError("Expected 60 core-site coordinates")
    ordered_sites = coordinates["site"].astype(int).tolist()
    thirds = {}
    for index, site in enumerate(ordered_sites):
        thirds[site] = THIRDS[min(index // 20, 2)]
    core["third"] = core["site"].map(thirds)
    core_ids = core["sample_id"].tolist()

    # ---- composition from the feature table ----------------------------
    tables = accumulate_feature_table(
        inputs["feature_table"], inputs["taxonomy"], core_ids, args.chunk_size
    )
    depth = tables["depth"]["reads"]
    genus_from_table = tables["genus"].drop(index="unclassified_genus", errors="ignore")
    shared_genera = genus_from_table.index.intersection(genus_cache.index)
    cache_core = genus_cache.loc[shared_genera, core_ids]
    max_abs_difference = float(
        (genus_from_table.loc[shared_genera, core_ids] - cache_core).abs().to_numpy().max()
    )
    if max_abs_difference > 0.5 or len(shared_genera) != len(genus_cache.index):
        raise ValueError(
            "Genus sums from the canonical feature table do not reproduce the cached "
            f"genus table (max |difference| {max_abs_difference}, "
            f"{len(shared_genera)} of {len(genus_cache.index)} genera matched)"
        )
    phylum_table = composition_table(
        tables["phylum"], depth, core, tables["phylum_asv_counts"], "phylum"
    )
    class_table = composition_table(
        tables["class"], depth, core, tables["class_asv_counts"], "class"
    )
    genus_table = composition_table(
        tables["genus"], depth, core, tables["genus_asv_counts"], "genus"
    )
    genus_phylum = tables["genus_phylum"]
    genus_table.insert(2, "phylum", genus_table["genus"].map(genus_phylum["Phylum"]).fillna(""))
    genus_table.insert(3, "class", genus_table["genus"].map(genus_phylum["Class"]).fillna(""))
    write_tsv(phylum_table, output / "phylum_composition.tsv")
    write_tsv(class_table, output / "class_composition.tsv")
    write_tsv(genus_table, output / "genus_composition.tsv")

    # Leading genera per compartment and per transect third.
    named = genus_table[genus_table["genus"] != "unclassified_genus"]
    leaders = []
    for column, label in (
        ("mean_relative_abundance", "all_core_profiles"),
        ("mean_surface", "surface"),
        ("mean_shallow_subsurface", "shallow_subsurface"),
        ("mean_root_adjacent", "root_adjacent"),
        ("mean_west_third", "west_third"),
        ("mean_central_third", "central_third"),
        ("mean_east_third", "east_third"),
    ):
        top = named.sort_values(column, ascending=False, kind="mergesort").head(10)
        for position, (_, row) in enumerate(top.iterrows(), start=1):
            leaders.append(
                {
                    "stratum": label,
                    "rank_in_stratum": position,
                    "genus": row["genus"],
                    "phylum": row["phylum"],
                    "mean_relative_abundance": float(row[column]),
                }
            )
    write_tsv(pd.DataFrame(leaders), output / "leading_genera_by_stratum.tsv")

    # ---- transect replacement -------------------------------------------
    grouped = grouped_counts(genus_cache, core, args.minimum_group_reads)
    ranked = rank_taxa(grouped, args.group_prevalence)
    primary = ranked[: args.top_genera]
    replacement, site_means = transect_replacement(
        grouped, primary, coordinates, thirds, args.pseudocount
    )
    replacement.insert(1, "phylum", replacement["genus"].map(genus_phylum["Phylum"]).fillna(""))
    share = genus_table.set_index("genus")
    for third in THIRDS:
        replacement[f"mean_relative_abundance_{third}_third"] = replacement["genus"].map(
            share[f"mean_{third}_third"]
        )
    write_tsv(replacement, output / "transect_replacement.tsv")
    site_means.round(6).reset_index().to_csv(
        output / "site_clr_means_200_genera.tsv", sep="\t", index=False, lineterminator="\n"
    )

    # ---- genus-set overlap after rarefaction ----------------------------
    eq_prevalence: dict[str, pd.Series] = {}
    eq_profiles: dict[str, int] = {}
    eq_named = genus_cache.loc[[g for g in genus_cache.index if g != "unclassified_genus"]]
    for third in THIRDS:
        ids = core.loc[core["third"] == third, "sample_id"].tolist()
        eq_prevalence[third], eq_profiles[third] = rarefied_prevalence(
            eq_named[ids], args.rarefaction_depth, args.rarefaction_draws, rng
        )
    eq_prevalence["all"], eq_profiles["all"] = rarefied_prevalence(
        eq_named[core_ids], args.rarefaction_depth, args.rarefaction_draws, rng
    )
    surface_ids = core.loc[core["position"] == "Surface", "sample_id"].tolist()
    eq_prevalence["surface"], eq_profiles["surface"] = rarefied_prevalence(
        eq_named[surface_ids], args.rarefaction_depth, args.rarefaction_draws, rng
    )

    pit_counts = pd.read_csv(inputs["atacama_pit_asv_table"], sep="\t", index_col=0)
    pit_tax = pd.read_csv(inputs["atacama_pit_taxonomy"], sep="\t", index_col=0)
    pit_depths = pd.read_csv(inputs["atacama_pit_depth_map"], sep="\t")
    prokaryote = pit_tax.index[pit_tax["Kingdom"].isin(["Bacteria", "Archaea"])]
    pit_counts = pit_counts.loc[pit_counts.index.intersection(prokaryote)]
    pit_genus_labels = pit_tax.loc[pit_counts.index, "Genus"].fillna("").astype(str).str.strip()
    pit_named = pit_counts.loc[pit_genus_labels != ""].groupby(
        pit_genus_labels[pit_genus_labels != ""].to_numpy()
    ).sum()
    pit_named = pit_named[pit_depths["sampleID"].tolist()]
    pit_prevalence, pit_profiles = rarefied_prevalence(
        pit_named, args.rarefaction_depth, args.rarefaction_draws, rng
    )
    upper_ids = pit_depths.loc[pit_depths["depth_cm"] <= 10, "sampleID"].tolist()
    pit_upper_prevalence, pit_upper_profiles = rarefied_prevalence(
        pit_named[upper_ids], args.rarefaction_depth, args.rarefaction_draws, rng
    )

    def detected(prevalence: pd.Series) -> set[str]:
        return set(prevalence.index[prevalence >= args.detection_prevalence])

    def leading(counts: pd.DataFrame) -> list[str]:
        relative = counts.div(counts.sum(axis=0), axis=1)
        return relative.mean(axis=1).sort_values(ascending=False, kind="mergesort").index.tolist()

    eq_top = {
        third: leading(eq_named[core.loc[core["third"] == third, "sample_id"].tolist()])
        for third in THIRDS
    }
    eq_top["all"] = leading(eq_named[core_ids])
    eq_top["surface"] = leading(eq_named[surface_ids])
    pit_top = leading(pit_named)
    pit_upper_top = leading(pit_named[upper_ids])

    within_note = (
        "Same primers, workflow and taxonomy; genera detected in at least "
        f"{args.detection_prevalence:.0%} of profiles after {args.rarefaction_draws} "
        f"subsamples to {args.rarefaction_depth} reads."
    )
    between_note = (
        "Different primers (V3-V4 versus V4), DNA fraction (total versus "
        "intracellular), depth range and site design; same SILVA 138.2 "
        "nomenclature and the same detection rule. Presence comparison only."
    )
    overlap_rows = [
        overlap_row(
            "Empty Quarter west third", detected(eq_prevalence["west"]),
            "Empty Quarter east third", detected(eq_prevalence["east"]),
            eq_top["west"], eq_top["east"], eq_profiles["west"], eq_profiles["east"], within_note,
        ),
        overlap_row(
            "Empty Quarter west third", detected(eq_prevalence["west"]),
            "Empty Quarter central third", detected(eq_prevalence["central"]),
            eq_top["west"], eq_top["central"], eq_profiles["west"], eq_profiles["central"], within_note,
        ),
        overlap_row(
            "Empty Quarter central third", detected(eq_prevalence["central"]),
            "Empty Quarter east third", detected(eq_prevalence["east"]),
            eq_top["central"], eq_top["east"], eq_profiles["central"], eq_profiles["east"], within_note,
        ),
        overlap_row(
            "Empty Quarter all core profiles", detected(eq_prevalence["all"]),
            "Atacama pit all depths (PRJEB39249)", detected(pit_prevalence),
            eq_top["all"], pit_top, eq_profiles["all"], pit_profiles, between_note,
        ),
        overlap_row(
            "Empty Quarter surface", detected(eq_prevalence["surface"]),
            "Atacama pit 2.5-10 cm (PRJEB39249)", detected(pit_upper_prevalence),
            eq_top["surface"], pit_upper_top, eq_profiles["surface"], pit_upper_profiles, between_note,
        ),
    ]
    write_tsv(pd.DataFrame(overlap_rows), output / "genus_set_overlap.tsv")

    # ---- direction of the environmental gradients along the route ---------
    climate = pd.read_csv(inputs["climate_site_summary"], sep="\t")
    ph = pd.read_csv(inputs["ph_sample_profile_join"], sep="\t")
    ph = ph[ph["disposition"] == "ADMITTED_MEASUREMENT"]
    ph_site = ph.groupby("site")["ph_value"].mean().rename("mean_ph").reset_index()
    gradient_frame = coordinates[["site", "transect_km"]].merge(
        climate[["site", "mean_air_temperature_c", "mean_monthly_rain_mm", "mean_relative_humidity_pct"]],
        on="site",
        how="left",
    ).merge(ph_site, on="site", how="left")
    gradient_rows = []
    for column, label in (
        ("mean_air_temperature_c", "49-month mean air temperature"),
        ("mean_monthly_rain_mm", "49-month mean monthly rainfall"),
        ("mean_relative_humidity_pct", "49-month mean relative humidity"),
        ("mean_ph", "archived-soil pH (admitted measurements, site mean)"),
    ):
        block = gradient_frame.dropna(subset=[column])
        rho, p_value = stats.spearmanr(block["transect_km"], block[column])
        west_sites = block[block["site"].map(thirds) == "west"]
        east_sites = block[block["site"].map(thirds) == "east"]
        gradient_rows.append(
            {
                "variable": column,
                "description": label,
                "n_sites": int(len(block)),
                "spearman_rho_route_position": float(rho),
                "p_value": float(p_value),
                "mean_west_third": float(west_sites[column].mean()),
                "mean_east_third": float(east_sites[column].mean()),
            }
        )
    site_gradients = pd.DataFrame(gradient_rows)
    write_tsv(site_gradients, output / "site_gradients.tsv")

    # ---- predicted pathway dominance -------------------------------------
    pathways = pd.read_csv(inputs["picrust_pathways"], sep="\t", index_col=0)
    descriptions = pd.read_csv(
        inputs["picrust_descriptions"], sep="\t", usecols=["pathway", "description"]
    ).set_index("pathway")["description"]
    pathway_profiles = [sid for sid in core_ids if sid in pathways.columns]
    pathways = pathways[pathway_profiles]
    relative_pathways = pathways.div(pathways.sum(axis=0), axis=1)
    dominance = pd.DataFrame(
        {
            "pathway": relative_pathways.index,
            "description": descriptions.reindex(relative_pathways.index).fillna("").to_numpy(),
            "mean_relative_abundance": relative_pathways.mean(axis=1).to_numpy(),
            "prevalence": (pathways > 0).mean(axis=1).to_numpy(),
        }
    )
    dominance["pathway_class"] = dominance["description"].map(classify_pathway)
    dominance["autotrophy_marker"] = dominance["pathway"].isin(AUTOTROPHY_PATHWAYS)
    dominance = dominance.sort_values(
        "mean_relative_abundance", ascending=False, kind="mergesort"
    ).reset_index(drop=True)
    dominance.insert(0, "rank", np.arange(1, len(dominance) + 1))
    for position in POSITIONS:
        ids = core.loc[core["position"] == position, "sample_id"]
        ids = [sid for sid in ids if sid in relative_pathways.columns]
        dominance[f"mean_{POSITION_LABELS[position]}"] = dominance["pathway"].map(
            relative_pathways[ids].mean(axis=1)
        )
    write_tsv(dominance, output / "pathway_dominance.tsv")

    class_share = (
        dominance.groupby("pathway_class")["mean_relative_abundance"]
        .agg(["sum", "count"])
        .rename(columns={"sum": "share_of_predicted_pathway_abundance", "count": "n_pathways"})
        .reset_index()
        .sort_values("share_of_predicted_pathway_abundance", ascending=False, kind="mergesort")
    )
    position_effects = pd.read_csv(inputs["picrust_position_effects"], sep="\t")
    geographic = pd.read_csv(inputs["picrust_geographic_correlations"], sep="\t")
    position_effects["pathway_class"] = position_effects["description"].map(classify_pathway)
    geographic["pathway_class"] = geographic["description"].map(classify_pathway)
    supported_rows = []
    for label, frame, sign_column in (
        ("route_correlation_supported", geographic, "spearman_rho"),
        ("compartment_contrast_supported", position_effects, "mean_clr_difference"),
    ):
        supported = frame[frame["supported_q_lt_0_05"]]
        for pathway_class, block in supported.groupby("pathway_class"):
            supported_rows.append(
                {
                    "family": label,
                    "pathway_class": pathway_class,
                    "n_supported": int(len(block)),
                    "n_positive": int((block[sign_column] > 0).sum()),
                    "n_negative": int((block[sign_column] < 0).sum()),
                }
            )
    write_tsv(class_share, output / "pathway_class_share.tsv")
    write_tsv(pd.DataFrame(supported_rows), output / "supported_pathways_by_class.tsv")

    autotrophy = dominance[dominance["autotrophy_marker"]][
        ["rank", "pathway", "description", "mean_relative_abundance", "prevalence"]
    ]
    write_tsv(autotrophy, output / "autotrophy_pathways.tsv")

    # ---- manifest, README, checksums ------------------------------------
    cyanobacteria = phylum_table.loc[phylum_table["phylum"] == "Cyanobacteriota"]
    top_phyla = phylum_table.head(8)
    replacement_supported = replacement[replacement["supported_q_lt_0_05"]]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "scope": (
            "Descriptive taxon and predicted-pathway context; presence-based "
            "genus-set overlap only; no feature-table merging across surveys "
            "and no mechanism claim."
        ),
        "status": "descriptive_context_complete",
        "cohort": {
            "quality_controlled_profiles": int(genus_cache.shape[1]),
            "core_site_profiles": int(len(core_ids)),
            "core_sites": 60,
            "asvs_in_feature_table": int(tables["n_asvs_total"].iloc[0, 0]),
            "genera_in_cache": int(genus_cache.shape[0]),
            "genus_sum_max_abs_difference_vs_cache": max_abs_difference,
            "picrust_profiles": int(len(pathway_profiles)),
            "pathways": int(len(dominance)),
            "grouped_profiles": int(len(grouped)),
            "primary_genera": int(len(primary)),
            "transect_thirds": {
                third: [site for site in ordered_sites if thirds[site] == third] for third in THIRDS
            },
            "rarefied_profiles": {**eq_profiles, "atacama_pit_all": pit_profiles, "atacama_pit_upper": pit_upper_profiles},
        },
        "parameters": {
            "seed": args.seed,
            "rarefaction_depth": args.rarefaction_depth,
            "rarefaction_draws": args.rarefaction_draws,
            "detection_prevalence": args.detection_prevalence,
            "minimum_group_reads": args.minimum_group_reads,
            "group_prevalence": args.group_prevalence,
            "top_genera": args.top_genera,
            "pseudocount": args.pseudocount,
            "pathway_class_rules": {
                label: pattern.pattern for label, pattern in PATHWAY_CLASS_RULES
            },
        },
        "headline": {
            "top_phyla": [
                {"phylum": row.phylum, "mean_relative_abundance": row.mean_relative_abundance}
                for row in top_phyla.itertuples()
            ],
            "cyanobacteriota_mean_relative_abundance": (
                float(cyanobacteria["mean_relative_abundance"].iloc[0]) if len(cyanobacteria) else None
            ),
            "leading_genera_all": leaders[:5],
            "route_supported_genera": int(len(replacement_supported)),
            "route_increasing_eastward": int((replacement_supported["direction"] == "increases_eastward").sum()),
            "route_decreasing_eastward": int((replacement_supported["direction"] == "decreases_eastward").sum()),
            "genus_set_overlap": overlap_rows,
            "site_gradients": site_gradients.to_dict(orient="records"),
            "pathway_class_share": class_share.to_dict(orient="records"),
            "autotrophy_pathways": autotrophy.to_dict(orient="records"),
        },
        "inputs": {
            name: {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for name, path in sorted(inputs.items())
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    readme = [
        "# Taxon and predicted-pathway context",
        "",
        "- Status: `descriptive_context_complete`",
        f"- Cohort: {len(core_ids)} core-site profiles (sites 1-60), {manifest['cohort']['asvs_in_feature_table']} ASVs, "
        f"{genus_cache.shape[0]} named genera",
        f"- Genus sums recomputed from the canonical feature table reproduce `analysis/v2/review/cache/genus_counts.tsv` "
        f"(max |difference| {max_abs_difference:g} reads)",
        "",
        "## Leading phyla (mean share of total reads)",
        "",
    ]
    for row in top_phyla.itertuples():
        readme.append(f"- {row.phylum}: {100 * row.mean_relative_abundance:.1f} % (prevalence {100 * row.prevalence:.0f} %)")
    readme += ["", "## Leading genera (all core profiles)", ""]
    for entry in leaders[:10]:
        readme.append(f"- {entry['genus']} ({entry['phylum']}): {100 * entry['mean_relative_abundance']:.2f} %")
    readme += [
        "",
        "## Transect replacement (200 primary genera, site-level CLR)",
        "",
        f"- Supported route correlations (BH q < 0.05 over 200 tests): {len(replacement_supported)} "
        f"({manifest['headline']['route_increasing_eastward']} increase eastward, "
        f"{manifest['headline']['route_decreasing_eastward']} decrease eastward)",
        "- Strongest eastward decreases: "
        + ", ".join(replacement.head(5)["genus"]),
        "- Strongest eastward increases: "
        + ", ".join(replacement.tail(5)["genus"][::-1]),
        "",
        "## Genus-set overlap (rarefied presence)",
        "",
    ]
    for row in overlap_rows:
        readme.append(
            f"- {row['set_a']} vs {row['set_b']}: {row['n_genera_a']} vs {row['n_genera_b']} genera, "
            f"{row['n_shared']} shared, Jaccard {row['jaccard']:.3f}; "
            f"{row['top50_a_detected_in_b']}/50 leading genera of A detected in B"
        )
    readme += ["", "## Environmental gradients along the route (Spearman rho with route position)", ""]
    for row in site_gradients.itertuples():
        readme.append(f"- {row.description}: rho = {row.spearman_rho_route_position:+.2f} (n = {row.n_sites}; west third {row.mean_west_third:.2f}, east third {row.mean_east_third:.2f})")
    readme += ["", "## Predicted pathway classes (share of predicted pathway abundance)", ""]
    for row in class_share.itertuples():
        readme.append(f"- {row.pathway_class}: {100 * row.share_of_predicted_pathway_abundance:.1f} % ({row.n_pathways} pathways)")
    readme += ["", "## Autotrophy-related pathways", ""]
    for row in autotrophy.itertuples():
        readme.append(f"- rank {row.rank}: {row.pathway} {row.description}: {100 * row.mean_relative_abundance:.2f} %")
    readme += [
        "",
        "## Permitted wording",
        "",
        "- Relative abundances are shares of reads and describe the marker-gene profile, not cell counts.",
        "- The pathway classes are keyword groupings of MetaCyc descriptions defined in `run_manifest.json`; they are not MetaCyc ontology classes.",
        "- The Atacama comparison is a presence-based overlap under a shared detection rule; it does not merge feature tables and does not correct for primer or DNA-fraction differences.",
        "- Route correlations describe which genera change along the transect; they do not identify a cause.",
        "",
    ]
    (output / "README.md").write_text("\n".join(readme), encoding="utf-8")

    checksum_lines = []
    for path in sorted(output.iterdir()):
        if path.name == "SHA256SUMS" or not path.is_file():
            continue
        checksum_lines.append(f"{sha256(path)}  {path.name}")
    (output / "SHA256SUMS").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print("\n".join(readme))


if __name__ == "__main__":
    main()
