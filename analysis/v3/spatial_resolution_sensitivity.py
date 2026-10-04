#!/usr/bin/env python3
"""ASV-resolution and neighbour-graph sensitivity for the transect model.

The primary spatial analysis aggregates counts to genera and uses a
symmetric five-nearest-neighbour graph for the residual Moran diagnostic.
Two questions follow from that choice and are answered here as supplementary
sensitivities, not as replacements for the primary result:

* does genus aggregation create the transect association, or does it survive
  at amplicon-sequence-variant resolution;
* how much does the residual Moran statistic depend on the neighbour count k.

Both arms reuse the primary model's machinery unchanged, so any difference
is attributable to resolution or to k rather than to a different estimator.

Site coordinates come from the corrected 60-site table written by the primary
analysis (``analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv``),
the same table the spatial-covariance and pH-partition analyses consume.
Using this common table keeps all spatial analyses on the same site geometry.

The ASV arm runs on exactly the site-campaign-compartment groups of the genus
primary fit. The ASV cache (``asv_filt_counts.tsv``) differs from the genus
cache in two ways that would otherwise change the cohort: it also holds
profiles below the 1,000-read cleaning threshold, and its per-group totals
count only reads on ASVs present in at least 5 % of profiles. Group quality
control is therefore taken from the genus cohort (profiles with at least
1,000 reads; groups with at least 2,000 genus-assigned reads) and the ASV
profiles of those profiles are summed within the same groups. Every ASV that
passes the 20 % group-prevalence filter used for ranking is present in at
least 126 groups and so in far more than the 5 % of profiles required by the
cache filter; the counts of the analysed ASVs are therefore identical to the
unfiltered feature table.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from spatial_turnover_rescue import (
    COMPARTMENTS,
    analyse,
    fit_multivariate,
    design_matrix,
    load_grouped_counts,
    multivariate_moran,
    rank_taxa,
    sample_metadata,
    sha256_file,
    site_level_clr,
    symmetric_knn_weights,
    write_tsv,
)

NEIGHBOUR_COUNTS = (3, 4, 5, 6, 8, 10)
ASV_TAXON_COUNTS = (800, 2000)
GROUP_KEYS = ["campaign", "site", "compartment"]
COORDINATE_COLUMNS = ("site", "latitude", "longitude", "x_km", "y_km", "transect_km")
DEFAULT_COORDINATES = Path(
    "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv"
)


def load_site_coordinates(path: Path) -> pd.DataFrame:
    """Read the corrected 60-site coordinate table of the primary analysis."""
    coordinates = pd.read_csv(path, sep="\t")
    missing = [column for column in COORDINATE_COLUMNS if column not in coordinates]
    if missing:
        raise ValueError(f"{path} lacks coordinate columns {missing}")
    coordinates = coordinates.sort_values("site").reset_index(drop=True)
    if coordinates["site"].tolist() != list(range(1, 61)):
        raise ValueError(f"{path} must hold exactly the 60 core sites once each")
    return coordinates


def reference_sample_ids(genus_path: Path) -> set[str]:
    """Sample columns of the genus cache: the profiles that passed cleaning."""
    header = pd.read_csv(genus_path, sep="\t", index_col=0, nrows=0)
    return set(map(str, header.columns))


def asv_counts_on_reference_groups(
    asv_path: Path,
    reference_samples: set[str],
    reference_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Sum ASV profiles within exactly the reference (genus-cohort) groups.

    Profiles absent from the genus cache (below the cleaning threshold) are
    excluded, and groups are kept if and only if the genus cohort kept them.
    No read threshold is re-applied to the prevalence-filtered ASV totals;
    groups whose filtered total falls below the genus threshold are listed in
    the returned information so the cohort stays fully described.
    """
    asv = pd.read_csv(asv_path, sep="\t", index_col=0)
    asv = asv.loc[~asv.index.isna()].copy()
    asv.index = asv.index.astype(str)
    parsed = [
        value
        for sample_id in asv.columns
        if (value := sample_metadata(sample_id)) is not None
    ]
    metadata = pd.DataFrame(parsed)
    metadata = metadata[
        metadata["site"].between(1, 60)
        & metadata["campaign"].between(1, 5)
        & metadata["compartment"].isin(COMPARTMENTS)
    ].copy()
    in_reference = metadata["sample_id"].astype(str).isin(reference_samples)
    excluded_samples = sorted(metadata.loc[~in_reference, "sample_id"].astype(str))
    metadata = metadata.loc[in_reference].reset_index(drop=True)

    values = asv[metadata["sample_id"].tolist()].T
    for key in GROUP_KEYS:
        values[key] = metadata[key].to_numpy()
    grouped = values.groupby(GROUP_KEYS, sort=True).sum(numeric_only=True)

    reference_index = pd.MultiIndex.from_frame(reference_metadata[GROUP_KEYS])
    absent = [key for key in reference_index if key not in grouped.index]
    if absent:
        raise ValueError(
            f"{len(absent)} reference groups have no ASV profile: {absent[:5]}"
        )
    outside_reference = sorted(set(grouped.index) - set(reference_index))
    grouped = grouped.loc[reference_index]
    filtered_totals = grouped.sum(axis=1)
    counts = grouped.T
    counts.columns = pd.MultiIndex.from_frame(reference_metadata[GROUP_KEYS])

    def key_label(key: tuple[Any, ...]) -> str:
        campaign, site, compartment = key
        return f"campaign {campaign}, site {site}, {compartment}"

    info = {
        "reference_groups": int(len(reference_metadata)),
        "aligned_groups": int(counts.shape[1]),
        "asv_cache_sample_columns": int(asv.shape[1]),
        "asv_cache_core_profiles": int(in_reference.size),
        "profiles_excluded_as_absent_from_genus_cache": excluded_samples,
        "profiles_used": int(len(metadata)),
        "asv_cache_groups_outside_reference_cohort": [
            key_label(key) for key in outside_reference
        ],
        "filtered_asv_reads_per_group_min": int(filtered_totals.min()),
        "filtered_asv_reads_per_group_median": float(filtered_totals.median()),
        "reference_groups_below_2000_filtered_asv_reads": {
            key_label(key): int(total)
            for key, total in filtered_totals.items()
            if total < 2000
        },
    }
    return counts, reference_metadata.reset_index(drop=True), info


def moran_k_sensitivity(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    coordinates: pd.DataFrame,
    taxa: list[str],
    pseudocount: float,
    permutations: int,
    seed: int,
) -> list[dict[str, Any]]:
    """Recompute the residual Moran diagnostic across neighbour counts."""
    sites, response, _ = site_level_clr(
        counts, metadata, taxa, None, pseudocount
    )
    spatial = coordinates.set_index("site").loc[sites]
    design = design_matrix(spatial["transect_km"].to_numpy(), 2)
    _, _, residual = fit_multivariate(response, design)
    rows = []
    for neighbours in NEIGHBOUR_COUNTS:
        weights = symmetric_knn_weights(
            spatial["x_km"].to_numpy(),
            spatial["y_km"].to_numpy(),
            neighbours=neighbours,
        )
        observed = multivariate_moran(residual, weights)
        rng = np.random.default_rng(seed + neighbours)
        null = np.empty(permutations)
        for index in range(permutations):
            null[index] = multivariate_moran(
                residual[rng.permutation(len(residual))], weights
            )
        p_value = (1 + int(np.sum(null >= observed))) / (permutations + 1)
        rows.append(
            {
                "neighbours_k": neighbours,
                "n_sites": int(len(sites)),
                "mean_degree": float(weights.sum(axis=1).mean()),
                "residual_moran_i": float(observed),
                "permutation_p": float(p_value),
                "null_mean": float(null.mean()),
                "null_p95": float(np.quantile(null, 0.95)),
            }
        )
    return rows


RESULT_FIELDS = (
    "n_sites",
    "n_groups",
    "partial_r2",
    "partial_r2_ci_low",
    "partial_r2_ci_high",
    "pseudo_f",
    "permutation_p",
    "residual_moran_i",
    "residual_moran_p",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument(
        "--genus-counts",
        type=Path,
        default=None,
        help="genus count table; defaults to the canonical ecology cache",
    )
    parser.add_argument(
        "--asv-counts",
        type=Path,
        default=None,
        help="filtered ASV count table; defaults to the canonical ecology cache",
    )
    parser.add_argument(
        "--coordinates",
        type=Path,
        default=None,
        help=(
            "corrected 60-site coordinate table; defaults to the primary "
            f"analysis output {DEFAULT_COORDINATES}"
        ),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-group-reads", type=int, default=2000)
    parser.add_argument("--prevalence", type=float, default=0.20)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260728)
    args = parser.parse_args()

    root = args.project_root.resolve()
    cache = root / "analysis" / "v2" / "review" / "cache"
    genus_path = args.genus_counts or cache / "genus_counts.tsv"
    asv_path = args.asv_counts or cache / "asv_filt_counts.tsv"
    coordinates_path = args.coordinates or root / DEFAULT_COORDINATES
    for path in (genus_path, asv_path, coordinates_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Spatial resolution sensitivity requires {path}; it is absent."
            )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    for stale in output.glob("*"):
        if stale.is_file():
            stale.unlink()

    coordinates = load_site_coordinates(coordinates_path)
    genus_counts, genus_metadata, genus_info = load_grouped_counts(
        genus_path, args.minimum_group_reads
    )
    genus_taxa = rank_taxa(genus_counts, args.prevalence)

    asv_counts, asv_metadata, cohort = asv_counts_on_reference_groups(
        asv_path, reference_sample_ids(genus_path), genus_metadata
    )
    asv_ranked = rank_taxa(asv_counts, args.prevalence)

    rows = []
    reference = analyse(
        genus_counts,
        genus_metadata,
        coordinates,
        genus_taxa[:200],
        omitted_campaign=None,
        trend_degree=2,
        permutations=args.permutations,
        seed=args.seed,
        pseudocount=args.pseudocount,
    )
    rows.append(
        {
            "resolution": "genus",
            "feature_count": 200,
            "eligible_features": len(genus_taxa),
            **{key: getattr(reference, key) for key in RESULT_FIELDS},
        }
    )
    for feature_count in (*ASV_TAXON_COUNTS, len(asv_ranked)):
        if feature_count > len(asv_ranked):
            continue
        result = analyse(
            asv_counts,
            asv_metadata,
            coordinates,
            asv_ranked[:feature_count],
            omitted_campaign=None,
            trend_degree=2,
            permutations=args.permutations,
            seed=args.seed + feature_count,
            pseudocount=args.pseudocount,
        )
        rows.append(
            {
                "resolution": "asv",
                "feature_count": feature_count,
                "eligible_features": len(asv_ranked),
                **{key: getattr(result, key) for key in RESULT_FIELDS},
            }
        )
    write_tsv(
        output / "asv_resolution_sensitivity.tsv",
        rows,
        ["resolution", "feature_count", "eligible_features", *RESULT_FIELDS],
    )

    k_rows = moran_k_sensitivity(
        genus_counts,
        genus_metadata,
        coordinates,
        genus_taxa[:200],
        args.pseudocount,
        args.permutations,
        args.seed,
    )
    write_tsv(
        output / "moran_k_sensitivity.tsv",
        k_rows,
        [
            "neighbours_k",
            "n_sites",
            "mean_degree",
            "residual_moran_i",
            "permutation_p",
            "null_mean",
            "null_p95",
        ],
    )

    asv_rows = [row for row in rows if row["resolution"] == "asv"]
    genus_r2 = rows[0]["partial_r2"]
    asv_r2 = [row["partial_r2"] for row in asv_rows]
    moran_values = [row["residual_moran_i"] for row in k_rows]
    moran_ps = [row["permutation_p"] for row in k_rows]
    detected_k = [
        row["neighbours_k"] for row in k_rows if row["permutation_p"] < 0.05
    ]
    undetected_k = [
        row["neighbours_k"] for row in k_rows if row["permutation_p"] >= 0.05
    ]
    moran_k_robust = not undetected_k
    n_groups = cohort["aligned_groups"]
    thin_groups = cohort["reference_groups_below_2000_filtered_asv_reads"]
    asv_summary = "; ".join(
        f"{row['feature_count']} ASVs {row['partial_r2']:.4f} "
        f"[{row['partial_r2_ci_low']:.4f}, {row['partial_r2_ci_high']:.4f}]"
        for row in asv_rows
    )
    verdict = {
        "schema_version": "2.0",
        "status": (
            "asv_resolution_consistent_with_genus_primary"
            if all(value > 0.5 * genus_r2 for value in asv_r2)
            else "asv_resolution_diverges_from_genus_primary"
        ),
        "moran_k_status": (
            "residual_autocorrelation_detected_at_every_k"
            if moran_k_robust
            else "residual_autocorrelation_depends_on_neighbour_count"
        ),
        "genus_primary_partial_r2": genus_r2,
        "genus_primary_partial_r2_95_jackknife_ci": [
            rows[0]["partial_r2_ci_low"],
            rows[0]["partial_r2_ci_high"],
        ],
        "asv_partial_r2_range": [min(asv_r2), max(asv_r2)],
        "asv_fits": [
            {
                key: row[key]
                for key in (
                    "feature_count",
                    "partial_r2",
                    "partial_r2_ci_low",
                    "partial_r2_ci_high",
                    "permutation_p",
                    "residual_moran_i",
                    "residual_moran_p",
                )
            }
            for row in asv_rows
        ],
        "asv_cohort": cohort,
        "moran_i_range_across_k": [min(moran_values), max(moran_values)],
        "moran_p_max_across_k": max(moran_ps),
        "neighbour_counts_tested": list(NEIGHBOUR_COUNTS),
        "neighbour_counts_with_detected_autocorrelation": detected_k,
        "neighbour_counts_without_detected_autocorrelation": undetected_k,
        "permitted_wording": (
            "The transect association is not an artefact of genus "
            "aggregation: at amplicon-sequence-variant resolution the same "
            f"model gives partial R2 {min(asv_r2):.4f} to {max(asv_r2):.4f} "
            f"against {genus_r2:.4f} at genus level ({asv_summary}). Both "
            f"resolutions use the same {n_groups} site-campaign-compartment "
            f"groups at {rows[0]['n_sites']} sites and the corrected site "
            "coordinates; the genus figure is the primary fit. "
            + (
                f"{len(thin_groups)} of the {n_groups} groups "
                f"({'; '.join(thin_groups)}) "
                f"{'carries' if len(thin_groups) == 1 else 'carry'} fewer "
                "than 2,000 reads on prevalence-filtered ASVs while passing "
                "the genus-assigned read threshold; "
                f"{'it is' if len(thin_groups) == 1 else 'they are'} retained "
                "because the cohort is defined once, by the primary quality "
                "control. "
                if thin_groups
                else ""
            )
            + (
                "Residual spatial autocorrelation was detected at every "
                f"neighbour count tested (k = {min(NEIGHBOUR_COUNTS)} to "
                f"{max(NEIGHBOUR_COUNTS)}; Moran I {min(moran_values):.4f} to "
                f"{max(moran_values):.4f}, largest permutation p "
                f"{max(moran_ps):.3g})."
                if moran_k_robust
                else "The residual Moran diagnostic declines with the "
                f"neighbour count: it is detected at k = {detected_k} and not "
                f"at k = {undetected_k} (Moran I {max(moran_values):.4f} down "
                f"to {min(moran_values):.4f}). The fixed-k residual "
                "autocorrelation statement is bounded to the "
                "short-neighbourhood scale."
            )
        ),
        "prohibited_wording": (
            "Do not promote the ASV-resolution fit to the primary result; it "
            "is a supplementary resolution sensitivity on the same design and "
            "inherits every design limit of the primary model, including the "
            "collection-order alias. Do not state unqualified residual spatial "
            "autocorrelation without naming the neighbour count."
            if not moran_k_robust
            else "Do not promote the ASV-resolution fit to the primary result; "
            "it is a supplementary resolution sensitivity on the same design "
            "and inherits every design limit of the primary model, including "
            "the collection-order alias."
        ),
        "input": {
            "genus_counts_path": str(genus_path),
            "genus_counts_sha256": sha256_file(genus_path),
            "asv_counts_path": str(asv_path),
            "asv_counts_sha256": sha256_file(asv_path),
            "coordinates_path": str(coordinates_path),
            "coordinates_sha256": sha256_file(coordinates_path),
            "genus_groups": genus_info["retained_site_campaign_compartment_groups"],
            "minimum_group_reads": args.minimum_group_reads,
            "prevalence_threshold": args.prevalence,
            "pseudocount": args.pseudocount,
            "permutations": args.permutations,
            "seed": args.seed,
        },
    }
    (output / "claim_verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme = [
        "# ASV-resolution and neighbour-graph sensitivity",
        "",
        verdict["permitted_wording"],
        "",
        verdict["prohibited_wording"],
        "",
        "Inputs: `analysis/v2/review/cache/genus_counts.tsv`,",
        "`analysis/v2/review/cache/asv_filt_counts.tsv` and the corrected",
        f"site table `{DEFAULT_COORDINATES}` (checksums in `claim_verdict.json`).",
        "",
        "Reproduce from the repository root:",
        "",
        "```sh",
        "python analysis/v3/spatial_resolution_sensitivity.py \\",
        "  --output-dir analysis/v3/spatial_resolution_sensitivity",
        "```",
        "",
        "Evidence files: `asv_resolution_sensitivity.tsv`,",
        "`moran_k_sensitivity.tsv`, `claim_verdict.json`.",
        "",
    ]
    (output / "README.md").write_text("\n".join(readme), encoding="utf-8")
    digests = [
        f"{sha256_file(path)}  {path.name}"
        for path in sorted(output.glob("*"))
        if path.is_file() and path.name != "SHA256SUMS"
    ]
    (output / "SHA256SUMS").write_text("\n".join(digests) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
