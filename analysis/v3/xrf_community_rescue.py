#!/usr/bin/env python3
"""Re-analyse laboratory XRF/community associations without calling XRF salinity.

The analysis uses all 725 available laboratory records (547 from Trips 1--4
and 178 from Trip 5), keeps field XRF out of the ecological join, and excludes
Trip-1-only sites 61--64 from repeated-campaign inference. Instrument zeros are
not interpreted as measured absence: the primary PCA is limited to elements
reported positive in at least 75% of laboratory records.

The primary elemental axis is calculated after within-campaign standardisation,
so a Trip-5 method shift cannot itself define the axis. Alpha-diversity
associations use site fixed effects and site-clustered errors. Community
composition uses a Bray-Curtis PCoA followed by a nested multivariate model;
the elemental-axis residual is permuted within site after removing campaign,
compartment, and site effects.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
import statsmodels
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA


SAMPLE_RE = re.compile(
    r"^(?:e\d+_)?(?P<prefix>[TFSV])?(?P<site>\d+)"
    r"(?P<compartment>PR|P|D|S)r(?P<rep>\d+)"
)
PREFIX_TRIP = {"": 1, "T": 2, "F": 3, "S": 4, "V": 5}
COMPARTMENT = {
    "D": "Deep",
    "S": "Surface",
    "P": "Rhizosphere",
    "PR": "Rhizosphere",
}
META_COLUMNS = {
    "SampleID",
    "SoilType",
    "Material",
    "Mode",
    "Diameter",
    "Method",
}
ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")
CORE_SITES = set(range(1, 61))
SEED = 20260723
LAB_XRF_TRIPS1_4 = "data/processed/geochemistry/xrf_lab_table_trips1-4.tsv"
LAB_XRF_TRIP5 = "data/processed/geochemistry/xrf_lab_table_filtered.tsv"
ALPHA_CACHE = "analysis/v2/review/cache/alpha.tsv"
GENUS_CACHE = "analysis/v2/review/cache/genus_counts.tsv"
GEODATA = [f"data/metadata/geodata/trip{trip}_geodata.tsv" for trip in range(1, 6)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def data_repository_lock(root: Path) -> dict[str, str]:
    """Return the pinned data-repository record, if the lock file exists."""
    path = root / "DATA_REPOSITORY.lock"
    if not path.is_file():
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(frame["field"], frame["value"]))


def parse_sample_id(
    sample_id: str,
    default_trip: int | None = None,
) -> tuple[int, int, str] | None:
    """Return the campaign, site, and compartment encoded by a sample ID."""
    match = SAMPLE_RE.match(str(sample_id).replace(" ", ""))
    if not match:
        return None
    prefix = match.group("prefix") or ""
    return (
        default_trip if default_trip is not None else PREFIX_TRIP[prefix],
        int(match.group("site")),
        COMPARTMENT[match.group("compartment")],
    )


def load_xrf_table(path: Path, default_trip: int | None) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    id_column = "SampleID" if "SampleID" in frame else frame.columns[0]
    keys = frame[id_column].map(
        lambda value: parse_sample_id(value, default_trip=default_trip)
    )
    frame = frame.loc[keys.notna()].copy()
    parsed = list(keys[keys.notna()])
    frame["Trip"] = [item[0] for item in parsed]
    frame["Site"] = [item[1] for item in parsed]
    frame["Type"] = [item[2] for item in parsed]
    return frame


def load_all_lab_xrf(root: Path) -> tuple[pd.DataFrame, list[str]]:
    t14 = load_xrf_table(root / LAB_XRF_TRIPS1_4, default_trip=None)
    t5 = load_xrf_table(root / LAB_XRF_TRIP5, default_trip=5)
    shared_elements = [
        column
        for column in t14.columns
        if (
            column in t5
            and column not in META_COLUMNS
            and column != "LE"
            and ELEMENT_RE.fullmatch(column)
        )
    ]
    columns = ["Trip", "Site", "Type", *shared_elements]
    frame = pd.concat([t14[columns], t5[columns]], ignore_index=True)
    for column in shared_elements:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = (
        frame.groupby(["Trip", "Site", "Type"], as_index=False)[shared_elements]
        .mean()
        .sort_values(["Trip", "Site", "Type"])
    )
    return frame, shared_elements


def fit_elemental_axis(
    frame: pd.DataFrame,
    candidates: list[str],
    detection_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int | str]]:
    detection = (frame[candidates].fillna(0) > 0).mean()
    retained = detection[detection >= detection_threshold].index.tolist()
    if len(retained) < 3:
        raise ValueError("Fewer than three well-detected XRF elements")

    transformed = np.log1p(frame[retained])
    within_trip = transformed.groupby(frame["Trip"]).transform(
        lambda values: (values - values.mean()) / values.std(ddof=0)
    )
    within_trip = within_trip.replace([np.inf, -np.inf], np.nan).fillna(0)
    model = PCA(n_components=min(len(retained), len(frame)))
    scores = model.fit_transform(within_trip)

    evaporite_elements = [
        element for element in ("Cl", "Na", "S", "Sr", "Br") if element in retained
    ]
    evaporite_reference = within_trip[evaporite_elements].mean(axis=1)
    orientation = float(
        np.corrcoef(scores[:, 0], evaporite_reference.to_numpy())[0, 1]
    )
    sign = -1.0 if np.isfinite(orientation) and orientation < 0 else 1.0

    result = frame[["Trip", "Site", "Type"]].copy()
    result["elemental_pc1"] = sign * scores[:, 0]

    global_z = (transformed - transformed.mean()) / transformed.std(ddof=0)
    global_z = global_z.replace([np.inf, -np.inf], np.nan).fillna(0)
    global_scores = PCA(n_components=1).fit_transform(global_z)[:, 0]
    if np.corrcoef(result["elemental_pc1"], global_scores)[0, 1] < 0:
        global_scores *= -1
    result["elemental_pc1_global_sensitivity"] = global_scores

    loadings = pd.DataFrame(
        {
            "element": retained,
            "positive_detection_fraction": detection[retained].to_numpy(),
            "pc1_loading": sign * model.components_[0],
        }
    ).sort_values("pc1_loading", ascending=False)
    info: dict[str, float | int | str] = {
        "laboratory_records": len(frame),
        "retained_elements": ",".join(retained),
        "retained_element_count": len(retained),
        "detection_threshold": detection_threshold,
        "pc1_variance_explained": float(model.explained_variance_ratio_[0]),
        "pc1_evaporite_reference_correlation": float(abs(orientation)),
        "within_vs_global_pc1_correlation": float(
            np.corrcoef(
                result["elemental_pc1"],
                result["elemental_pc1_global_sensitivity"],
            )[0, 1]
        ),
    }
    return result, loadings, info


def load_coordinates(root: Path) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for trip, relative in enumerate(GEODATA, start=1):
        frame = pd.read_csv(root / relative, sep="\t")
        frame["Site"] = pd.to_numeric(frame["Site"], errors="coerce")
        frame = frame.dropna(subset=["Site", "Latitude", "Longitude"])
        for _, row in frame.iterrows():
            rows.append(
                {
                    "Trip": trip,
                    "Site": int(row["Site"]),
                    "Latitude": float(row["Latitude"]),
                    "Longitude": float(row["Longitude"]),
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["Trip", "Site"])


def load_community(
    root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    alpha = pd.read_csv(root / ALPHA_CACHE, sep="\t", index_col=0)
    counts = pd.read_csv(root / GENUS_CACHE, sep="\t", index_col=0)
    shared = counts.columns.intersection(alpha.index)
    counts = counts[shared]
    metadata = alpha.loc[shared, ["Trip", "Site", "Type", "shannon", "depth"]].copy()
    metadata["Trip"] = metadata["Trip"].astype(int)
    metadata["Site"] = metadata["Site"].astype(int)

    keyed = counts.T.join(metadata[["Trip", "Site", "Type"]])
    community = keyed.groupby(["Trip", "Site", "Type"]).sum(numeric_only=True)
    community = community.loc[community.sum(axis=1) >= 2000]
    community = community.loc[
        community.index.get_level_values("Site").isin(CORE_SITES)
    ]
    alpha_aggregate = (
        metadata.groupby(["Trip", "Site", "Type"], as_index=False)
        .agg(
            shannon=("shannon", "mean"),
            sequencing_depth=("depth", "sum"),
        )
    )
    alpha_aggregate = alpha_aggregate[
        alpha_aggregate["Site"].isin(CORE_SITES)
    ]
    return community, alpha_aggregate


def pcoa_coordinates(
    distance: np.ndarray,
    variance_fraction: float,
) -> tuple[np.ndarray, float]:
    """Return positive-axis PCoA coordinates up to a variance threshold."""
    n = distance.shape[0]
    center = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * center @ np.square(distance) @ center
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > max(1e-12, eigenvalues[0] * 1e-10)
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]
    cumulative = np.cumsum(eigenvalues) / eigenvalues.sum()
    axes = int(np.searchsorted(cumulative, variance_fraction) + 1)
    coordinates = eigenvectors[:, :axes] * np.sqrt(eigenvalues[:axes])
    return coordinates, float(cumulative[axes - 1])


def design_matrix(frame: pd.DataFrame) -> np.ndarray:
    categorical = frame[["Trip", "Type", "Site"]].astype(str)
    encoded = pd.get_dummies(categorical, drop_first=True, dtype=float)
    return np.column_stack([np.ones(len(frame)), encoded.to_numpy()])


def residual_sse(design: np.ndarray, response: np.ndarray) -> float:
    coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
    residual = response - design @ coefficients
    return float(np.square(residual).sum())


def nested_multivariate_test(
    response: np.ndarray,
    reduced: np.ndarray,
    predictor: np.ndarray,
    sites: np.ndarray,
    permutations: int,
    seed: int = SEED,
) -> dict[str, float | int]:
    """Test one predictor after a reduced model using within-site permutations."""
    predictor = np.asarray(predictor, dtype=float)
    predictor_fit = reduced @ np.linalg.lstsq(
        reduced, predictor, rcond=None
    )[0]
    predictor_residual = predictor - predictor_fit
    full = np.column_stack([reduced, predictor])
    sse_reduced = residual_sse(reduced, response)
    sse_full = residual_sse(full, response)
    effect_ss = max(0.0, sse_reduced - sse_full)
    residual_df = len(predictor) - np.linalg.matrix_rank(full)
    observed_f = effect_ss / (sse_full / residual_df)

    groups = {
        site: np.flatnonzero(sites == site) for site in np.unique(sites)
    }
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        permuted = predictor_residual.copy()
        for indices in groups.values():
            permuted[indices] = rng.permutation(permuted[indices])
        candidate = predictor_fit + permuted
        candidate_full = np.column_stack([reduced, candidate])
        candidate_sse = residual_sse(candidate_full, response)
        candidate_effect = max(0.0, sse_reduced - candidate_sse)
        candidate_f = candidate_effect / (candidate_sse / residual_df)
        exceed += candidate_f >= observed_f
    return {
        "pseudo_f": float(observed_f),
        "p": (exceed + 1) / (permutations + 1),
        "partial_r2": effect_ss / sse_reduced,
        "residual_df": int(residual_df),
        "permutations": permutations,
    }


def analyse(
    root: Path,
    output: Path,
    permutations: int,
    detection_threshold: float,
) -> dict[str, object]:
    xrf, elements = load_all_lab_xrf(root)
    if len(xrf) != 725:
        raise ValueError(f"Expected 725 laboratory XRF records; found {len(xrf)}")
    axis, loadings, axis_info = fit_elemental_axis(
        xrf,
        elements,
        detection_threshold,
    )
    community, alpha = load_community(root)
    coordinates = load_coordinates(root)

    analysis = (
        axis.merge(alpha, on=["Trip", "Site", "Type"], how="inner")
        .merge(coordinates, on=["Trip", "Site"], how="left")
        .sort_values(["Trip", "Site", "Type"])
    )
    analysis["log_sequencing_depth"] = np.log1p(
        analysis["sequencing_depth"]
    )
    fixed = smf.ols(
        "shannon ~ elemental_pc1 + log_sequencing_depth + "
        "C(Trip) + C(Type) + C(Site)",
        data=analysis,
    ).fit(cov_type="cluster", cov_kwds={"groups": analysis["Site"]})
    spatial = smf.ols(
        "shannon ~ elemental_pc1 + log_sequencing_depth + C(Trip) + "
        "C(Type) + Latitude + Longitude + I(Latitude ** 2) + "
        "I(Longitude ** 2)",
        data=analysis,
    ).fit(cov_type="cluster", cov_kwds={"groups": analysis["Site"]})
    global_sensitivity = smf.ols(
        "shannon ~ elemental_pc1_global_sensitivity + "
        "log_sequencing_depth + C(Trip) + C(Type) + C(Site)",
        data=analysis,
    ).fit(cov_type="cluster", cov_kwds={"groups": analysis["Site"]})

    alpha_rows = []
    for model_name, model, term in (
        ("site_fixed_primary", fixed, "elemental_pc1"),
        ("geographic_trend_sensitivity", spatial, "elemental_pc1"),
        (
            "global_standardisation_sensitivity",
            global_sensitivity,
            "elemental_pc1_global_sensitivity",
        ),
    ):
        interval = model.conf_int().loc[term]
        alpha_rows.append(
            {
                "model": model_name,
                "n_observations": int(model.nobs),
                "n_sites": analysis["Site"].nunique(),
                "estimate": float(model.params[term]),
                "std_error": float(model.bse[term]),
                "ci_low": float(interval.iloc[0]),
                "ci_high": float(interval.iloc[1]),
                "p": float(model.pvalues[term]),
            }
        )
    alpha_results = pd.DataFrame(alpha_rows)

    community_keys = pd.DataFrame(
        community.index.tolist(),
        columns=["Trip", "Site", "Type"],
    )
    community_keys["row_index"] = np.arange(len(community_keys))
    joined_keys = (
        axis.merge(community_keys, on=["Trip", "Site", "Type"], how="inner")
        .sort_values(["Trip", "Site", "Type"])
        .reset_index(drop=True)
    )
    counts = community.iloc[joined_keys["row_index"].to_numpy()]
    relative = counts.div(counts.sum(axis=1), axis=0)
    distance = squareform(pdist(relative.to_numpy(), metric="braycurtis"))
    reduced = design_matrix(joined_keys)
    multivariate_rows = []
    for fraction in (0.80, 0.95):
        pcoa, retained = pcoa_coordinates(distance, fraction)
        result = nested_multivariate_test(
            pcoa,
            reduced,
            joined_keys["elemental_pc1"].to_numpy(),
            joined_keys["Site"].to_numpy(),
            permutations=permutations,
        )
        multivariate_rows.append(
            {
                "pcoa_target_fraction": fraction,
                "pcoa_retained_fraction": retained,
                "pcoa_axes": pcoa.shape[1],
                "n_observations": len(joined_keys),
                "n_sites": joined_keys["Site"].nunique(),
                **result,
            }
        )
    multivariate = pd.DataFrame(multivariate_rows)

    alpha_primary = alpha_results[
        alpha_results["model"] == "site_fixed_primary"
    ].iloc[0]
    beta_primary = multivariate[
        multivariate["pcoa_target_fraction"] == 0.80
    ].iloc[0]
    alpha_supported = bool(alpha_primary["p"] < 0.05)
    beta_supported = bool(
        beta_primary["p"] < 0.05 and (multivariate["p"] < 0.05).all()
    )
    if alpha_supported and beta_supported:
        status = "elemental_axis_association"
    elif beta_supported:
        status = "composition_only_elemental_axis_association"
    elif alpha_supported:
        status = "alpha_only_elemental_axis_association"
    else:
        status = "not_supported"

    output.mkdir(parents=True, exist_ok=True)
    loadings.to_csv(output / "elemental_pc1_loadings.tsv", sep="\t", index=False)
    axis.to_csv(output / "laboratory_xrf_axis.tsv", sep="\t", index=False)
    analysis.to_csv(
        output / "xrf_alpha_analysis_table.tsv",
        sep="\t",
        index=False,
    )
    alpha_results.to_csv(
        output / "xrf_alpha_models.tsv",
        sep="\t",
        index=False,
    )
    multivariate.to_csv(
        output / "xrf_community_models.tsv",
        sep="\t",
        index=False,
    )
    summary: dict[str, object] = {
        "schema_version": "1.0",
        "status": status,
        "permitted_wording": (
            "Laboratory XRF defines an elemental/mineralogical gradient "
            "associated with community structure after campaign, compartment, "
            "and site adjustment."
            if beta_supported
            else "No robust adjusted community association with the primary "
            "laboratory-XRF elemental axis was detected."
        ),
        "prohibited_wording": (
            "Do not call the axis salinity or infer a salinity mechanism until "
            "independent EC/soluble-ion calibration is available."
        ),
        "counts": {
            "laboratory_xrf_records": len(xrf),
            "trips1_4_records": int((xrf["Trip"] <= 4).sum()),
            "trip5_records": int((xrf["Trip"] == 5).sum()),
            "alpha_joined_observations": len(analysis),
            "community_joined_observations": len(joined_keys),
            "core_sites": joined_keys["Site"].nunique(),
        },
        "axis": axis_info,
        "alpha_primary": alpha_primary.to_dict(),
        "community_primary": beta_primary.to_dict(),
        "field_xrf_role": (
            "Excluded from the ecological join; retained as a separate "
            "Trip-5 site-level diagnostic because no physical aliquot ID exists."
        ),
    }
    (output / "xrf_claim_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    readme = [
        "# Laboratory XRF/community claim rescue",
        "",
        f"- Status: `{status}`",
        f"- Laboratory records: {len(xrf)} "
        f"({summary['counts']['trips1_4_records']} Trips 1–4; "
        f"{summary['counts']['trip5_records']} Trip 5)",
        f"- Community joins: {len(joined_keys)} core-site observations",
        f"- Primary alpha p: {alpha_primary['p']:.4g}",
        f"- Primary multivariate p: {beta_primary['p']:.4g}",
        f"- Primary multivariate partial R2: {beta_primary['partial_r2']:.4g}",
        "",
        str(summary["permitted_wording"]),
        "",
        str(summary["prohibited_wording"]),
        "",
        "Field XRF is not pooled with laboratory XRF. Sites 61–64 remain in "
        "the data release but are excluded from repeated-campaign inference.",
        "",
    ]
    (output / "README.md").write_text("\n".join(readme))
    write_run_manifest(
        root,
        output,
        {
            "permutations": permutations,
            "detection_threshold": detection_threshold,
            "seed": SEED,
        },
    )
    return summary


def write_run_manifest(
    root: Path, output: Path, parameters: dict[str, object]
) -> None:
    """Record outputs, input digests and the pinned data revision."""
    inputs = {
        relative: {
            "sha256": sha256_file(root / relative),
            "bytes": (root / relative).stat().st_size,
        }
        for relative in (
            LAB_XRF_TRIPS1_4,
            LAB_XRF_TRIP5,
            ALPHA_CACHE,
            GENUS_CACHE,
            *GEODATA,
        )
    }
    names = sorted(
        path.name
        for path in output.iterdir()
        if path.is_file() and path.name not in {"run_manifest.json", "SHA256SUMS"}
    )
    manifest = {
        "artifact_root": "analysis/v3/xrf_community_rescue",
        "files": [
            {
                "bytes": (output / name).stat().st_size,
                "path": name,
                "sha256": sha256_file(output / name),
            }
            for name in names
        ],
        "inputs": inputs,
        "data_repository": data_repository_lock(root),
        "parameters": parameters,
        "provenance": {
            "script": "analysis/v3/xrf_community_rescue.py",
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
                "scikit_learn": sklearn.__version__,
                "statsmodels": statsmodels.__version__,
            },
            "command": "python analysis/v3/xrf_community_rescue.py "
            + " ".join(sys.argv[1:]).replace(str(root) + "/", ""),
        },
        "schema_version": "1.1",
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    checksums = [*names, "run_manifest.json"]
    (output / "SHA256SUMS").write_text(
        "".join(
            f"{sha256_file(output / name)}  {name}\n" for name in sorted(checksums)
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--detection-threshold", type=float, default=0.75)
    args = parser.parse_args()
    analyse(
        args.project_root.resolve(),
        args.output_dir.resolve(),
        permutations=args.permutations,
        detection_threshold=args.detection_threshold,
    )


if __name__ == "__main__":
    main()
