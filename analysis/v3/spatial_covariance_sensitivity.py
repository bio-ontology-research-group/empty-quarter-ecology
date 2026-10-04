#!/usr/bin/env python3
"""Spatial-covariance rotation sensitivity on the frozen 60-site composition.

The covariance family is set from coordinates and group membership before the
observed response is tested. It propagates campaign/compartment centring and
unequal site averaging. Fixed-covariance rotation inference assumes a Gaussian
matrix-normal null with separable row/column covariance; the column covariance
may be singular and is conditioned out. The finite-model p envelope is scoped
to this explicit family. No covariance parameter is fitted to obtain a p-value.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
from scipy import linalg, stats
from scipy.spatial.distance import pdist, squareform

import spatial_turnover_rescue as original

ROOT = Path(__file__).resolve().parents[2]
SEED = 202609091
KERNELS = ("exponential", "matern_3_2", "matern_5_2")
FRACTIONS = (0.01, 0.1, 0.5, 0.9)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def averaging_operators(metadata: pd.DataFrame, sites: np.ndarray):
    """Exact linear operator for the original centring/averaging algorithm."""
    n, g = len(sites), len(metadata)
    incidence = (metadata["site"].to_numpy()[:, None] == sites[None, :]).astype(float)
    average = incidence.T / incidence.sum(axis=0)[:, None]
    center = np.eye(g)
    for indices in metadata.groupby(["campaign", "compartment"], sort=True).indices.values():
        center[np.ix_(indices, indices)] -= 1 / len(indices)
    transform = average @ center
    return transform, transform @ incidence, transform @ transform.T


def kernel(distance: np.ndarray, kind: str, length: float) -> np.ndarray:
    scaled = distance / length
    if kind == "exponential":
        return np.exp(-scaled)
    if kind == "matern_3_2":
        x = np.sqrt(3) * scaled
        return (1 + x) * np.exp(-x)
    if kind == "matern_5_2":
        x = np.sqrt(5) * scaled
        return (1 + x + x * x / 3) * np.exp(-x)
    raise ValueError(kind)


def unit_trace(covariance: np.ndarray) -> np.ndarray:
    return covariance / (np.trace(covariance) / len(covariance))


def covariance_square_root(covariance: np.ndarray) -> np.ndarray:
    """Simulate on the same known covariance support used by the test."""
    values, vectors = np.linalg.eigh(covariance)
    tolerance = len(covariance) * np.finfo(float).eps * values[-1]
    if values[0] < -tolerance:
        raise ValueError("Simulation covariance is not positive semidefinite")
    keep = values > tolerance
    return (vectors[:, keep] * np.sqrt(values[keep])) @ vectors[:, keep].T


def propagated_covariance(distance, shared_operator, noise, kind, length, fraction):
    spatial = shared_operator @ kernel(distance, kind, length) @ shared_operator.T
    result = (1 - fraction) * unit_trace(spatial) + fraction * unit_trace(noise)
    return (result + result.T) / 2


def null_projection(covariance: np.ndarray, design: np.ndarray):
    """Whiten known covariance support, then remove any remaining intercept.

    Propagated centring imposes one deterministic weighted-sum constraint.
    Its zero eigenvalue is a known design constraint, not an estimated/ridged
    covariance component. The propagated intercept is identically zero.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = len(covariance) * np.finfo(float).eps * eigenvalues[-1]
    if eigenvalues[0] < -tolerance:
        raise ValueError("Covariance is not positive semidefinite")
    keep = eigenvalues > tolerance
    whitening = (eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])).T
    nuisance = whitening @ design[:, :1]
    basis = linalg.null_space(nuisance.T) if np.linalg.norm(nuisance) > 1e-10 else np.eye(keep.sum())
    transform = basis.T @ whitening
    tested = transform @ design[:, 1:]
    q, _ = np.linalg.qr(tested, mode="reduced")
    rank = np.linalg.matrix_rank(tested)
    if rank != design.shape[1] - 1:
        raise ValueError("Route design rank changed after whitening")
    return transform, q[:, :rank]


def random_subspaces(rng, draws: int, dimension: int, rank: int):
    values = rng.normal(size=(draws, dimension, rank))
    # Gaussian QR produces a uniform unoriented subspace; column signs cancel
    # in the squared projection statistic used here.
    return np.linalg.qr(values, mode="reduced")[0]


def rotation_test(response, transform, tested, draws, rng):
    whitened = transform @ response
    tss = np.square(whitened).sum()
    observed = np.square(tested.T @ whitened).sum() / tss
    spaces = random_subspaces(rng, draws, len(whitened), tested.shape[1])
    projection = np.einsum("bnq,np->bqp", spaces, whitened, optimize=True)
    null = np.square(projection).sum(axis=(1, 2)) / tss
    exceed = int(np.count_nonzero(null >= observed - 1e-14))
    p = (exceed + 1) / (draws + 1)
    interval = stats.binomtest(exceed, draws).proportion_ci(method="exact")
    return {"whitened_route_r2": float(observed), "rotation_p": p,
            "null_exceedances": exceed, "rotations": draws,
            "monte_carlo_tail_ci_low": float(interval.low),
            "monte_carlo_tail_ci_high": float(interval.high)}


def plan_definition(coords, metadata, observed_draws, calibration_replicates, calibration_draws):
    distance = squareform(pdist(coords[["x_km", "y_km"]]))
    nearest = distance + np.diag(np.full(len(distance), np.inf))
    diameter = float(distance.max())
    lengths = [float(np.median(nearest.min(axis=1))),
               *[diameter * f for f in (0.05, 0.1, 0.2, 0.5, 1.0)]]
    models = [{"model": "group_noise_only", "kernel": "group_noise", "range_km": None,
               "group_noise_fraction": 1.0}]
    for kind in KERNELS:
        for index, length in enumerate(lengths):
            for fraction in FRACTIONS:
                models.append({"model": f"{kind}_range{index}_noise{fraction:g}",
                               "kernel": kind, "range_km": length,
                               "group_noise_fraction": fraction})
    return {
        "schema_version": "spatial-covariance-sensitivity-v1.0.0",
        "seed": SEED, "site_count": len(coords), "group_count": len(metadata),
        "range_rule": "median nearest-site distance; diameter times .05,.1,.2,.5,1",
        "noise_fraction_definition": "fraction of mean marginal site variance contributed by propagated independent group noise",
        "range_km": lengths, "diameter_km": diameter,
        "kernels": {"exponential": "exp(-d/l)", "matern_3_2": "(1+sqrt(3)d/l)exp(-sqrt(3)d/l)",
                    "matern_5_2": "(1+sqrt(5)d/l+5d²/(3l²))exp(-sqrt(5)d/l)"},
        "models": models, "reference": "iid_site_identity_covariance",
        "observed_rotations": observed_draws,
        "calibration_replicates": calibration_replicates, "calibration_rotations": calibration_draws,
        "calibration_cases": ["iid_sites", "group_noise_only", "matched_exp_medium",
                              "matched_matern_smooth", "off_grid_range_noise",
                              "nonstationary_site_field", "nonseparable_two_site_fields",
                              "heavy_tailed_group_innovations"],
        "calibration_column_covariance": "three response factors with eigenvalues 1,.3,.1, mixed by a fixed orthogonal matrix",
        "nuisance": "whitened intercept; original route design has linear and quadratic terms",
        "model_scope": "Gaussian group/site fields sharing one cross-genus covariance; original selected genus set and CLR transformation held fixed",
        "coverage": "A averages groups within site, B centres campaign-by-compartment groups; H=ABZ and V=ABB' A' propagate unequal group membership",
        "scope_of_envelope": "maximum p across the finite propagated-covariance family; a model sensitivity envelope, not a universal spatial p-value",
        "calibration_uncertainty": "exact binomial intervals for rejection frequencies; independent null datasets and rotation draws for each replicate",
        "references": ["https://doi.org/10.1007/s11222-005-4789-5",
                       "https://doi.org/10.1093/bioinformatics/btab063"],
    }


def covariance_family(plan, distance, shared, noise):
    covariances = []
    for model in plan["models"]:
        if model["kernel"] == "group_noise":
            covariance = unit_trace(noise)
        else:
            covariance = propagated_covariance(distance, shared, noise, model["kernel"],
                                                model["range_km"], model["group_noise_fraction"])
        # Reject ill-conditioned/indefinite models explicitly; never repair the
        # covariance silently or drop a model after seeing its probability.
        eigenvalues = np.linalg.eigvalsh(covariance)
        tolerance = len(covariance) * np.finfo(float).eps * eigenvalues[-1]
        if eigenvalues[0] < -tolerance or np.count_nonzero(eigenvalues > tolerance) != len(covariance) - 1:
            raise ValueError(f"Covariance does not have the known rank-59 support: {model['model']}")
        covariances.append(covariance)
    return covariances


def calibration_responses(case, plan, distance, transform, shared, noise, count, rng):
    """Generate Gaussian and deliberate misspecification nulls with no route mean."""
    n = len(distance)
    diameter = plan["diameter_km"]
    def prop(kind, length, fraction):
        return propagated_covariance(distance, shared, noise, kind, length, fraction)
    covariance = None
    if case == "iid_sites":
        covariance = np.eye(n)
    elif case == "group_noise_only":
        covariance = unit_trace(noise)
    elif case == "matched_exp_medium":
        covariance = prop("exponential", diameter * .2, .5)
    elif case == "matched_matern_smooth":
        covariance = prop("matern_3_2", diameter * .5, .1)
    elif case == "off_grid_range_noise":
        covariance = prop("exponential", diameter / 3, .35)
    elif case == "nonstationary_site_field":
        weight = np.linspace(0, 1, n)
        spatial = (np.outer(weight, weight) * kernel(distance, "matern_3_2", diameter * .5)
                   + np.outer(1 - weight, 1 - weight) * kernel(distance, "exponential", diameter * .05))
        diagonal = np.sqrt(np.diag(spatial))
        spatial /= np.outer(diagonal, diagonal)
        covariance = .5 * unit_trace(shared @ spatial @ shared.T) + .5 * unit_trace(noise)
    if covariance is not None:
        values = np.einsum("ij,rjp->rip", covariance_square_root(covariance),
                           rng.normal(size=(count, n, 3)))
    elif case == "nonseparable_two_site_fields":
        values = np.empty((count, n, 3))
        for axis, (kind, length, fraction) in enumerate((
            ("exponential", diameter * .05, .1),
            ("matern_3_2", diameter * .5, .5),
            ("matern_5_2", diameter * .1, .9))):
            values[:, :, axis] = rng.normal(size=(count, n)) @ covariance_square_root(prop(kind, length, fraction)).T
    elif case == "heavy_tailed_group_innovations":
        values = np.einsum("ng,rgp->rnp", transform,
                           rng.standard_t(5, size=(count, transform.shape[1], 3)) / np.sqrt(5 / 3))
    else:
        raise ValueError(case)
    # An arbitrary non-diagonal cross-response covariance. Rotation preserves
    # it jointly; no responses are permuted or simulated independently in a test.
    mixing, _ = np.linalg.qr(np.array([[1., 2., 1.], [2., -1., 3.], [1., 1., -2.]]))
    values = (values * np.sqrt([1., .3, .1])) @ mixing
    return values, covariance


def calibrate(plan, distance, averaging, shared, noise, design, family):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 2]))
    count, draws = plan["calibration_replicates"], plan["calibration_rotations"]
    n = len(distance)
    propagated_design = shared @ design
    transforms = [null_projection(np.eye(n), design)] + [null_projection(c, propagated_design) for c in family]
    rows, replicate_rows = [], []
    for case in plan["calibration_cases"]:
        responses, oracle = calibration_responses(case, plan, distance, averaging, shared, noise, count, rng)
        probabilities = np.empty((count, len(transforms)))
        # Reuse subspaces across covariance models, but give each simulated
        # dataset its own independent subspaces. This keeps model comparisons
        # paired without sharing Monte Carlo draws across calibration replicates.
        spaces = np.stack([random_subspaces(rng, draws, n - 1, design.shape[1] - 1)
                           for _ in range(count)])
        for model_index, (whitening, tested) in enumerate(transforms):
            y = np.einsum("ij,rjp->rip", whitening, responses, optimize=True)
            observed = np.square(np.einsum("nq,rnp->rqp", tested, y)).sum(axis=(1, 2))
            projection = np.einsum("rbnq,rnp->rbqp", spaces, y, optimize=True)
            null = np.square(projection).sum(axis=(2, 3))
            probabilities[:, model_index] = (1 + (null >= observed[:, None] - 1e-14).sum(axis=1)) / (draws + 1)
        comparisons = {"iid_site_rotation": probabilities[:, 0],
                       "propagated_group_noise_rotation": probabilities[:, 1],
                       "finite_family_maximum": probabilities[:, 1:].max(axis=1)}
        if oracle is not None:
            whitening, tested = null_projection(oracle, design if case == "iid_sites" else propagated_design)
            y = np.einsum("ij,rjp->rip", whitening, responses, optimize=True)
            observed = np.square(np.einsum("nq,rnp->rqp", tested, y)).sum(axis=(1, 2))
            null = np.square(np.einsum("rbnq,rnp->rbqp", spaces, y, optimize=True)).sum(axis=(2, 3))
            comparisons["generating_covariance_rotation"] = (1 + (null >= observed[:, None] - 1e-14).sum(axis=1)) / (draws + 1)
        # Reproduce the original unrestricted site-label null as a comparator.
        centered = responses - responses.mean(axis=1, keepdims=True)
        original_q = np.linalg.qr(design, mode="reduced")[0][:, 1:]
        observed = np.square(np.einsum("nq,rnp->rqp", original_q, centered)).sum(axis=(1, 2))
        exceedances = np.zeros(count, dtype=int)
        for _ in range(draws):
            permutations = np.stack([rng.permutation(n) for _ in range(count)])
            permuted = np.take_along_axis(centered, permutations[:, :, None], axis=1)
            statistic = np.square(np.einsum("nq,rnp->rqp", original_q, permuted)).sum(axis=(1, 2))
            exceedances += statistic >= observed - 1e-14
        comparisons["unrestricted_site_permutation"] = (1 + exceedances) / (draws + 1)
        for method, pvalues in comparisons.items():
            rejects = int((pvalues <= .05).sum())
            ci = stats.binomtest(rejects, count).proportion_ci(method="exact")
            rows.append({"null_case": case, "method": method, "replicates": count,
                         "rotations_per_replicate": draws, "rejections_at_0_05": rejects,
                         "rejection_fraction": rejects / count,
                         "binomial_ci_low": ci.low, "binomial_ci_high": ci.high})
            replicate_rows.extend({"null_case": case, "method": method, "replicate": i, "p": float(p)}
                                  for i, p in enumerate(pvalues))
        print(f"calibrated {case}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(replicate_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--rotations", type=int, default=9999)
    parser.add_argument("--calibration-replicates", type=int, default=500)
    parser.add_argument("--calibration-rotations", type=int, default=199)
    args = parser.parse_args()
    root, out = args.project_root.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sources = {"genus_counts": root / "analysis/v2/review/cache/genus_counts.tsv",
               "coordinates": root / "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv",
               "original_analysis": Path(original.__file__).resolve()}
    counts, metadata, info = original.load_grouped_counts(sources["genus_counts"], 2000)
    coords = pd.read_csv(sources["coordinates"], sep="\t").sort_values("site")
    sites = coords["site"].to_numpy()
    plan = plan_definition(coords, metadata, args.rotations, args.calibration_replicates, args.calibration_rotations)
    plan["inputs"] = {name: {"path": str(path.relative_to(root)), "sha256": digest(path)}
                      for name, path in sources.items()}
    plan_path = out / "frozen_design.json"
    if args.plan_only:
        if plan_path.exists():
            raise FileExistsError("Design already frozen; do not overwrite")
        plan["frozen_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(plan_path, plan)
        print(f"Froze geometry/coverage design before response tests: {plan_path}")
        return
    frozen = json.loads(plan_path.read_text())
    if {k: v for k, v in frozen.items() if k != "frozen_at_utc"} != plan:
        raise ValueError("Inputs or parameters differ from frozen design")
    addendum_path = out / "frozen_design_addendum.json"
    addendum = json.loads(addendum_path.read_text())
    if addendum["known_covariance_rank"] != len(sites) - 1:
        raise ValueError("Known support addendum inconsistent with site design")
    taxa = original.rank_taxa(counts, .20)[:200]
    response_sites, response, groups = original.site_level_clr(counts, metadata, taxa, None, .5)
    if not np.array_equal(response_sites, sites):
        raise ValueError("Site order mismatch")
    design = original.design_matrix(coords["transect_km"].to_numpy(), 2)
    if np.linalg.matrix_rank(design) != 3:
        raise ValueError("Original quadratic design rank is not three")
    averaging, shared, noise = averaging_operators(metadata, sites)
    logged = np.log(counts.loc[taxa].T.to_numpy() + .5)
    clr = logged - logged.mean(axis=1, keepdims=True)
    if not np.allclose(averaging @ clr, response, atol=1e-12):
        raise AssertionError("Linear operator does not reproduce original centring/averaging")
    distance = squareform(pdist(coords[["x_km", "y_km"]]))
    covariances = covariance_family(plan, distance, shared, noise)
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 1]))
    original_r2 = original.fit_multivariate(response, design)[0]
    output_rows = []
    definitions = [{"model": "iid_site_reference", "kernel": "iid_site", "range_km": None,
                    "group_noise_fraction": None}] + plan["models"]
    for model, covariance in zip(definitions, [np.eye(len(sites))] + covariances):
        is_reference = model["model"] == "iid_site_reference"
        transform, tested = null_projection(covariance, design if is_reference else shared @ design)
        result = rotation_test(response, transform, tested, args.rotations, rng)
        eigenvalues = np.linalg.eigvalsh(covariance)
        positive = eigenvalues[eigenvalues > len(sites) * np.finfo(float).eps * eigenvalues[-1]]
        output_rows.append({**model, "site_count": len(sites), "group_count": groups,
                            "tested_rank": tested.shape[1], "nuisance_rank_after_preprocessing": int(is_reference),
                            "covariance_support_rank": len(positive),
                            "rotation_dimension": transform.shape[0],
                            "covariance_support_condition_number": float(positive[-1] / positive[0]), **result})
    observed = pd.DataFrame(output_rows)
    observed.to_csv(out / "spatial_covariance_tests.tsv", sep="\t", index=False)
    pd.DataFrame(response, index=sites, columns=taxa).rename_axis("site").to_csv(out / "site_clr_response.tsv", sep="\t")
    metadata.to_csv(out / "group_cohort.tsv", sep="\t", index=False)
    counts_per_site = metadata.groupby("site").size().rename("n_groups").reset_index()
    counts_per_site.to_csv(out / "site_group_counts.tsv", sep="\t", index=False)
    print(observed[["model", "rotation_p"]].to_string(index=False), flush=True)
    calibration, calibration_replicates = calibrate(plan, distance, averaging, shared, noise, design, covariances)
    calibration.to_csv(out / "null_calibration.tsv", sep="\t", index=False)
    calibration_replicates.to_csv(out / "null_calibration_replicates.tsv", sep="\t", index=False)
    family_rows = observed.iloc[1:]
    summary = {"original_descriptive_r2": original_r2,
               "model_count": len(family_rows), "family_p_min": float(family_rows.rotation_p.min()),
               "family_p_max": float(family_rows.rotation_p.max()),
               "reference_iid_rotation_p": float(observed.iloc[0].rotation_p),
               "groups_per_site_min": int(counts_per_site.n_groups.min()),
               "groups_per_site_max": int(counts_per_site.n_groups.max()),
               "scope": plan["scope_of_envelope"], "input_info": info}
    write_json(out / "summary.json", summary)
    outputs = sorted(p for p in out.iterdir() if p.suffix in (".tsv", ".json") and p.name != "run_manifest.json")
    write_json(out / "run_manifest.json", {"finished_at_utc": datetime.now(timezone.utc).isoformat(),
               "generator": {"path": str(Path(__file__).resolve().relative_to(root)), "sha256": digest(Path(__file__))},
               "inputs": plan["inputs"], "frozen_design_sha256": digest(plan_path),
               "frozen_design_addendum_sha256": digest(addendum_path),
               "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
               "outputs": {p.name: digest(p) for p in outputs}})
    (out / "SHA256SUMS").write_text("".join(f"{digest(p)}  {p.name}\n" for p in sorted(out.iterdir()) if p.is_file() and p.name != "SHA256SUMS"))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
