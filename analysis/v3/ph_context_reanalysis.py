#!/usr/bin/env python3
"""Recompute site-level pH context from independently aggregated assay groups.

Consumes the corrected campaign/site/compartment pH join and the existing
200-genus site CLR response. Site pH is the equal-weight mean of assay-group
means. This is a descriptive site-level association, not a specimen match.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from scipy import stats
from claim_rescue import bh_fdr


def site_ph(groups):
    keys = ["trip", "site", "compartment"]
    if groups.duplicated(keys).any():
        raise ValueError("Repeated assay-group keys would change site weighting")
    if groups.ph.isna().any():
        raise ValueError("pH group means must be observed")
    return groups.groupby("site").ph.mean()


def run(groups_path, clr_path, reference_path, output):
    groups = pd.read_csv(groups_path, sep="\t")
    ph = site_ph(groups)
    clr = pd.read_csv(clr_path, sep="\t", index_col="site")
    if clr.index.duplicated().any():
        raise ValueError("Duplicate site CLR responses")
    sites = ph.index.intersection(clr.index).sort_values()
    reference = pd.read_csv(reference_path, sep="\t").set_index("genus")
    if set(clr.columns) != set(reference.index):
        raise ValueError("Genus selection differs from the route analysis")
    rows = []
    for genus in clr.columns:
        rho, p = stats.spearmanr(ph.loc[sites], clr.loc[sites, genus])
        rows.append(dict(genus=genus, phylum=reference.loc[genus, "phylum"],
                         n_sites=len(sites), spearman_rho_site_ph=float(rho), p_value=float(p),
                         route_supported=bool(reference.loc[genus, "supported_q_lt_0_05"]),
                         route_rho=reference.loc[genus, "spearman_rho_route_position"]))
    results = pd.DataFrame(rows)
    results["q_bh_200"] = bh_fdr(results.p_value)
    results["supported_q_lt_0_05"] = results.q_bh_200 < .05
    results = results.sort_values("spearman_rho_site_ph")
    coordinates = groups.groupby("site").transect_km.agg(["first", "nunique"])
    if (coordinates["nunique"] != 1).any():
        raise ValueError("Route position differs within a site")
    rho, p = stats.spearmanr(ph.loc[sites], coordinates.loc[sites, "first"])
    summary = dict(n_sites=len(sites), n_groups=len(groups),
                   n_assays=int(groups.n_ph_specimens.sum()),
                   site_pH_weighting="equal weight per campaign/site/compartment assay mean",
                   pH_route_rho=float(rho), pH_route_p=float(p),
                   supported_genera=int(results.supported_q_lt_0_05.sum()),
                   supported_also_route=int((results.supported_q_lt_0_05 & results.route_supported).sum()),
                   null_assumption="ordinary site-level rank correlation; spatial dependence unadjusted",
                   inputs={str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in (groups_path, clr_path, reference_path)},
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output.mkdir(parents=True, exist_ok=True)
    results.to_csv(output / "ph_genus_correlations.tsv", sep="\t", index=False)
    ph.rename("mean_ph").to_csv(output / "site_ph.tsv", sep="\t")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("groups", "clr", "reference", "output"):
        parser.add_argument("--" + argument, type=Path, required=True)
    args = parser.parse_args()
    run(args.groups, args.clr, args.reference, args.output)
