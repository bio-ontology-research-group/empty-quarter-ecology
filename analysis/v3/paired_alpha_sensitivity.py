#!/usr/bin/env python3
"""Compare marginal and campaign-matched site-level alpha contrasts.

The original estimates average each compartment over its available campaigns
before subtraction. This sensitivity forms differences within campaign and
site first, then averages those differences within site. Both endpoint-specific
and common-profile cohorts are reported; the original outputs are preserved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from claim_rescue import bootstrap_mean_ci, bh_fdr

METRICS = ("shannon", "richness_rare", "normalized_shannon")
CONTRASTS = (("Deep", "Surface"), ("Rhizosphere", "Surface"),
             ("Rhizosphere", "Deep"))


def site_differences(grouped, metric, first, second, matched):
    if matched:
        wide = grouped.pivot(index=["Trip", "Site"], columns="Type", values=metric)
        pair = wide.reindex(columns=[first, second]).dropna()
        differences = pair[first] - pair[second]
        return differences.groupby(level="Site").mean(), len(pair)
    wide = grouped.groupby(["Site", "Type"])[metric].mean().unstack()
    pair = wide.reindex(columns=[first, second]).dropna()
    return pair[first] - pair[second], None


def run(alpha_path: Path, output: Path, resamples: int, seed: int):
    alpha = pd.read_csv(alpha_path, sep="\t", index_col=0)
    alpha["Type"] = alpha["Type"].replace({"Rhizo": "Rhizosphere"})
    alpha = alpha[alpha.Site.between(1, 60) & alpha.Trip.between(1, 5)].copy()
    alpha["normalized_shannon"] = alpha.shannon / np.log(alpha.richness_rare)
    rows, differences = [], []
    for cohort in ("endpoint_specific", "common_profiles"):
        eligible = alpha if cohort == "endpoint_specific" else alpha.dropna(subset=list(METRICS))
        grouped = eligible.groupby(["Trip", "Site", "Type"], as_index=False)[list(METRICS)].mean()
        for metric in METRICS:
            for matched in (False, True):
                for first, second in CONTRASTS:
                    values, n_pairs = site_differences(grouped, metric, first, second, matched)
                    low, high = bootstrap_mean_ci(values.to_numpy(), resamples=resamples, seed=seed)
                    p = 1.0 if np.all(values == 0) else float(stats.wilcoxon(values).pvalue)
                    common = dict(cohort=cohort, metric=metric,
                                  estimand="campaign_matched" if matched else "marginal",
                                  contrast=f"{first}-{second}")
                    rows.append(dict(**common, n_sites=len(values), n_campaign_site_pairs=n_pairs,
                                     mean=float(values.mean()), ci_low=low, ci_high=high,
                                     wilcoxon_p=p))
                    differences.extend(dict(**common, site=int(site), difference=float(value))
                                       for site, value in values.items())
    result = pd.DataFrame(rows)
    # Each sensitivity/cohort/metric has the three compartment contrasts.
    result["q_three_contrasts"] = result.groupby(
        ["cohort", "metric", "estimand"], group_keys=False
    )["wilcoxon_p"].transform(lambda values: bh_fdr(values))
    output.mkdir(parents=True, exist_ok=True)
    result.to_csv(output / "paired_alpha_estimands.tsv", sep="\t", index=False)
    pd.DataFrame(differences).to_csv(output / "site_differences.tsv", sep="\t", index=False)
    finite = alpha.dropna(subset=["normalized_shannon", "depth"])
    depth_correlation = stats.spearmanr(finite.normalized_shannon, finite.depth)
    manifest = {
        "source": str(alpha_path),
        "source_sha256": hashlib.sha256(alpha_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "bootstrap_resamples": resamples, "seed": seed,
        "bootstrap_unit": "site", "interval": "percentile 95%",
        "wilcoxon_null": "symmetric distribution of site differences about zero",
        "adjustment_family": "three contrasts per metric, cohort and estimand",
        "normalized_shannon_range": [float(finite.normalized_shannon.min()),
                                     float(finite.normalized_shannon.max())],
        "normalized_shannon_depth_spearman": float(depth_correlation.statistic),
        "scope": "descriptive sensitivity; retains original marginal estimates",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(result.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=20260723)
    args = parser.parse_args()
    run(args.alpha, args.output, args.resamples, args.seed)
