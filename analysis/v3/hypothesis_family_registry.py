#!/usr/bin/env python3
"""Enumerate reported BH families and verify their adjusted probabilities.

The original result files retain the scientific identifiers. This registry
adds an explicit family identifier, membership, source hash and adjustment
method for the main paper and retained sensitivity tables. Single omnibus
tests and rainfall's maximum-statistic search are documented separately.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from claim_rescue import bh_fdr

SPECS = [
    ("alpha_original", "results/paired_compartment_effects.tsv", "p", "q_within_metric", ["metric"]),
    ("composition", "compartment_composition/compartment_location_results.tsv", "permutation_p", "q_primary_family", []),
    ("climate_alpha", "environment_associations/climate_alpha_correlations.tsv", "p_value", "q_global_9", []),
    ("climate_genus", "environment_associations/climate_genus_correlations.tsv", "p_value", "q_global_600", []),
    ("field_weather", "environment_associations/field_weather_adjusted_models.tsv", "p_value", "q_global_9", []),
    ("taxon_route", "taxon_context_corrected_20260909/transect_replacement.tsv", "p_value", "q_bh_200", []),
    ("compartment_genus", "biology_context_corrected_20260909/compartment_genus_family.tsv", "sign_flip_p", "q_bh_600", []),
    ("landform_genus", "biology_context_corrected_20260909/landform_genus_contrasts.tsv", "mannwhitney_p", "q_bh_200", ["adjustment"]),
    ("dune_route", "biology_context_corrected_20260909/route_genus_correlations_dune_sites.tsv", "p_value", "q_bh_200", []),
    ("dune_climate", "biology_context_corrected_20260909/climate_diversity_dune_sensitivity.tsv", "p_dune_44", "q_bh_9_dune_44", []),
    ("ph_genus", "biology_context_corrected_20260909/ph_genus_correlations.tsv", "p_value", "q_bh_200", []),
    ("xrf_genus", "biology_context_corrected_20260909/xrf_axis_genus_correlations.tsv", "p_value", "q_bh_200", []),
    ("pathway_profile", "picrust2_ecology/position_profile_tests.tsv", "permutation_p", "q_primary_three", []),
    ("pathway_compartment", "picrust2_ecology/pathway_position_effects.tsv", "p_value", "q_global_600", []),
    ("pathway_route", "picrust2_ecology/pathway_geographic_correlations.tsv", "p_value", "q_global_200", []),
    ("marker_compartment", "trait_genes/picrust_trait_compartment_contrasts.tsv", "sign_flip_p", "q_bh", []),
    ("marker_route", "trait_genes/picrust_trait_summary.tsv", "p_route", "q_bh_route", []),
    ("normalized_alpha", "evenness_decomposition/paired_contrasts.tsv", "p", "q_within_metric", ["metric"]),
    ("paired_alpha_sensitivity", "paired_alpha_sensitivity_20260909/paired_alpha_estimands.tsv", "wilcoxon_p", "q_three_contrasts", ["cohort", "metric", "estimand"]),
]


def run(root, output):
    families, members = [], []
    for prefix, relative, p_column, q_column, grouping in SPECS:
        path = root / relative
        frame = pd.read_csv(path, sep="\t")
        frame = frame[frame[q_column].notna()].copy()
        iterator = frame.groupby(grouping, dropna=False, sort=True) if grouping else [((), frame)]
        for key, part in iterator:
            key = key if isinstance(key, tuple) else (key,)
            family = prefix + (":" + ":".join(map(str, key)) if key else "")
            calculated = bh_fdr(part[p_column])
            if not np.allclose(calculated, part[q_column], atol=1e-8, rtol=1e-7):
                raise ValueError(f"Published BH values disagree with membership in {family}")
            families.append(dict(family=family, n_tests=len(part), method="Benjamini-Hochberg",
                                 source=relative, source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                 p_column=p_column, q_column=q_column,
                                 filter=dict(zip(grouping, key))))
            id_columns = [c for c in ("trip", "metric", "comparison", "contrast", "genus", "pathway",
                                      "trait", "climate_variable", "weather_variable", "response", "diversity_measure")
                          if c in part]
            for index, row in part.iterrows():
                identity = {column: str(row[column]) for column in id_columns}
                members.append(dict(family=family, source=relative, source_data_row=int(index)+1,
                                    hypothesis=json.dumps(identity, sort_keys=True),
                                    p=float(row[p_column]), q=float(row[q_column])))
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(members).to_csv(output / "hypothesis_members.tsv", sep="\t", index=False)
    registry = dict(schema_version="1.0", families=families,
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    other_tests={"omnibus_tests": "Reported individually with their permutation scheme",
                                 "rainfall_search": "Maximum statistic over 3 endpoints x 60 peak lags; shared circular shifts within expedition; 19,999 draws; seed 20260804",
                                 "distance_slope_contrasts": "Maximum statistic over two compartment slope contrasts; 9,999 whole-site-label permutations"})
    (output / "hypothesis_families.json").write_text(json.dumps(registry, indent=2) + "\n")
    print(f"Verified {len(families)} BH families and {len(members)} memberships")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.analysis_root, args.output)
