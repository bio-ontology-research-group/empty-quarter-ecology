#!/usr/bin/env python3
"""Aitchison ordination of the grouped genus profiles used by the spatial test.

The ordination is descriptive. It uses the same cohort and genus set as
``spatial_turnover_rescue.py``:

1. sequencing replicates are summed within campaign x site x compartment and
   groups below the read threshold are discarded;
2. the 200 most abundant eligible genera are transformed with a centred log
   ratio (CLR) after adding the pseudocount; and
3. the CLR profiles are centred by genus and decomposed by singular value
   decomposition, which is a principal component analysis in Aitchison
   geometry.

No campaign or compartment effect is removed, so the axes show all sources of
compositional variation. The ordination supports no inference by itself.

Run: python analysis/v3/aitchison_ordination.py --output-dir analysis/v3/aitchison_ordination
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

import spatial_turnover_rescue as spatial


SCHEMA_VERSION = "1.0"
TAXON_COUNT = 200
COMPONENTS = 2


def aitchison_pca(
    counts: pd.DataFrame,
    taxa: list[str],
    pseudocount: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return group scores, genus loadings and variance fractions."""
    matrix = counts.loc[taxa].T.to_numpy(dtype=float)
    logged = np.log(matrix + pseudocount)
    clr = logged - logged.mean(axis=1, keepdims=True)
    centred = clr - clr.mean(axis=0, keepdims=True)
    left, singular, right = np.linalg.svd(centred, full_matrices=False)
    variance = np.square(singular) / np.square(singular).sum()
    return left * singular, right, variance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--counts", type=Path, default=None)
    parser.add_argument("--site-coordinates", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-group-reads", type=int, default=2000)
    parser.add_argument("--prevalence", type=float, default=0.20)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()

    root = args.project_root.resolve()
    counts_path = (
        args.counts or root / "analysis/v2/review/cache/genus_counts.tsv"
    ).resolve()
    coordinates_path = (
        args.site_coordinates
        or root
        / "analysis/v3/spatial_turnover_rescue/results/site_coordinates.tsv"
    ).resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    counts, metadata, input_info = spatial.load_grouped_counts(
        counts_path, args.minimum_group_reads
    )
    ranked_taxa = spatial.rank_taxa(counts, args.prevalence)
    if len(ranked_taxa) < TAXON_COUNT:
        raise ValueError(
            f"At least {TAXON_COUNT} eligible genera are required; "
            f"found {len(ranked_taxa)}"
        )
    taxa = ranked_taxa[:TAXON_COUNT]
    scores, loadings, variance = aitchison_pca(counts, taxa, args.pseudocount)

    coordinates = pd.read_csv(coordinates_path, sep="\t")
    table = metadata.merge(
        coordinates[["site", "transect_km"]],
        on="site",
        how="left",
        validate="many_to_one",
    )
    if table["transect_km"].isna().any():
        raise ValueError("A retained group has no transect coordinate")

    # Singular vectors have arbitrary sign. PC1 is oriented to increase
    # eastwards; later axes make their largest-magnitude loading positive.
    transect = table["transect_km"].to_numpy(dtype=float)
    for component in range(COMPONENTS):
        if component == 0:
            flip = np.corrcoef(scores[:, 0], transect)[0, 1] < 0
        else:
            row = loadings[component]
            flip = row[np.argmax(np.abs(row))] < 0
        if flip:
            scores[:, component] *= -1
            loadings[component] *= -1
    for component in range(COMPONENTS):
        table[f"pc{component + 1}"] = scores[:, component]

    spatial.write_tsv(
        output / "ordination_scores.tsv",
        table.to_dict("records"),
        ["campaign", "site", "compartment", "transect_km"]
        + [f"pc{component + 1}" for component in range(COMPONENTS)],
    )
    spatial.write_tsv(
        output / "ordination_loadings.tsv",
        (
            {
                "genus": genus,
                **{
                    f"pc{component + 1}": float(loadings[component, index])
                    for component in range(COMPONENTS)
                },
            }
            for index, genus in enumerate(taxa)
        ),
        ["genus"] + [f"pc{component + 1}" for component in range(COMPONENTS)],
    )
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": "descriptive_ordination_only",
        "method": (
            "Principal component analysis of centred log-ratio genus "
            "profiles (Aitchison geometry); no campaign or compartment "
            "effect removed"
        ),
        "n_groups": int(len(table)),
        "n_sites": int(table["site"].nunique()),
        "n_genera": TAXON_COUNT,
        "pseudocount": args.pseudocount,
        "prevalence_threshold": args.prevalence,
        "variance_explained": {
            f"pc{component + 1}": float(variance[component])
            for component in range(COMPONENTS)
        },
        "pc1_transect_spearman_rho": float(
            pd.Series(scores[:, 0]).corr(
                pd.Series(transect), method="spearman"
            )
        ),
        "input": input_info,
        "input_sha256": {
            "genus_counts": spatial.sha256_file(counts_path),
            "site_coordinates": spatial.sha256_file(coordinates_path),
        },
        "provenance": {
            "script": "analysis/v3/aitchison_ordination.py",
            "script_sha256": spatial.sha256_file(Path(__file__).resolve()),
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
        },
    }
    (output / "ordination_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
