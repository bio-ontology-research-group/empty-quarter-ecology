#!/usr/bin/env python3
"""Regenerate the supplementary PMA paired-richness figure from endpoint results."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source = root / "analysis/v3/pma_endpoint_results/pma_pair_endpoints.tsv"
    data = pd.read_csv(source, sep="\t")
    assert len(data) == 9 and set(data.rarefaction_depth) == {123897}
    columns = ["untreated_expected_rarefied_richness", "treated_expected_rarefied_richness"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(5.2, 3.6), layout="constrained")
    for group, colour in zip(("C1R", "C2R", "C2S"), ("#0072B2", "#009E73", "#D55E00")):
        values = data.loc[data.pair_id.str.startswith(group), columns].to_numpy()
        assert len(values) == 3
        for row in values:
            ax.plot([0, 1], row, "o-", color=colour, alpha=.4, lw=1.2, ms=5)
        ax.plot([0, 1], values.mean(axis=0), "D-", color=colour, lw=2.5, ms=7, label=group+" mean")
    ax.set_xticks([0,1], ["Untreated", "PMA treated"])
    ax.set(ylabel="Expected ASV richness", xlim=(-.08,1.08))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "fig_pma_richness.pdf"
    fig.savefig(output, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {"input": str(source.relative_to(root)), "input_sha256": digest(source),
                "generator": str(Path(__file__).resolve().relative_to(root)),
                "generator_sha256": digest(Path(__file__)), "output_sha256": digest(output),
                "definition": "Nine aliquot pairs from three parent soil tubes at two campsites; expected ASV richness at 123897 reads; group means are descriptive."}
    (args.output_dir / "pma_richness_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    main()
