from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3.measured_function_summary import (
    collapse_picrust_columns,
    denominator_audit,
    pathway_screen_rows,
    profile_correlations,
    write_outputs,
)


ROOT = Path(__file__).resolve().parents[1]


def test_picrust_duplicate_columns_are_summed_then_normalized():
    source = pd.DataFrame(
        {
            "c1a": [1.0, 3.0],
            "c1b": [2.0, 0.0],
            "c2": [2.0, 2.0],
        },
        index=["K01601", "K00001"],
    )
    metadata = pd.DataFrame(
        {
            "picrust2_col": ["c1a", "c1b", "c2"],
            "sample_id": ["s1", "s1", "s2"],
        }
    )
    collapsed, audit = collapse_picrust_columns(source, metadata)
    assert list(collapsed.columns) == ["s1", "s2"]
    assert np.allclose(collapsed["s1"], [0.5, 0.5])
    assert np.allclose(collapsed["s2"], [0.5, 0.5])
    s1 = next(row for row in audit if row["sample_id"] == "s1")
    assert s1["n_source_picrust_columns"] == 2
    assert s1["aggregation_rule"] == "sum_then_within_sample_normalize"


def test_profile_and_marker_correlations_are_separate():
    measured = pd.DataFrame(
        {
            "s1": [0.1, 0.6, 0.3],
            "s2": [0.2, 0.2, 0.6],
            "s3": [0.6, 0.3, 0.1],
        },
        index=["K01601", "K00001", "K00002"],
    )
    predicted = pd.DataFrame(
        {
            "s1": [0.2, 0.5, 0.3],
            "s2": [0.3, 0.1, 0.6],
            "s3": [0.4, 0.5, 0.1],
        },
        index=measured.index,
    )
    sample_rows, metrics, diagnostics = profile_correlations(
        measured,
        predicted,
        {"s1": "Surface", "s2": "Deep", "s3": "Rhizosphere"},
    )
    assert len(sample_rows) == 3
    marker = next(
        row
        for row in metrics
        if row["metric"] == "marker_across_sample_spearman"
    )
    assert marker["feature"] == "K01601"
    assert marker["n_samples"] == 3
    assert diagnostics["rubisco_measured_nonzero_samples"] == 3


def test_pathway_count_and_denominator_are_bounded():
    pathway = pd.DataFrame(
        {
            "genome": ["g1", "g2", "g3"],
            "CBB": [1, 1, 0],
            "rTCA": [0, 0, 0],
            "WL": [0, 0, 0],
            "3HP": [0, 0, 0],
            "HB": [0, 0, 1],
        }
    )
    rows = pathway_screen_rows(pathway)
    assert next(row for row in rows if row["pathway"] == "CBB")[
        "positive_genome_records"
    ] == 2
    filtered = pd.DataFrame({"genome": ["g1sta", "g2sta"]})
    audit = denominator_audit(pathway, filtered)
    assert audit["pathway_positive_labels_matched_to_filtered_table"] == 2
    assert audit["pathway_positive_labels_not_matched_to_filtered_table"] == 1
    assert audit["total_screened_genomes"] is None
    assert audit["genome_fraction_reportable"] is False


def _write_toy_inputs(base: Path) -> tuple[Path, Path, Path]:
    measured_dir = base / "measured"
    measured_dir.mkdir(parents=True)
    measured = pd.DataFrame(
        {
            "s1": [1.0, 4.0, 2.0],
            "s2": [2.0, 1.0, 5.0],
            "s3": [5.0, 2.0, 1.0],
        },
        index=["K01601", "K00001", "K00002"],
    )
    measured.index.name = "KO"
    measured.to_csv(
        measured_dir / "measured_ko_by_sample.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "compartment": ["Surface", "Deep", "Rhizosphere"],
        }
    ).to_csv(
        measured_dir / "measured_marker_by_sample.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        {
            "genome": ["g1", "g2"],
            "CBB": [1, 0],
            "rTCA": [0, 0],
            "WL": [0, 0],
            "3HP": [0, 0],
            "HB": [0, 1],
            "lineage": ["p__A", "p__B"],
        }
    ).to_csv(
        measured_dir / "genome_cfix_taxonomy.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame({"genome": ["g1sta"]}).to_csv(
        measured_dir / "filtered_genomes.tsv",
        sep="\t",
        index=False,
    )
    picrust = pd.DataFrame(
        {
            "c1a": [1.0, 4.0, 2.0],
            "c1b": [0.5, 2.0, 1.0],
            "c2": [2.0, 1.0, 5.0],
            "c3": [4.0, 2.0, 1.0],
        },
        index=measured.index,
    )
    picrust.index.name = "function"
    picrust_path = base / "picrust.tsv"
    picrust.to_csv(picrust_path, sep="\t")
    metadata_path = base / "metadata.tsv"
    pd.DataFrame(
        {
            "picrust2_col": ["c1a", "c1b", "c2", "c3"],
            "sample_id": ["s1", "s1", "s2", "s3"],
        }
    ).to_csv(metadata_path, sep="\t", index=False)
    return measured_dir, picrust_path, metadata_path


def test_output_bundle_is_byte_deterministic(tmp_path):
    measured_dir, picrust_path, metadata_path = _write_toy_inputs(
        tmp_path / "inputs"
    )
    first = tmp_path / "first"
    second = tmp_path / "second"
    for output in (first, second):
        summary = write_outputs(
            ROOT,
            measured_dir,
            picrust_path,
            metadata_path,
            output,
        )
        assert summary["cbb_joint_marker_screen"][
            "positive_genome_records"
        ] == 1
        assert summary["genome_denominator_audit"][
            "genome_fraction_reportable"
        ] is False
    first_files = sorted(path.name for path in first.iterdir())
    second_files = sorted(path.name for path in second.iterdir())
    assert first_files == second_files
    for name in first_files:
        assert (first / name).read_bytes() == (second / name).read_bytes()
