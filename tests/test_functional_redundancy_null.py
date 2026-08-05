from pathlib import Path

import numpy as np

from analysis.v3.functional_redundancy_null import (
    functional_profiles,
    normalise_genome,
    run_null,
    sample_compartment,
)


def test_name_and_compartment_parsing():
    assert normalise_genome("/tmp/a.fasta.gz") == "a"
    assert normalise_genome("b_assembly") == "b_assembly"
    assert sample_compartment("F12PRr3") == "Rhizosphere"
    assert sample_compartment("V2Dr1") == "Deep"


def test_functional_projection_is_row_normalised():
    abundance = np.array([[0.75, 0.25], [0.2, 0.8]], dtype=float)
    annotation = np.array([[2, 0, 1], [0, 3, 1]], dtype=float)
    projected = functional_profiles(abundance, annotation, False)
    assert np.allclose(projected.sum(axis=1), 1)


def test_null_is_deterministic():
    abundance = np.array(
        [
            [0.9, 0.1, 0.0],
            [0.8, 0.2, 0.0],
            [0.0, 0.2, 0.8],
            [0.0, 0.1, 0.9],
        ]
    )
    annotation = np.eye(3)
    compartments = ["Surface", "Surface", "Deep", "Deep"]
    first = run_null(
        abundance,
        annotation,
        compartments,
        permutations=19,
        genome_normalised=False,
        seed=3,
    )
    second = run_null(
        abundance,
        annotation,
        compartments,
        permutations=19,
        genome_normalised=False,
        seed=3,
    )
    assert [row.keys() for row in first] == [row.keys() for row in second]
    for first_row, second_row in zip(first, second):
        for key in first_row:
            if isinstance(first_row[key], float):
                assert np.allclose(
                    first_row[key],
                    second_row[key],
                    equal_nan=True,
                )
            else:
                assert first_row[key] == second_row[key]


def test_unestimable_groups_do_not_receive_false_significance():
    abundance = np.eye(3)
    rows = run_null(
        abundance,
        np.eye(3),
        ["Surface", "Deep", "Deep"],
        permutations=19,
        genome_normalised=False,
        seed=3,
    )
    surface = next(row for row in rows if row["group"] == "Surface")
    assert np.isnan(surface["lower_tail_p"])
