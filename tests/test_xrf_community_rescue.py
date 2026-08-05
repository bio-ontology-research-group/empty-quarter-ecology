import numpy as np

from analysis.v3.xrf_community_rescue import (
    nested_multivariate_test,
    parse_sample_id,
    pcoa_coordinates,
)


def test_parse_sample_id_handles_all_campaigns_and_compartments():
    assert parse_sample_id("10Dr2") == (1, 10, "Deep")
    assert parse_sample_id("T2Sr1") == (2, 2, "Surface")
    assert parse_sample_id("F3PRr2") == (3, 3, "Rhizosphere")
    assert parse_sample_id("V4Dr1", default_trip=5) == (5, 4, "Deep")


def test_pcoa_coordinates_recover_positive_euclidean_axes():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0], [1.0, 2.0]])
    distance = np.sqrt(
        np.square(points[:, None, :] - points[None, :, :]).sum(axis=2)
    )
    coordinates, retained = pcoa_coordinates(distance, 0.95)
    assert coordinates.shape == (4, 2)
    assert retained == 1.0


def test_nested_multivariate_test_is_deterministic():
    predictor = np.array([-1.0, 1.0, -0.8, 0.8, -1.2, 1.2])
    response = np.column_stack([predictor, predictor**3])
    reduced = np.column_stack(
        [np.ones(6), np.repeat([0.0, 1.0, 2.0], 2)]
    )
    sites = np.repeat([1, 2, 3], 2)
    first = nested_multivariate_test(
        response, reduced, predictor, sites, permutations=49, seed=7
    )
    second = nested_multivariate_test(
        response, reduced, predictor, sites, permutations=49, seed=7
    )
    assert first == second
    assert first["partial_r2"] > 0.9
