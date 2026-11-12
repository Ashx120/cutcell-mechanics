import numpy as np
import pytest

from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.solver import Material
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


def test_elastic_plane_stress_matches_existing_constitutive_matrix():
    material = J2Material()
    strain = np.array([1e-4, -2e-5, 3e-5])
    response = update_plane_stress(strain, MaterialState(), material, .01)
    elastic_matrix = Material(material.young, material.poisson).D
    np.testing.assert_allclose(response.stress, elastic_matrix @ strain, atol=1e-10)
    np.testing.assert_allclose(response.tangent, elastic_matrix, atol=1e-10)
    np.testing.assert_allclose(response.strain_zz, -material.poisson / (1-material.poisson) * strain[:2].sum())
    assert abs(response.stress_zz) < 1e-8


@pytest.mark.parametrize("exponent", [.5, 1., 2.])
def test_plane_stress_tangent_with_history(exponent):
    material = J2Material(relaxation_time=100, rate_exponent=exponent)
    previous = update_plane_stress([.003, -.0005, .002], MaterialState(), material, .1).state
    original_history = previous.plastic_strain.copy()
    strain = np.array([.004, -.001, .003])
    response = update_plane_stress(strain, previous, material, .02)
    numerical = np.empty((3, 3))
    for column in range(3):
        perturbation = np.eye(3)[column] * 1e-7
        plus = update_plane_stress(strain + perturbation, previous, material, .02)
        minus = update_plane_stress(strain - perturbation, previous, material, .02)
        numerical[:, column] = (plus.stress - minus.stress) / 2e-7
    assert abs(response.stress_zz) <= 1e-8
    assert response.state.accumulated_strain > previous.accumulated_strain
    assert np.linalg.norm(response.tangent-numerical) / np.linalg.norm(response.tangent) < 1e-6
    np.testing.assert_array_equal(previous.plastic_strain, original_history)
