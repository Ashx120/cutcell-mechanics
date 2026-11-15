import numpy as np

from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.thermal_expansion import update_thermal_plane_stress
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


def test_free_expansion_is_stress_free_in_all_directions():
    response = update_thermal_plane_stress([.0006, .0006, 0], MaterialState(),
                                          J2Material(), .1, temperature=70.)
    np.testing.assert_allclose(response.stress, 0., atol=1e-10)
    np.testing.assert_allclose(response.strain_zz, .0006, atol=1e-12)
    assert response.state.accumulated_strain == 0


def test_restrained_expansion_matches_plane_stress_elastic_solution():
    material = J2Material()
    response = update_thermal_plane_stress([0, 0, 0], MaterialState(), material, .1, temperature=70.)
    np.testing.assert_allclose(response.stress, [-180., -180., 0.], atol=1e-9)
    np.testing.assert_allclose(response.strain_zz, .0006 * (1+material.poisson)/(1-material.poisson))
    assert abs(response.stress_zz) < 1e-8
    assert response.state.accumulated_strain == 0


def test_zero_temperature_change_recovers_original_plane_stress_update():
    material = J2Material(relaxation_time=100)
    previous = MaterialState()
    strain = [.004, -.001, .002]
    expected = update_plane_stress(strain, previous, material, .1)
    actual = update_thermal_plane_stress(strain, previous, material, .1, temperature=20.)
    np.testing.assert_allclose(actual.stress, expected.stress)
    np.testing.assert_allclose(actual.tangent, expected.tangent)
    np.testing.assert_allclose(actual.state.plastic_strain, expected.state.plastic_strain)


def test_thermal_tangent_is_total_strain_derivative_at_fixed_temperature():
    material, previous = J2Material(relaxation_time=100), MaterialState()
    strain = np.array([.004, -.001, .002])
    response = update_thermal_plane_stress(strain, previous, material, .1, temperature=70.)
    numerical = np.empty((3, 3))
    for column in range(3):
        perturbation = np.eye(3)[column] * 1e-7
        plus = update_thermal_plane_stress(strain + perturbation, previous, material, .1, temperature=70.)
        minus = update_thermal_plane_stress(strain - perturbation, previous, material, .1, temperature=70.)
        numerical[:, column] = (plus.stress - minus.stress) / 2e-7
    assert np.linalg.norm(numerical-response.tangent) / np.linalg.norm(response.tangent) < 1e-6
