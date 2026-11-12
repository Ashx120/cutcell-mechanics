import numpy as np
import pytest

from immersed_elasticity.viscoplasticity import J2Material, MaterialState, update_material


def shear_strain(gamma):
    strain = np.zeros((3, 3))
    strain[0, 1] = strain[1, 0] = gamma / 2
    return strain


def test_elastic_shear_and_hydrostatic_loading():
    material = J2Material()
    previous = MaterialState()
    stress, state = update_material(shear_strain(1e-4), previous, material, .01)
    np.testing.assert_allclose(stress[0, 1], material.shear_modulus * 1e-4)
    assert state.accumulated_strain == 0
    stress, state = update_material(.01 * np.eye(3), previous, material, .01)
    np.testing.assert_allclose(stress, .03 * material.bulk_modulus * np.eye(3))
    assert state.accumulated_strain == 0


def test_perfect_plastic_m1_matches_closed_form_backward_euler_update():
    material = J2Material(saturation_increment=0, rate_exponent=1)
    gamma, dt = .01, .1
    q_trial = np.sqrt(3) * material.shear_modulus * gamma
    exact_increment = (q_trial - material.initial_yield) / (
        3 * material.shear_modulus + material.relaxation_time * material.initial_yield / dt)
    previous = MaterialState()
    stress, state = update_material(shear_strain(gamma), previous, material, dt)
    np.testing.assert_allclose(state.accumulated_strain, exact_increment, rtol=1e-11)
    expected_tau = (q_trial - 3 * material.shear_modulus * exact_increment) / np.sqrt(3)
    np.testing.assert_allclose(stress[0, 1], expected_tau, rtol=1e-11)
    np.testing.assert_array_equal(previous.plastic_strain, np.zeros((3, 3)))
    assert previous.accumulated_strain == 0
    np.testing.assert_allclose(np.trace(state.plastic_strain), 0, atol=1e-14)


def test_hold_relaxes_stress_and_elastic_unloading_retains_history():
    material = J2Material()
    strain = shear_strain(.01)
    stress, state = update_material(strain, MaterialState(), material, .01)
    initial_tau, initial_p = stress[0, 1], state.accumulated_strain
    for _ in range(100):
        stress, state = update_material(strain, state, material, .01)
    assert 0 < stress[0, 1] < initial_tau
    assert state.accumulated_strain > initial_p
    relaxed_tau = stress[0, 1]
    unloaded, final_state = update_material(shear_strain(.0099), state, material, .01)
    np.testing.assert_allclose(unloaded[0, 1], relaxed_tau - material.shear_modulus * .0001)
    np.testing.assert_array_equal(final_state.plastic_strain, state.plastic_strain)
    assert final_state.accumulated_strain == state.accumulated_strain


@pytest.mark.parametrize("exponent", [.5, 1., 2.])
def test_analytical_tangent_with_nonzero_history_and_mixed_loading(exponent):
    from immersed_elasticity.viscoplasticity import (
        engineering_strain_basis, stress_vector, update_with_tangent,
    )

    material = J2Material(relaxation_time=100, rate_exponent=exponent)
    _, previous = update_material(shear_strain(.01), MaterialState(), material, .1)
    strain = np.array([[.004, .006, .001], [.006, -.001, .002], [.001, .002, 0.]])
    _, state, tangent = update_with_tangent(strain, previous, material, .01)
    assert state.accumulated_strain > previous.accumulated_strain
    numerical = np.empty((6, 6))
    for column in range(6):
        perturbation = 1e-7 * engineering_strain_basis(column)
        plus, _ = update_material(strain + perturbation, previous, material, .01)
        minus, _ = update_material(strain - perturbation, previous, material, .01)
        numerical[:, column] = (stress_vector(plus) - stress_vector(minus)) / 2e-7
    assert np.linalg.norm(tangent - numerical) / np.linalg.norm(tangent) < 1e-6
    np.testing.assert_allclose(tangent, tangent.T, rtol=1e-12, atol=1e-9)
