import numpy as np
import pytest

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.solver import structured_mesh, symmetry_constraints
from immersed_elasticity.viscoplasticity import J2Material
from scripts.thermal_viscoplastic_rectangle import run_check
from immersed_elasticity.solver import Material
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermoelastic import solve_thermoelastic
from immersed_elasticity.temperature_properties import with_temperature_dependent_elasticity, with_temperature_dependent_hardening
from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.viscoplasticity import MaterialState


@pytest.mark.parametrize("mode,property_temperature", [("paper-reference", 20.), ("paper-temperature", 400.)])
def test_thermal_hardening_matches_plastic_material_point_and_retains_history(mode, property_temperature):
    mesh = structured_mesh(1, 1)
    base = J2Material(relaxation_time=100.)
    model = NonlinearModel(mesh, base, elasticity="paper-temperature", hardening=mode)
    material = with_temperature_dependent_hardening(
        with_temperature_dependent_elasticity(base, 400.), property_temperature)
    previous = MaterialState()
    histories = model.initial_state().history
    for strain in (.003, .004):
        displacement = (strain * mesh.nodes).ravel()
        _, _, histories, stress = model.evaluate(displacement, histories, .1,
                                               nodal_temperature=np.full(len(mesh.nodes), 400.),
                                               expansion_coefficient=0.)
        reference = update_plane_stress([strain, strain, 0.], previous, material, .1)
        previous = reference.state
        np.testing.assert_allclose(stress, np.tile(reference.stress, (model.point_count, 1)), atol=1e-8)
        np.testing.assert_allclose([h.accumulated_strain for h in histories], previous.accumulated_strain, atol=1e-12)
    assert previous.accumulated_strain > 0
    assert model.material.initial_yield == base.initial_yield
    assert model.material.relaxation_time == base.relaxation_time


def test_temperature_dependent_hardening_requires_temperature_input():
    model = NonlinearModel(structured_mesh(1, 1), J2Material(), hardening="paper-temperature")
    state = model.initial_state()
    with pytest.raises(ValueError, match="requires nodal temperatures"):
        model.evaluate(state.displacement, state.history, .1)


@pytest.mark.parametrize("mode,property_temperature", [("paper-reference", 20.), ("paper-temperature", 400.)])
def test_paper_elastic_modes_match_constant_temperature_elastic_stress(mode, property_temperature):
    mesh = structured_mesh(2, 2)
    model = NonlinearModel(mesh, J2Material(), elasticity=mode)
    displacement = (1e-4 * mesh.nodes).ravel()
    _, _, _, stress = model.evaluate(displacement, model.initial_state().history, .1,
                                     nodal_temperature=np.full(len(mesh.nodes), 400.),
                                     expansion_coefficient=0.)
    material = with_temperature_dependent_elasticity(J2Material(), property_temperature)
    expected = Material(material.young, material.poisson).D @ np.array([1e-4, 1e-4, 0.])
    np.testing.assert_allclose(stress, np.tile(expected, (model.point_count, 1)), atol=1e-10)
    assert model.material.initial_yield == 250.
    assert model.material.saturation_increment == 100.


def test_temperature_dependent_global_tangent_holds_temperature_and_history_fixed():
    mesh = structured_mesh(2, 2)
    model = NonlinearModel(mesh, J2Material(), elasticity="paper-temperature")
    state = model.initial_state()
    temperature = 20+50*mesh.nodes[:, 0]
    displacement = (1e-5*mesh.nodes).ravel()
    direction = np.random.default_rng(4).normal(size=model.ndofs)
    kwargs = dict(nodal_temperature=temperature, expansion_coefficient=0.)
    _, tangent, _, _ = model.evaluate(displacement, state.history, .1, **kwargs)
    step = 1e-8
    plus = model.evaluate(displacement+step*direction, state.history, .1, **kwargs)[0]
    minus = model.evaluate(displacement-step*direction, state.history, .1, **kwargs)[0]
    np.testing.assert_allclose((plus-minus)/(2*step), tangent@direction, rtol=1e-9, atol=1e-5)


def test_temperature_dependent_mode_requires_temperature_input():
    model = NonlinearModel(structured_mesh(1, 1), J2Material(), elasticity="paper-temperature")
    state = model.initial_state()
    with pytest.raises(ValueError, match="requires nodal temperatures"):
        model.evaluate(state.displacement, state.history, .1)


def test_cut_plate_mild_heating_matches_linear_thermoelasticity():
    mesh = structured_mesh(4, 4, radius=1.)
    x, y = mesh.nodes.T
    boundary = {int(i): 20. for i in np.flatnonzero(np.isclose(x, 0) & (y >= 1))}
    boundary.update({int(i): 120. for i in np.flatnonzero(np.isclose(x, 4))})
    heat = solve_stationary_heat(mesh, .045, 1., boundary)
    constraints = symmetry_constraints(mesh)
    linear = solve_thermoelastic(mesh, Material(), heat.temperature, constraints)
    model = NonlinearModel(mesh, J2Material())
    state = model.initial_state()
    for scale in np.linspace(.1, 1., 10):
        result = solve_increment(model, state, constraints, .1,
                                 nodal_temperature=20+scale*(heat.temperature-20))
        state = result.state
    np.testing.assert_allclose(result.state.displacement, linear.displacement, atol=1e-12)
    np.testing.assert_allclose(result.stress[model.physical], linear.stress[linear.physical], atol=1e-7)
    assert all(h.accumulated_strain == 0 for h in result.state.history)


def test_cut_plate_strong_heating_creates_physical_history_only():
    mesh = structured_mesh(4, 4, radius=1.)
    x, y = mesh.nodes.T
    boundary = {int(i): 20. for i in np.flatnonzero(np.isclose(x, 0) & (y >= 1))}
    boundary.update({int(i): 620. for i in np.flatnonzero(np.isclose(x, 4))})
    heat = solve_stationary_heat(mesh, .045, 1., boundary)
    model = NonlinearModel(mesh, J2Material(relaxation_time=100))
    previous = model.initial_state()
    state = previous
    for scale in np.linspace(.05, 1., 20):
        result = solve_increment(model, state, symmetry_constraints(mesh), .05,
                                 nodal_temperature=20+scale*(heat.temperature-20))
        state = result.state
    plastic = np.array([h.accumulated_strain for h in result.state.history])
    assert plastic[model.physical].max() > 0
    np.testing.assert_array_equal(plastic[~model.physical], 0.)
    assert all(h.accumulated_strain == 0 for h in previous.history)
    assert result.free_residual_norm < 1e-5


def test_restrained_thermal_history_matches_point_and_relaxes_during_hold():
    snapshots, stress_error, force_error, residual = run_check()
    assert stress_error < 1e-8 and force_error < 1e-8 and residual < 1e-8
    _, _, ramp_stress, ramp_p = snapshots[0]
    _, _, hold_stress, hold_p = snapshots[1]
    assert ramp_stress < hold_stress < 0.
    assert 0 < ramp_p < hold_p


def test_uniform_heating_in_nonlinear_model_allows_stress_free_expansion():
    mesh = structured_mesh(2, 2)
    model = NonlinearModel(mesh, J2Material())
    previous = model.initial_state()
    result = solve_increment(model, previous, symmetry_constraints(mesh), .1,
                             nodal_temperature=np.full(len(mesh.nodes), 70.))
    np.testing.assert_allclose(result.state.displacement.reshape(-1, 2), .0006*mesh.nodes, atol=1e-14)
    np.testing.assert_allclose(result.stress, 0., atol=1e-8)
    assert all(h.accumulated_strain == 0 for h in result.state.history)
    assert all(h.accumulated_strain == 0 for h in previous.history)
    np.testing.assert_array_equal(previous.displacement, 0.)


def test_nonlinear_thermal_input_requires_one_finite_temperature_per_node():
    mesh = structured_mesh(1, 1)
    model = NonlinearModel(mesh, J2Material())
    with pytest.raises(ValueError, match="one finite temperature"):
        solve_increment(model, model.initial_state(), symmetry_constraints(mesh), .1,
                        nodal_temperature=[np.nan]*len(mesh.nodes))
