import numpy as np
import pytest

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.solver import structured_mesh
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


def test_uniform_shear_equilibrium_matches_point_through_loading_and_hold():
    mesh = structured_mesh(2, 2, length=2, height=1)
    material = J2Material(relaxation_time=100)
    model = NonlinearModel(mesh, material)
    state, point = model.initial_state(), MaterialState()
    x, y = mesh.nodes.T
    boundary = np.flatnonzero(np.isclose(x, 0) | np.isclose(x, 2) | np.isclose(y, 0) | np.isclose(y, 1))
    for gamma in (.003, .006, .006, .002):
        constraints = {2*n+c: gamma*y[n] if c == 0 else 0. for n in boundary for c in (0, 1)}
        previous_displacement = state.displacement.copy()
        previous_p = [h.accumulated_strain for h in state.history]
        result = solve_increment(model, state, constraints, .1)
        reference = update_plane_stress([0, 0, gamma], point, material, .1)
        np.testing.assert_allclose(result.stress, np.tile(reference.stress, (model.point_count, 1)), atol=1e-6)
        np.testing.assert_allclose(result.state.displacement[0::2], gamma*y, atol=1e-11)
        np.testing.assert_allclose([h.accumulated_strain for h in result.state.history],
                                   reference.state.accumulated_strain, atol=1e-11)
        np.testing.assert_allclose(result.residual.reshape(-1, 2).sum(axis=0), 0, atol=1e-7)
        np.testing.assert_array_equal(state.displacement, previous_displacement)
        assert [h.accumulated_strain for h in state.history] == previous_p
        assert result.free_residual_norm < 1e-6
        state, point = result.state, reference.state


def test_failed_increment_leaves_previous_state_unchanged():
    mesh = structured_mesh(2, 2, length=2, height=1)
    model = NonlinearModel(mesh, J2Material())
    state = model.initial_state()
    constraints = {2*i: .01 for i in np.flatnonzero(np.isclose(mesh.nodes[:, 1], 1))}
    with pytest.raises(RuntimeError, match="not committed"):
        solve_increment(model, state, constraints, .1, max_iterations=0)
    np.testing.assert_array_equal(state.displacement, np.zeros(model.ndofs))
    assert all(h.accumulated_strain == 0 for h in state.history)


@pytest.mark.parametrize("method", ["quadtree", "moment-fitting"])
def test_cut_plate_elastic_limit_matches_v1_solver(method):
    from immersed_elasticity.solver import Material, assemble, solve, symmetry_constraints

    mesh = structured_mesh(4, 4, radius=1)
    material = J2Material(initial_yield=1e9)
    model = NonlinearModel(mesh, material, depth=2, method=method)
    constraints = symmetry_constraints(mesh)
    right = np.flatnonzero(np.isclose(mesh.nodes[:, 0], mesh.length))
    constraints.update({2*n: 1e-4 for n in right})
    result = solve_increment(model, model.initial_state(), constraints, .1)
    stiffness, _ = assemble(mesh, Material(), depth=2)
    exact, _ = solve(stiffness, np.zeros(model.ndofs), constraints)
    np.testing.assert_allclose(result.state.displacement, exact, atol=1e-12)


def test_fictitious_domain_remains_elastic_while_physical_material_yields():
    from immersed_elasticity.solver import symmetry_constraints

    mesh = structured_mesh(4, 4, radius=1)
    model = NonlinearModel(mesh, J2Material(relaxation_time=100), depth=2)
    constraints = symmetry_constraints(mesh)
    constraints.update({2*n: .012 for n in np.flatnonzero(np.isclose(mesh.nodes[:, 0], mesh.length))})
    result = solve_increment(model, model.initial_state(), constraints, .1)
    p = np.array([h.accumulated_strain for h in result.state.history])
    assert p[model.physical].max() > 0
    np.testing.assert_array_equal(p[~model.physical], 0)


def test_separate_region_fitting_reduces_physical_history_points():
    mesh = structured_mesh(4, 4, radius=1)
    material = J2Material()
    original = NonlinearModel(mesh, material, depth=3)
    fitted = NonlinearModel(mesh, material, depth=3, method="moment-fitting")
    assert fitted.point_count < original.point_count
    assert fitted.physical.sum() < original.physical.sum()
    for physical in (True, False):
        original_mass = sum(np.concatenate([w for _, _, w in original.elements])[original.physical == physical])
        fitted_mass = sum(np.concatenate([w for _, _, w in fitted.elements])[fitted.physical == physical])
        np.testing.assert_allclose(fitted_mass, original_mass, rtol=1e-10)


def test_nonlinear_fit_fallback_retains_original_rules(monkeypatch):
    import immersed_elasticity.nonlinear_solver as solver

    def fail(*args, **kwargs):
        raise RuntimeError("Failed fit")

    mesh = structured_mesh(4, 4, radius=1)
    original = NonlinearModel(mesh, J2Material(), depth=2)
    monkeypatch.setattr(solver, "compress_rule", fail)
    with pytest.warns(RuntimeWarning, match="fallback"):
        fitted = NonlinearModel(mesh, J2Material(), depth=2, method="moment-fitting")
    np.testing.assert_array_equal(fitted.points, original.points)
    np.testing.assert_array_equal(fitted.physical, original.physical)
    assert fitted.fallback_cells > 0
