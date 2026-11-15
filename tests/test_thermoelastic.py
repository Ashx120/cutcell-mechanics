import numpy as np
import pytest

from immersed_elasticity.solver import Material, structured_mesh, symmetry_constraints
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermoelastic import solve_thermoelastic


def test_heat_solution_drives_exact_stress_free_uniform_expansion():
    mesh = structured_mesh(5, 2, length=10, height=2)
    boundary = np.flatnonzero(np.isclose(mesh.nodes[:, 0], 0) | np.isclose(mesh.nodes[:, 0], 10))
    heat = solve_stationary_heat(mesh, .045, 1., {int(i): 70. for i in boundary})
    result = solve_thermoelastic(mesh, Material(), heat.temperature, symmetry_constraints(mesh))
    np.testing.assert_allclose(result.displacement.reshape(-1, 2), .0006 * mesh.nodes, atol=1e-15)
    np.testing.assert_allclose(result.strain, np.tile([.0006, .0006, 0], (len(result.points), 1)), atol=1e-15)
    np.testing.assert_allclose(result.stress, 0., atol=1e-9)
    np.testing.assert_allclose(result.residual, 0., atol=1e-9)


def test_full_inplane_restraint_recovers_biaxial_compression():
    mesh = structured_mesh(3, 2)
    right = np.flatnonzero(np.isclose(mesh.nodes[:, 0], mesh.length))
    top = np.flatnonzero(np.isclose(mesh.nodes[:, 1], mesh.height))
    constraints = symmetry_constraints(mesh)
    constraints.update({2 * int(i): 0. for i in right})
    constraints.update({2 * int(i) + 1: 0. for i in top})
    result = solve_thermoelastic(mesh, Material(), np.full(len(mesh.nodes), 70.),
                                constraints)
    np.testing.assert_allclose(result.stress, np.tile([-180., -180., 0.], (len(result.points), 1)), atol=1e-10)
    np.testing.assert_allclose(result.displacement, 0., atol=1e-15)
    np.testing.assert_allclose(result.residual[2 * right].sum(), -180 * mesh.height, atol=1e-9)
    np.testing.assert_allclose(result.residual[2 * top + 1].sum(), -180 * mesh.length, atol=1e-9)
    free = np.setdiff1d(np.arange(len(result.residual)), list(constraints))
    np.testing.assert_allclose(result.residual[free], 0., atol=1e-9)


def test_temperature_interpolates_bilinear_nodal_field_at_quadrature_points():
    mesh = structured_mesh(3, 2)
    def field(points):
        x, y = points.T
        return 20 + 3*x + 2*y + x*y
    result = solve_thermoelastic(mesh, Material(), field(mesh.nodes), symmetry_constraints(mesh))
    np.testing.assert_allclose(result.temperature, field(result.points), atol=1e-12)


def test_free_linear_temperature_gradient_converges_to_compatible_zero_stress_field():
    errors, stresses = [], []
    for nx, ny in ((5, 2), (10, 4), (20, 8)):
        mesh = structured_mesh(nx, ny, length=10, height=2)
        x, y = mesh.nodes.T
        boundary = {int(i): 20. + 10*x[i] for i in np.flatnonzero(np.isclose(x, 0) | np.isclose(x, 10))}
        heat = solve_stationary_heat(mesh, .045, 1., boundary)
        anchor = int(np.flatnonzero(np.isclose(x, 10) & np.isclose(y, 0))[0])
        result = solve_thermoelastic(mesh, Material(), heat.temperature, {0: 0., 1: 0., 2*anchor+1: 0.})
        exact = 12e-6 * np.column_stack((5*(x*x-y*y), 10*x*y))
        np.testing.assert_allclose(heat.temperature, 20 + 10*x, atol=1e-10)
        errors.append(np.max(np.abs(result.displacement.reshape(-1, 2) - exact)))
        stresses.append(np.sqrt(np.mean(result.stress**2)))
    assert errors[1] < .6 * errors[0] and errors[2] < .6 * errors[1]
    assert stresses[1] < .7 * stresses[0] and stresses[2] < .7 * stresses[1]


def test_cut_plate_thermal_response_is_linear_and_in_equilibrium():
    mesh = structured_mesh(8, 8, radius=1.)
    x, y = mesh.nodes.T
    boundary = {int(i): 20. for i in np.flatnonzero(np.isclose(x, 0) & (y >= 1))}
    boundary.update({int(i): 120. for i in np.flatnonzero(np.isclose(x, 4))})
    heat = solve_stationary_heat(mesh, .045, 1., boundary)
    constraints = symmetry_constraints(mesh)
    result = solve_thermoelastic(mesh, Material(), heat.temperature, constraints)
    doubled = solve_thermoelastic(mesh, Material(), 20+2*(heat.temperature-20), constraints)
    np.testing.assert_allclose(doubled.displacement, 2*result.displacement, atol=1e-13)
    np.testing.assert_allclose(doubled.stress, 2*result.stress, atol=1e-8)
    free = np.setdiff1d(np.arange(len(result.residual)), list(constraints))
    np.testing.assert_allclose(result.residual[free], 0, atol=1e-9)
    np.testing.assert_allclose(result.residual.reshape(-1, 2).sum(axis=0), 0, atol=1e-9)
    np.testing.assert_array_equal(result.physical, heat.physical)
    assert np.max(np.abs(result.stress[result.physical])) > 1.


@pytest.mark.parametrize("temperature", [[70.], [np.nan] * 4])
def test_invalid_nodal_temperatures_rejected(temperature):
    mesh = structured_mesh(1, 1)
    with pytest.raises(ValueError, match="one finite temperature"):
        solve_thermoelastic(mesh, Material(), temperature, symmetry_constraints(mesh))
