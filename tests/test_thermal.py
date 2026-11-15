import numpy as np
import pytest

from immersed_elasticity.solver import structured_mesh
from immersed_elasticity.thermal import solve_stationary_heat


def test_rectangular_conduction_matches_exact_profile_flux_and_heat_balance():
    mesh = structured_mesh(5, 2, length=10, height=2)
    x = mesh.nodes[:, 0]
    left, right = np.flatnonzero(np.isclose(x, 0)), np.flatnonzero(np.isclose(x, 10))
    prescribed = {int(i): 20. for i in left}
    prescribed.update({int(i): 120. for i in right})
    result = solve_stationary_heat(mesh, .045, 1., prescribed)
    np.testing.assert_allclose(result.temperature, 20 + 10*x, atol=1e-11)
    np.testing.assert_allclose(result.flux[:, 0], -.45, atol=1e-12)
    np.testing.assert_allclose(result.flux[:, 1], 0, atol=1e-12)
    np.testing.assert_allclose(result.heat_residual[left].sum(), -.9, atol=1e-12)
    np.testing.assert_allclose(result.heat_residual[right].sum(), .9, atol=1e-12)
    free = np.setdiff1d(np.arange(len(x)), list(prescribed))
    np.testing.assert_allclose(result.heat_residual[free], 0, atol=1e-12)
    doubled = solve_stationary_heat(mesh, .045, 2., prescribed)
    np.testing.assert_allclose(doubled.temperature, result.temperature)
    np.testing.assert_allclose(doubled.heat_residual, 2 * result.heat_residual, atol=1e-12)


def test_uniform_temperature_has_zero_flux():
    mesh = structured_mesh(3, 2)
    boundary = {int(i): 75. for i in np.flatnonzero(np.isclose(mesh.nodes[:, 0], 0))}
    result = solve_stationary_heat(mesh, .045, 1., boundary)
    np.testing.assert_allclose(result.temperature, 75., atol=1e-11)
    np.testing.assert_allclose(result.flux, 0, atol=1e-12)


@pytest.mark.parametrize("conductivity,thickness", [(0, 1), (-1, 1), (.045, 0), (np.nan, 1)])
def test_invalid_thermal_properties(conductivity, thickness):
    with pytest.raises(ValueError, match="positive finite"):
        solve_stationary_heat(structured_mesh(1, 1), conductivity, thickness, {0: 20})


def test_thermal_needs_a_temperature_constraint():
    with pytest.raises(ValueError, match="prescribed temperature"):
        solve_stationary_heat(structured_mesh(1, 1), .045, 1., {})


def test_cut_plate_heat_balance_linearity_and_conductivity_scaling():
    mesh = structured_mesh(8, 8, radius=1.)
    x, y = mesh.nodes.T
    left = np.flatnonzero(np.isclose(x, 0) & (y >= 1))
    right = np.flatnonzero(np.isclose(x, 4))
    constraints = {int(i): 20. for i in left}
    constraints.update({int(i): 120. for i in right})
    result = solve_stationary_heat(mesh, .045, 1., constraints)
    free = np.setdiff1d(np.arange(len(x)), list(constraints))
    np.testing.assert_allclose(result.heat_residual[free], 0, atol=1e-11)
    np.testing.assert_allclose(result.heat_residual[list(constraints)].sum(), 0, atol=1e-11)
    assert result.heat_residual[right].sum() > 0
    assert result.heat_residual[left].sum() < 0
    assert np.any(~result.physical) and np.any(result.physical)
    doubled = solve_stationary_heat(mesh, .09, 1., constraints)
    np.testing.assert_allclose(doubled.temperature, result.temperature, atol=1e-9)
    np.testing.assert_allclose(doubled.heat_residual, 2*result.heat_residual, atol=1e-11)
    uniform = solve_stationary_heat(mesh, .045, 1., {i: 75. for i in constraints})
    np.testing.assert_allclose(uniform.temperature, 75., atol=1e-8)
    np.testing.assert_allclose(uniform.flux, 0., atol=1e-10)
