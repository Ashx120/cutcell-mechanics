"""Compare a uniform-shear FEM rectangle with a material-point history."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.solver import structured_mesh
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


def main():
    mesh = structured_mesh(2, 2, length=2, height=1)
    material = J2Material(relaxation_time=100)
    model = NonlinearModel(mesh, material, thickness=1)
    state = model.initial_state()
    point_history = MaterialState()
    x, y = mesh.nodes.T
    boundary = np.flatnonzero(np.isclose(x, 0) | np.isclose(x, 2) | np.isclose(y, 0) | np.isclose(y, 1))
    top = np.flatnonzero(np.isclose(y, 1))
    times = np.linspace(0, 4, 81)
    gammas = np.interp(times, [0, 1, 3, 4], [0, .01, .01, 0])
    largest_error = largest_residual = 0.
    maximum_iterations = 0
    for step in range(1, len(times)):
        gamma = gammas[step]
        dt = times[step] - times[step-1]
        constraints = {}
        for node in boundary:
            constraints[2 * node] = gamma * y[node]
            constraints[2 * node + 1] = 0.
        result = solve_increment(model, state, constraints, dt)
        state = result.state  # Commit only after global equilibrium converges.
        reference = update_plane_stress([0., 0., gamma], point_history, material, dt)
        point_history = reference.state
        top_force = result.residual[2 * top].sum()
        expected_force = reference.stress[2] * mesh.length  # thickness = 1 mm
        largest_error = max(largest_error, abs(top_force - expected_force))
        largest_residual = max(largest_residual, result.free_residual_norm)
        maximum_iterations = max(maximum_iterations, result.iterations)
        if np.any(np.isclose(times[step], [1, 3, 4])):
            print(f"t={times[step]:g} s: FEM top force={top_force:+.6f} N; point prediction={expected_force:+.6f} N")
    print(f"Largest force difference across history: {largest_error:.3e} N")
    print(f"Largest free-DOF residual norm: {largest_residual:.3e} N")
    print(f"Maximum global Newton corrections per step: {maximum_iterations}")
    if largest_error > 1e-5:
        raise RuntimeError("Rectangle disagrees with material-point prediction")


if __name__ == "__main__":
    main()
