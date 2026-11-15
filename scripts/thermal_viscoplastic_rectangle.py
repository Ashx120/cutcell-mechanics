"""Restrained thermal ramp and hold: FEM against a material-point history."""

from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.solver import structured_mesh, symmetry_constraints
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermal_expansion import update_thermal_plane_stress
from immersed_elasticity.viscoplasticity import J2Material, MaterialState

REFERENCE_TEMPERATURE, FINAL_TEMPERATURE = 20., 170.  # deg C
STEPS_PER_SECOND = 20


def run_check(steps_per_second=STEPS_PER_SECOND, *, return_history=False):
    if (not isinstance(steps_per_second, (int, np.integer))
            or isinstance(steps_per_second, bool) or steps_per_second < 1):
        raise ValueError("Require a positive integer steps_per_second.")
    mesh = structured_mesh(2, 2, length=2., height=1.)
    material = J2Material(relaxation_time=100.)
    model = NonlinearModel(mesh, material)
    state, point_history = model.initial_state(), MaterialState()
    x, y = mesh.nodes.T
    left = np.flatnonzero(np.isclose(x, 0.))
    right = np.flatnonzero(np.isclose(x, mesh.length))
    top = np.flatnonzero(np.isclose(y, mesh.height))
    constraints = symmetry_constraints(mesh)
    constraints.update({2*int(i): 0. for i in right})
    constraints.update({2*int(i)+1: 0. for i in top})
    # Solve a uniform final thermal field once. Scaling its rise gives uniform
    # equilibrium snapshots; no transient heat equation is claimed here.
    heat = solve_stationary_heat(mesh, .045, 1.,
                                {int(i): FINAL_TEMPERATURE for i in np.r_[left, right]})
    largest_stress_error = largest_force_error = largest_residual = 0.
    snapshots = []
    history = [(0., 0., 0.)]  # time, mean sigma_xx, maximum accumulated strain.
    dt = 1 / steps_per_second
    for step in range(1, 3*steps_per_second + 1):
        time = step * dt
        scale = min(time, 1.)
        temperature = REFERENCE_TEMPERATURE + scale*(heat.temperature-REFERENCE_TEMPERATURE)
        result = solve_increment(model, state, constraints, dt, nodal_temperature=temperature,
                                 reference_temperature=REFERENCE_TEMPERATURE)
        state = result.state
        point = update_thermal_plane_stress(
            [0., 0., 0.], point_history, material, dt,
            temperature=REFERENCE_TEMPERATURE+scale*(FINAL_TEMPERATURE-REFERENCE_TEMPERATURE),
            reference_temperature=REFERENCE_TEMPERATURE)
        point_history = point.state
        largest_stress_error = max(largest_stress_error, np.max(np.abs(result.stress-point.stress)))
        force = result.residual[2*right].sum()
        largest_force_error = max(largest_force_error, abs(force-point.stress[0]*mesh.height))
        largest_residual = max(largest_residual, result.free_residual_norm)
        history.append((time, float(result.stress[:, 0].mean()),
                        max(h.accumulated_strain for h in state.history)))
        if step in (steps_per_second, 3*steps_per_second):
            snapshots.append((time, float(temperature.mean()), float(result.stress[:, 0].mean()),
                              max(h.accumulated_strain for h in state.history)))
    summary = snapshots, largest_stress_error, largest_force_error, largest_residual
    return (*summary, np.asarray(history)) if return_history else summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps-per-second", type=int, default=STEPS_PER_SECOND)
    parser.add_argument("--compare-time-steps", action="store_true",
                        help="Compare dt, dt/2, dt/4 at shared history times.")
    args = parser.parse_args()
    if args.steps_per_second < 1:
        parser.error("--steps-per-second must be positive")
    print("0-1 s: uniform heating 20 -> 170 C; 1-3 s: temperature hold.")
    print("All four edges restrained in their normal direction; constant material properties.")
    print(f"Time step: {1 / args.steps_per_second:g} s")
    snapshots, stress_error, force_error, residual, history = run_check(
        args.steps_per_second, return_history=True)
    for time, temperature, stress, p in snapshots:
        print(f"t={time:g} s: T={temperature:.1f} C; sigma_xx={stress:.6f} MPa; accumulated strain={p:.8f}")
    print(f"Largest FEM/material-point stress difference: {stress_error:.3e} MPa")
    print(f"Largest right-force difference: {force_error:.3e} N")
    print(f"Maximum free residual norm: {residual:.3e} N")
    print("Prescribed equilibrium temperature snapshots; no transient conduction or plastic heating.")
    if args.compare_time_steps:
        for refinement in (2, 4):
            count = args.steps_per_second * refinement
            finer = run_check(count, return_history=True)[-1]
            shared = finer[::2]
            np.testing.assert_allclose(history[:, 0], shared[:, 0], atol=1e-12, rtol=0)
            differences = np.abs(shared[:, 1] - history[:, 1])
            index = int(np.argmax(differences))
            print(f"dt={2/count:g} vs {1/count:g} s: largest shared-time stress difference="
                  f"{differences[index]:.6f} MPa at t={history[index, 0]:.3f} s; "
                  f"largest accumulated-strain difference={np.max(np.abs(shared[:, 2]-history[:, 2])):.3e}")
            history = finer


if __name__ == "__main__":
    main()
