"""Follow-along check: heat solution -> interpolated T -> free plate expansion."""

from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.solver import Material, structured_mesh, symmetry_constraints
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermoelastic import solve_thermoelastic

# Edit this small checkpoint here. The integrated playground comes later.
LENGTH, HEIGHT = 10., 2.  # mm
TEMPERATURE, REFERENCE_TEMPERATURE = 70., 20.  # degrees Celsius
EXPANSION_COEFFICIENT = 12e-6  # per K


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("free", "restrained"), default="free")
    args = parser.parse_args()
    mesh = structured_mesh(5, 2, length=LENGTH, height=HEIGHT)
    material = Material()
    left = np.flatnonzero(np.isclose(mesh.nodes[:, 0], 0.))
    right = np.flatnonzero(np.isclose(mesh.nodes[:, 0], LENGTH))
    heat = solve_stationary_heat(mesh, .045, material.thickness,
                                {int(i): TEMPERATURE for i in np.r_[left, right]})
    # ux=0 on x=0 and uy=0 on y=0 allow uniform expansion about the origin.
    # These rollers remove rigid motion without blocking the exact expansion.
    constraints = symmetry_constraints(mesh)
    top = np.flatnonzero(np.isclose(mesh.nodes[:, 1], HEIGHT))
    if args.case == "restrained":
        # Normal rollers on all four edges block in-plane expansion while
        # leaving tangential motion and interior displacement DOFs free.
        constraints.update({2 * int(i): 0. for i in right})
        constraints.update({2 * int(i) + 1: 0. for i in top})
    result = solve_thermoelastic(
        mesh, material, heat.temperature, constraints,
        reference_temperature=REFERENCE_TEMPERATURE,
        expansion_coefficient=EXPANSION_COEFFICIENT)
    expansion = EXPANSION_COEFFICIENT * (TEMPERATURE - REFERENCE_TEMPERATURE)
    exact = expansion * mesh.nodes if args.case == "free" else np.zeros_like(mesh.nodes)
    exact_stress = (np.zeros(3) if args.case == "free" else
                    np.array([-material.young * expansion / (1 - material.poisson)] * 2 + [0.]))
    displacement = result.displacement.reshape(-1, 2)
    free = np.setdiff1d(np.arange(len(result.residual)), list(constraints))
    print(f"Case: {args.case}")
    print(f"Temperature range at integration points: {result.temperature.min():.6f} to "
          f"{result.temperature.max():.6f} deg C")
    print(f"Right-edge ux: {displacement[right, 0].mean():.9f} mm; exact: {exact[right, 0].mean():.9f} mm")
    print(f"Top-edge uy: {displacement[top, 1].mean():.9f} mm; exact: {exact[top, 1].mean():.9f} mm")
    print(f"Maximum nodal displacement error: {np.max(np.abs(displacement - exact)):.3e} mm")
    print(f"Maximum absolute stress component: {np.max(np.abs(result.stress)):.3e} MPa")
    print(f"Mean stress [xx, yy, xy]: {result.stress.mean(axis=0)} MPa; exact: {exact_stress}")
    print(f"Maximum stress error: {np.max(np.abs(result.stress - exact_stress)):.3e} MPa")
    if args.case == "restrained":
        print(f"Right-edge x reaction: {result.residual[2 * right].sum():.6f} N; "
              f"exact: {exact_stress[0] * HEIGHT * material.thickness:.6f} N")
        print(f"Top-edge y reaction: {result.residual[2 * top + 1].sum():.6f} N; "
              f"exact: {exact_stress[1] * LENGTH * material.thickness:.6f} N")
    print(f"Maximum free-DOF residual: {np.max(np.abs(result.residual[free])):.3e} N")
    print("Constant properties, elastic plane stress, one-way coupling; no mechanical loads.")


if __name__ == "__main__":
    main()
