"""Stationary temperature -> thermal expansion -> elastic quarter-plate response."""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from immersed_elasticity.solver import Material, structured_mesh
from immersed_elasticity.boundary_conditions import thermal_constraints
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermoelastic import solve_thermoelastic

# This checkpoint's inputs are all here; integrated playground support follows.
GEOMETRY = dict(length=4., height=4., radius=1.)
MATERIAL = Material(young=210000., poisson=.3, thickness=1.)
CONDUCTIVITY = .045  # W/(mm K)
REFERENCE_TEMPERATURE, LEFT_TEMPERATURE, RIGHT_TEMPERATURE = 20., 20., 120.
EXPANSION_COEFFICIENT = 12e-6  # per K


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, default=16)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--alpha", type=float, default=1e-6)
    parser.add_argument("--right-temperature", type=float, default=RIGHT_TEMPERATURE,
                        help="Hot-edge temperature in degrees Celsius (default: 120).")
    parser.add_argument("--supports", choices=("symmetry", "minimal"), default="symmetry",
                        help="Edge rollers or three constraints removing rigid motion only.")
    parser.add_argument("--output", type=Path, default=Path("results/walkthrough/v4-thermoelastic-plate"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if not np.isfinite(args.right_temperature):
        parser.error("--right-temperature must be finite")
    if args.output.exists() and (not args.output.is_dir() or
                                 (any(args.output.iterdir()) and not args.overwrite)):
        parser.error("Choose an empty output directory or use --overwrite.")
    mesh = structured_mesh(args.mesh, args.mesh, **GEOMETRY)
    x, y = mesh.nodes.T
    left = np.flatnonzero(np.isclose(x, 0.) & (y >= mesh.radius))
    right = np.flatnonzero(np.isclose(x, mesh.length))
    prescribed = {int(i): LEFT_TEMPERATURE for i in left}
    prescribed.update({int(i): args.right_temperature for i in right})
    heat = solve_stationary_heat(mesh, CONDUCTIVITY, MATERIAL.thickness, prescribed,
                                depth=args.depth, alpha=args.alpha)
    constraints = thermal_constraints(mesh, args.supports)
    result = solve_thermoelastic(
        mesh, MATERIAL, heat.temperature, constraints, depth=args.depth, alpha=args.alpha,
        reference_temperature=REFERENCE_TEMPERATURE, expansion_coefficient=EXPANSION_COEFFICIENT)
    free = np.setdiff1d(np.arange(len(result.residual)), list(constraints))
    displacement = result.displacement.reshape(-1, 2)
    mask = result.physical
    stress = result.stress[mask]
    von_mises = np.sqrt(stress[:, 0]**2 - stress[:, 0]*stress[:, 1] + stress[:, 1]**2 + 3*stress[:, 2]**2)
    print(f"Active cells: {len(mesh.cells)}; mechanical DOFs: {len(result.displacement)}")
    if args.supports == "symmetry":
        print("Supports: ux=0 on x=0; uy=0 on y=0. Right/top edges mechanically unloaded.")
    else:
        print("Supports: top-left ux=uy=0; top-right uy=0; rigid motion removed only.")
    print(f"Reservoir temperatures: left={LEFT_TEMPERATURE:g} C; right={args.right_temperature:g} C; "
          f"stress-free reference={REFERENCE_TEMPERATURE:g} C")
    print(f"Hot-edge heat input: {heat.heat_residual[right].sum():.9f} W")
    print(f"Mean right-edge ux: {displacement[right, 0].mean():.9f} mm")
    print(f"Maximum physical-point von Mises stress: {von_mises.max():.6f} MPa")
    peak = result.points[mask][int(np.argmax(von_mises))]
    print(f"Sampled peak location: x={peak[0]:.6f}, y={peak[1]:.6f} mm")
    print(f"Maximum constrained-DOF reaction: {np.max(np.abs(result.residual[list(constraints)])):.3e} N")
    print(f"Maximum free mechanical residual: {np.max(np.abs(result.residual[free])):.3e} N")
    print(f"Net support force [x, y]: {result.residual.reshape(-1, 2).sum(axis=0)} N")
    print("Elastic diagnostic only: no yielding or temperature-dependent properties.")
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "nodes.csv", np.column_stack((mesh.nodes, heat.temperature, displacement)),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C,ux_mm,uy_mm")
    np.savetxt(args.output / "physical_samples.csv",
               np.column_stack((result.points[mask], result.temperature[mask], stress, von_mises)),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C,sigma_xx_MPa,sigma_yy_MPa,tau_xy_MPa,von_mises_MPa")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for ax, values, title, label in zip(axes, (result.temperature[mask], von_mises),
                                      ("Temperature", "Elastic von Mises stress"), ("deg C", "MPa")):
        field = ax.scatter(*result.points[mask].T, c=values, s=10, cmap="inferno")
        angle = np.linspace(0, np.pi/2, 100)
        ax.plot(mesh.radius*np.cos(angle), mesh.radius*np.sin(angle), color="black")
        ax.set(xlim=(0, mesh.length), ylim=(0, mesh.height), aspect="equal",
               xlabel="x [mm]", ylabel="y [mm]", title=title)
        fig.colorbar(field, ax=ax, label=label)
    fig.savefig(args.output / "fields.png", dpi=160)
    plt.close(fig)
    print(f"Saved nodal results, physical samples, and fields.png to {args.output}")


if __name__ == "__main__":
    main()
