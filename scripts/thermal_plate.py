"""Stationary heat flow around an approximately insulated quarter-circle hole."""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from immersed_elasticity.solver import element_kinematics, structured_mesh
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.quadrature import cut_rule

GEOMETRY = dict(length=4., height=4., radius=1.)
CONDUCTIVITY, THICKNESS = .045, 1.
LEFT_TEMPERATURE, RIGHT_TEMPERATURE = 20., 120.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, default=16)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--alpha", type=float, default=1e-6)
    parser.add_argument("--output", type=Path, default=Path("results/walkthrough/v4-thermal-plate"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or
                                 (any(args.output.iterdir()) and not args.overwrite)):
        parser.error("Choose an empty output directory or use --overwrite.")
    mesh = structured_mesh(args.mesh, args.mesh, **GEOMETRY)
    x, y = mesh.nodes.T
    left = np.flatnonzero(np.isclose(x, 0.) & (y >= mesh.radius))
    right = np.flatnonzero(np.isclose(x, mesh.length))
    prescribed = {int(i): LEFT_TEMPERATURE for i in left}
    prescribed.update({int(i): RIGHT_TEMPERATURE for i in right})
    result = solve_stationary_heat(mesh, CONDUCTIVITY, THICKNESS, prescribed,
                                   depth=args.depth, alpha=args.alpha)
    free = np.setdiff1d(np.arange(len(x)), list(prescribed))
    cold, hot = result.heat_residual[left].sum(), result.heat_residual[right].sum()
    print(f"Active cells: {len(mesh.cells)}; temperature DOFs: {len(x)}; integration points: {len(result.points)}")
    print(f"Fictitious conductivity factor: {args.alpha:g}")
    print(f"Cold-edge reservoir heat input: {cold:+.9f} W")
    print(f"Hot-edge reservoir heat input: {hot:+.9f} W")
    print(f"Net reservoir heat input: {cold+hot:.3e} W")
    print(f"Maximum free-node heat residual: {np.max(np.abs(result.heat_residual[free])):.3e} W")
    print("Top/bottom boundaries insulated; hole insulation approximated by weak fictitious conductivity.")
    # Reuse the heat solver's point ordering; only physical samples are plotted.
    sample_temperature = []
    for cell in mesh.cells:
        xy = mesh.nodes[cell]
        # Each cell occupies a contiguous block of quadrature points.
        points, _ = cut_rule((xy[0, 0], xy[1, 0], xy[0, 1], xy[3, 1]), mesh.radius, args.depth, args.alpha)
        shape, _ = element_kinematics(xy, points)
        sample_temperature.extend(shape @ result.temperature[cell])
    sample_temperature = np.asarray(sample_temperature)
    mask = result.physical
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "temperature.csv", np.column_stack((mesh.nodes, result.temperature)),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C")
    np.savetxt(args.output / "physical_samples.csv",
               np.column_stack((result.points[mask], sample_temperature[mask], result.flux[mask])),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C,qx_W_per_mm2,qy_W_per_mm2")
    fig, ax = plt.subplots(figsize=(7, 6), layout="constrained")
    field = ax.scatter(*result.points[mask].T, c=sample_temperature[mask], s=9, cmap="inferno",
                       vmin=LEFT_TEMPERATURE, vmax=RIGHT_TEMPERATURE)
    angle = np.linspace(0, np.pi/2, 100)
    ax.plot(mesh.radius*np.cos(angle), mesh.radius*np.sin(angle), color="black")
    ax.set(xlim=(0, mesh.length), ylim=(0, mesh.height), aspect="equal",
           xlabel="x [mm]", ylabel="y [mm]", title="Temperature at physical integration points")
    fig.colorbar(field, ax=ax, label="Temperature [deg C]")
    fig.savefig(args.output / "temperature.png", dpi=160)
    plt.close(fig)
    print(f"Saved nodal temperatures, physical samples, and plot to {args.output}")


if __name__ == "__main__":
    main()
