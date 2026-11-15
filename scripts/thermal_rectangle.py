"""First v4 checkpoint: a rectangle with an exact stationary temperature profile."""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from immersed_elasticity.solver import structured_mesh
from immersed_elasticity.thermal import solve_stationary_heat

# Illustrative constant properties. Conductivity in W/(mm K), dimensions in mm.
GEOMETRY = dict(nx=5, ny=2, length=10., height=2., radius=0.)
CONDUCTIVITY = 0.045  # Equivalent to 45 W/(m K).
THICKNESS = 1.0
LEFT_TEMPERATURE = 20.0  # Degrees Celsius.
RIGHT_TEMPERATURE = 120.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/walkthrough/v4-thermal-rectangle"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        if not args.output.is_dir():
            parser.error("Output must be a directory.")
        if any(args.output.iterdir()) and not args.overwrite:
            parser.error("Output folder is nonempty. Choose a new folder or use --overwrite.")
    mesh = structured_mesh(**GEOMETRY)
    x = mesh.nodes[:, 0]
    left = np.flatnonzero(np.isclose(x, 0))
    right = np.flatnonzero(np.isclose(x, mesh.length))
    prescribed = {int(node): LEFT_TEMPERATURE for node in left}
    prescribed.update({int(node): RIGHT_TEMPERATURE for node in right})
    result = solve_stationary_heat(mesh, CONDUCTIVITY, THICKNESS, prescribed)
    exact = LEFT_TEMPERATURE + (RIGHT_TEMPERATURE - LEFT_TEMPERATURE) * x / mesh.length
    exact_flux = -CONDUCTIVITY * (RIGHT_TEMPERATURE - LEFT_TEMPERATURE) / mesh.length
    free = np.setdiff1d(np.arange(len(mesh.nodes)), list(prescribed))
    print(f"Nodal temperature range: {result.temperature.min():.3f} to {result.temperature.max():.3f} deg C")
    print(f"Maximum nodal temperature error: {abs(result.temperature - exact).max():.3e} K")
    print(f"Mean q_x: {result.flux[:, 0].mean():+.6f} W/mm^2; exact: {exact_flux:+.6f} W/mm^2")
    print(f"Left reservoir heat input: {result.heat_residual[left].sum():+.6f} W")
    print(f"Right reservoir heat input: {result.heat_residual[right].sum():+.6f} W")
    print(f"Maximum free-node heat residual: {abs(result.heat_residual[free]).max():.3e} W")
    print(f"Net reservoir heat input: {result.heat_residual[list(prescribed)].sum():.3e} W")
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "temperature.csv", np.column_stack((mesh.nodes, result.temperature, exact)),
               delimiter=",", comments="", header="x_mm,y_mm,temperature_C,exact_temperature_C")
    np.savetxt(args.output / "flux.csv", np.column_stack((result.points, result.flux)),
               delimiter=",", comments="", header="x_mm,y_mm,qx_W_per_mm2,qy_W_per_mm2")
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), layout="constrained")
    field = axes[0].tricontourf(*mesh.nodes.T, result.temperature, levels=20, cmap="inferno")
    axes[0].set(xlabel="x [mm]", ylabel="y [mm]", aspect="equal", title="Stationary temperature field")
    fig.colorbar(field, ax=axes[0], label="Temperature [deg C]")
    midline = np.isclose(mesh.nodes[:, 1], mesh.height / 2)
    axes[1].plot([0, mesh.length], [LEFT_TEMPERATURE, RIGHT_TEMPERATURE], label="Analytical profile")
    axes[1].plot(x[midline], result.temperature[midline], "o", label="FEM midline nodes")
    axes[1].set(xlabel="x [mm]", ylabel="Temperature [deg C]")
    axes[1].legend()
    axes[1].grid(alpha=.25)
    fig.savefig(args.output / "temperature.png", dpi=160)
    plt.close(fig)
    print(f"Saved temperature/flux CSV files and temperature.png to {args.output}")


if __name__ == "__main__":
    main()
