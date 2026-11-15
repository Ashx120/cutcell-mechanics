"""Free expansion under a linear temperature gradient: analytical convergence check."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.solver import Material, structured_mesh
from immersed_elasticity.thermal import solve_stationary_heat
from immersed_elasticity.thermoelastic import solve_thermoelastic

LENGTH, HEIGHT = 10., 2.  # mm
LEFT_TEMPERATURE, RIGHT_TEMPERATURE = 20., 120.  # degrees Celsius
REFERENCE_TEMPERATURE = 20.
EXPANSION_COEFFICIENT = 12e-6
MESHES = ((5, 2), (10, 4), (20, 8))


def main():
    print("Linear temperature gradient; constant properties; mechanically free edges.")
    print("Fix origin ux,uy and bottom-right uy to remove rigid motion only.")
    print(" mesh     max T error [K]    max nodal u error [mm]    RMS stress component [MPa]")
    for nx, ny in MESHES:
        mesh = structured_mesh(nx, ny, length=LENGTH, height=HEIGHT)
        x, y = mesh.nodes.T
        left = np.flatnonzero(np.isclose(x, 0.))
        right = np.flatnonzero(np.isclose(x, LENGTH))
        temperatures = {int(i): LEFT_TEMPERATURE for i in left}
        temperatures.update({int(i): RIGHT_TEMPERATURE for i in right})
        heat = solve_stationary_heat(mesh, .045, 1., temperatures)
        origin = int(np.flatnonzero(np.isclose(x, 0.) & np.isclose(y, 0.))[0])
        anchor = int(np.flatnonzero(np.isclose(x, LENGTH) & np.isclose(y, 0.))[0])
        result = solve_thermoelastic(
            mesh, Material(), heat.temperature,
            {2 * origin: 0., 2 * origin + 1: 0., 2 * anchor + 1: 0.},
            reference_temperature=REFERENCE_TEMPERATURE,
            expansion_coefficient=EXPANSION_COEFFICIENT)
        gradient = (RIGHT_TEMPERATURE - LEFT_TEMPERATURE) / LENGTH
        offset = LEFT_TEMPERATURE - REFERENCE_TEMPERATURE
        # sym(grad u) = alpha*(offset + gradient*x)*I; engineering shear=0.
        # This compatible quadratic displacement has exactly zero elastic stress.
        exact_u = EXPANSION_COEFFICIENT * np.column_stack(
            (offset*x + .5*gradient*(x*x-y*y), (offset+gradient*x)*y))
        temperature_error = np.max(np.abs(heat.temperature - (LEFT_TEMPERATURE + gradient*x)))
        displacement_error = np.max(np.abs(result.displacement.reshape(-1, 2) - exact_u))
        # Equal-area cells and identical Gauss weights make this an area RMS.
        stress_rms = np.sqrt(np.mean(result.stress**2))
        print(f"{nx:3d}x{ny:<3d} {temperature_error:18.3e} {displacement_error:27.3e} {stress_rms:29.6f}")
    print("Exact continuum stress is zero. Q4 displacement approximation leaves a mesh error.")
    print("A temperature gradient alone does not necessarily produce thermal stress.")


if __name__ == "__main__":
    main()
