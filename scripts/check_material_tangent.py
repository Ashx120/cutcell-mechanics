"""Compare the analytical material tangent with central finite differences."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.viscoplasticity import (
    J2Material, MaterialState, engineering_strain_basis, stress_vector,
    update_material, update_with_tangent,
)


def main():
    material = J2Material(relaxation_time=100)
    dt = .01
    previous = MaterialState()
    for label, gamma in (("Elastic", .0001), ("Viscoplastic", .01)):
        strain = engineering_strain_basis(3) * gamma
        _, _, tangent = update_with_tangent(strain, previous, material, dt)
        for step in (1e-6, 1e-7, 1e-8):
            numerical = np.empty((6, 6))
            for column in range(6):
                perturbation = step * engineering_strain_basis(column)
                plus, _ = update_material(strain + perturbation, previous, material, dt)
                minus, _ = update_material(strain - perturbation, previous, material, dt)
                numerical[:, column] = (stress_vector(plus) - stress_vector(minus)) / (2 * step)
            error = np.linalg.norm(tangent - numerical) / np.linalg.norm(tangent)
            print(f"{label}: strain perturbation={step:.0e}, relative tangent difference={error:.3e}")
            if error > 1e-5:
                raise RuntimeError("Material tangent check failed")
    print("Both evaluations in each perturbation pair use the same previous history.")


if __name__ == "__main__":
    main()
