"""Check thickness stress and the condensed material tangent."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.plane_stress import update_plane_stress
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


def main():
    material = J2Material(relaxation_time=100)
    previous = MaterialState()
    for label, strain in (("Elastic", np.array([1e-4, 0., 0.])),
                          ("Viscoplastic", np.array([.004, -.0005, .002]))):
        response = update_plane_stress(strain, previous, material, .01)
        numerical = np.empty((3, 3))
        for column in range(3):
            perturbation = np.eye(3)[column] * 1e-7
            plus = update_plane_stress(strain + perturbation, previous, material, .01)
            minus = update_plane_stress(strain - perturbation, previous, material, .01)
            numerical[:, column] = (plus.stress - minus.stress) / 2e-7
        error = np.linalg.norm(response.tangent - numerical) / np.linalg.norm(response.tangent)
        print(f"{label}: epsilon_zz={response.strain_zz:.6e}, sigma_zz={response.stress_zz:.3e} MPa")
        print(f"  Local Newton corrections: {response.iterations}; relative tangent difference: {error:.3e}")
        if abs(response.stress_zz) > 1e-8 or error > 1e-6:
            raise RuntimeError("Plane-stress verification failed")


if __name__ == "__main__":
    main()
