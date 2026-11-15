"""Inspect the paper's elastic fits and isolate their material-point effect."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.temperature_properties import with_temperature_dependent_elasticity
from immersed_elasticity.viscoplasticity import J2Material, MaterialState, update_material


def main():
    reference = J2Material(relaxation_time=100.)
    print("Temperature-dependent elastic properties from Appendix D (Table 2):")
    print(" T [C]     K [MPa]     G [MPa]     E [MPa]       nu    end-loading tau [MPa]")
    for temperature in (20., 400., 700.):
        material = with_temperature_dependent_elasticity(reference, temperature)
        state = MaterialState()
        # Identical mechanical shear strain history at each fixed temperature.
        # Normal strains here are mechanical strains relative to free expansion.
        for gamma in np.linspace(0, .01, 101)[1:]:
            strain = np.zeros((3, 3))
            strain[0, 1] = strain[1, 0] = gamma / 2
            stress, state = update_material(strain, state, material, .01)
        print(f"{temperature:6.1f} {material.bulk_modulus:11.3f} {material.shear_modulus:11.3f} "
              f"{material.young:11.3f} {material.poisson:8.5f} {stress[0, 1]:24.3f}")
    print("Yield, hardening, and viscosity parameters are held fixed to isolate elastic softening.")
    print("These are material-point comparisons; thermal FEM coupling comes later.")


if __name__ == "__main__":
    main()
