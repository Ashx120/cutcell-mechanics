"""Inspect thermal hardening softening independently of initial yield."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from immersed_elasticity.temperature_properties import with_temperature_dependent_hardening
from immersed_elasticity.viscoplasticity import J2Material


def main():
    reference = J2Material(initial_yield=302.3605)
    print("Hardening from equations (18)-(20) and Appendix D, Table 2:")
    print(" T [C]    Q(T) [MPa]    yield at p=0    yield at p=0.02    saturation [MPa]")
    for temperature in (20., 400., 700.):
        material = with_temperature_dependent_hardening(reference, temperature)
        print(f"{temperature:6.1f} {material.saturation_increment:13.6f} "
              f"{material.yield_stress(0):15.6f} {material.yield_stress(.02):18.6f} "
              f"{material.initial_yield + material.saturation_increment:19.6f}")
    print("Initial yield is deliberately fixed at 302.3605 MPa; viscosity is unchanged.")
    print("These are equivalent yield stresses, not shear stresses or viscous overstresses.")


if __name__ == "__main__":
    main()
