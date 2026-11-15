"""Verify free and in-plane-restrained expansion against elastic formulas."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from immersed_elasticity.thermal_expansion import update_thermal_plane_stress
from immersed_elasticity.viscoplasticity import J2Material, MaterialState


# Illustrative constants; these cases remain below the 250 MPa yield threshold.
MATERIAL = J2Material(relaxation_time=100)
REFERENCE_TEMPERATURE = 20.0
TEMPERATURE = 70.0
EXPANSION_COEFFICIENT = 12e-6  # 1/K, distinct from fictitious stiffness alpha.


def main():
    thermal_strain = EXPANSION_COEFFICIENT * (TEMPERATURE - REFERENCE_TEMPERATURE)
    options = dict(temperature=TEMPERATURE, reference_temperature=REFERENCE_TEMPERATURE,
                   expansion_coefficient=EXPANSION_COEFFICIENT)
    free = update_thermal_plane_stress([thermal_strain, thermal_strain, 0],
                                      MaterialState(), MATERIAL, 1.0, **options)
    restrained = update_thermal_plane_stress([0., 0., 0.], MaterialState(), MATERIAL, 1.0, **options)
    expected_stress = -MATERIAL.young * thermal_strain / (1 - MATERIAL.poisson)
    np.testing.assert_allclose(free.stress, 0., atol=1e-8)
    np.testing.assert_allclose(free.strain_zz, thermal_strain, atol=1e-12)
    np.testing.assert_allclose(restrained.stress, [expected_stress, expected_stress, 0.], atol=1e-8)
    assert free.state.accumulated_strain == restrained.state.accumulated_strain == 0
    print(f"Temperature rise: {TEMPERATURE - REFERENCE_TEMPERATURE:g} K")
    print(f"Free normal strain in each direction: {thermal_strain:.6e}")
    print(f"Free expansion stress [xx, yy, xy]: {free.stress} MPa")
    print(f"Restrained in-plane stress [xx, yy, xy]: {restrained.stress} MPa")
    print(f"Exact restrained sigma_xx = sigma_yy: {expected_stress:.6f} MPa")
    print(f"Total thickness strain when restrained in-plane: {restrained.strain_zz:.6e}")
    print(f"Restrained sigma_zz: {restrained.stress_zz:.3e} MPa")
    print("Both cases remain elastic; accumulated viscoplastic strain is zero.")


if __name__ == "__main__":
    main()
