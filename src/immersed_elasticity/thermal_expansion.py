"""Isotropic thermal expansion at fixed temperature during a material update.

The expansion coefficient is constant here and expressed per kelvin. Temperatures
may be in degrees Celsius because only T - T_reference enters thermal strain.
Material properties are supplied by the caller and remain fixed in this step.
"""

from dataclasses import replace

import numpy as np

from .plane_stress import update_plane_stress


def update_thermal_plane_stress(total_strain, previous, material, dt, *, temperature,
                                reference_temperature=20.0, expansion_coefficient=12e-6):
    """Return plane-stress response to total strain including thermal expansion.

    Strain order is [exx, eyy, engineering gxy]. The existing 3D material update
    receives total-minus-thermal strain. The reported strain_zz is restored to
    total thickness strain. The returned tangent differentiates stress with
    respect to total in-plane strain while temperature and previous history
    remain fixed. Plastic history excludes thermal strain.
    """
    strain = np.asarray(total_strain, dtype=float)
    if strain.shape != (3,) or not np.all(np.isfinite(strain)):
        raise ValueError("Require finite total strain [epsilon_xx, epsilon_yy, gamma_xy].")
    if not np.all(np.isfinite([temperature, reference_temperature, expansion_coefficient])):
        raise ValueError("Require finite temperatures and thermal expansion coefficient.")
    thermal_strain = expansion_coefficient * (temperature - reference_temperature)
    mechanical_strain = strain - np.array([thermal_strain, thermal_strain, 0.0])
    response = update_plane_stress(mechanical_strain, previous, material, dt)
    return replace(response, strain_zz=response.strain_zz + thermal_strain)
