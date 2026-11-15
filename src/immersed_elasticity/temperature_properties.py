"""Temperature-dependent elastic fits from the source paper's Appendix D.

Temperatures are in degrees Celsius and stress/moduli are in MPa. The absolute
reference 293.15 K becomes 20 deg C; temperature differences have the same value
in kelvin. Evaluation is limited to 0-1500 deg C, the plotted range in Appendix C.
These fits describe the paper's steel parameters, not a new material calibration.
"""

from dataclasses import replace

import numpy as np


def with_temperature_dependent_elasticity(reference_material, temperature):
    """Replace only elastic parameters using the paper's bulk/shear fits.

    Initial yield, hardening, and viscosity remain those supplied in the reference
    material. This isolates the elastic-property checkpoint. The returned material
    is evaluated at a fixed temperature; no history is reset or transferred here.
    """
    if not np.isfinite(temperature) or not 0 <= temperature <= 1500:
        raise ValueError("This property fit is restricted to finite temperatures in 0-1500 deg C.")
    delta_temperature = temperature - 20.0

    def fitted_modulus(reference, omega, chi, phi):
        factor = 1 - 2 * omega / (np.pi - 2 * phi) * (
            np.arctan(chi * delta_temperature + np.tan(phi)) - phi)
        return reference * factor

    # Table 2 lists K0 and G0 in GPa; convert to MPa here.
    bulk = fitted_modulus(163_195.2, .8592410, .004509029, -1.338357)
    shear = fitted_modulus(80_891.53, 1.265276, .002136860, -1.071127)
    young = 9 * bulk * shear / (3 * bulk + shear)
    poisson = (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear))
    return replace(reference_material, young=float(young), poisson=float(poisson))


def with_temperature_dependent_hardening(reference_material, temperature):
    """Replace Q and b using equations (18)-(20) and Appendix D, Table 2.

    Q(T) = 206.3425 * Phi_chi(T) MPa and b = 1 / 0.02784959.
    Elasticity, initial yield, and viscosity remain supplied by the caller.
    In particular, this is not a complete temperature-dependent yield model:
    Appendix D does not tabulate the separate initial-yield fit coefficients.
    """
    if not np.isfinite(temperature) or not 0 <= temperature <= 1500:
        raise ValueError("This property fit is restricted to finite temperatures in 0-1500 deg C.")
    omega, chi, phi = -39.22882, -.009042314, -1.461198
    factor = 1 - 2 * omega / (np.pi - 2 * phi) * (
        np.arctan(chi * (temperature - 20.) + np.tan(phi)) - phi)
    return replace(reference_material, saturation_increment=float(206.3425 * factor),
                   hardening_rate=1 / .02784959)
