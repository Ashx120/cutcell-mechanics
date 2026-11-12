"""Isothermal, small-strain 3D J2 Perzyna material-point update.

Tensor shear components are used, so epsilon_xy = engineering gamma_xy / 2.
Stress is in MPa, time in seconds, and strain is dimensionless. This implements
the fixed-temperature form of the paper's equations (16), (18)-(22), using
sigma_y = sigma_y0 + Q * (1 - exp(-b*p)).
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class J2Material:
    young: float = 210_000.0
    poisson: float = 0.3
    initial_yield: float = 250.0
    saturation_increment: float = 100.0
    hardening_rate: float = 20.0
    relaxation_time: float = 1.0
    rate_exponent: float = 1.0

    def __post_init__(self):
        values = list(vars(self).values())
        if (not np.all(np.isfinite(values)) or self.young <= 0
                or not -1 < self.poisson < 0.5 or self.initial_yield <= 0
                or self.saturation_increment < 0 or self.hardening_rate < 0
                or self.relaxation_time <= 0 or self.rate_exponent <= 0):
            raise ValueError("Require finite elastic parameters, positive yield/time/exponent, and nonnegative hardening.")

    @property
    def shear_modulus(self):
        return self.young / (2 * (1 + self.poisson))

    @property
    def bulk_modulus(self):
        return self.young / (3 * (1 - 2 * self.poisson))

    def yield_stress(self, accumulated_strain):
        return self.initial_yield + self.saturation_increment * (-np.expm1(-self.hardening_rate * accumulated_strain))


@dataclass(frozen=True)
class MaterialState:
    plastic_strain: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    accumulated_strain: float = 0.0

    def __post_init__(self):
        strain = np.array(self.plastic_strain, dtype=float, copy=True)
        if (strain.shape != (3, 3) or not np.all(np.isfinite(strain))
                or not np.allclose(strain, strain.T, rtol=0, atol=1e-12)
                or abs(np.trace(strain)) > 1e-12
                or not np.isfinite(self.accumulated_strain) or self.accumulated_strain < 0):
            raise ValueError("Require symmetric traceless plastic strain and finite nonnegative accumulated strain.")
        strain.setflags(write=False)
        object.__setattr__(self, "plastic_strain", strain)


def equivalent_stress(stress):
    deviator = stress - np.trace(stress) * np.eye(3) / 3
    return float(np.sqrt(1.5 * np.sum(deviator * deviator)))


def update_material(total_strain, previous, material, dt):
    """Return stress and a new trial state without changing the previous state.

    Backward Euler and associative J2 flow reduce to a scalar equation for
    delta_p. A bracketed solve stays between zero and the full deviatoric return.
    The caller decides when to commit the returned state. Use update_with_tangent
    to additionally obtain its analytical derivative.
    """
    strain = np.asarray(total_strain, dtype=float)
    if (strain.shape != (3, 3) or not np.all(np.isfinite(strain))
            or not np.allclose(strain, strain.T, rtol=0, atol=1e-12)
            or not np.isfinite(dt) or dt <= 0):
        raise ValueError("Require a finite symmetric 3x3 strain tensor and dt > 0.")

    elastic_trial = strain - previous.plastic_strain
    volumetric_strain = np.trace(elastic_trial)
    trial_deviator = 2 * material.shear_modulus * (elastic_trial - volumetric_strain * np.eye(3) / 3)
    pressure_part = material.bulk_modulus * volumetric_strain * np.eye(3)
    trial_equivalent = float(np.sqrt(1.5 * np.sum(trial_deviator**2)))
    if trial_equivalent <= material.yield_stress(previous.accumulated_strain):
        return pressure_part + trial_deviator, previous

    def residual(delta_p):
        q = trial_equivalent - 3 * material.shear_modulus * delta_p
        threshold = material.yield_stress(previous.accumulated_strain + delta_p)
        overstress = max(q / threshold - 1, 0.0)
        with np.errstate(over="raise", invalid="raise"):
            increment = dt / material.relaxation_time * np.power(overstress, 1 / material.rate_exponent)
        return delta_p - increment

    delta_p = brentq(residual, 0, trial_equivalent / (3 * material.shear_modulus), xtol=1e-14, rtol=1e-12)
    direction = trial_deviator / trial_equivalent
    plastic_strain = previous.plastic_strain + 1.5 * delta_p * direction
    state = MaterialState(plastic_strain, previous.accumulated_strain + delta_p)
    stress = pressure_part + trial_deviator - 3 * material.shear_modulus * delta_p * direction
    return stress, state


# Stress order: xx, yy, zz, xy, yz, xz. Strain order uses engineering shear.
COMPONENTS = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))


def engineering_strain_basis(column):
    basis = np.zeros((3, 3))
    i, j = COMPONENTS[column]
    if i == j:
        basis[i, j] = 1.0
    else:
        basis[i, j] = basis[j, i] = 0.5
    return basis


def stress_vector(stress):
    return np.array([stress[i, j] for i, j in COMPONENTS])


def update_with_tangent(total_strain, previous, material, dt):
    """Return stress, trial state, and the analytical backward-Euler 6x6 tangent.

    Differentiate with previous history and dt held fixed. Shear columns use
    engineering strain, so an elastic shear diagonal entry equals G. The
    derivative is branch-specific at the elastic/viscoplastic transition.
    """
    stress, state = update_material(total_strain, previous, material, dt)
    elastic_trial = np.asarray(total_strain) - previous.plastic_strain
    trial_deviator = 2 * material.shear_modulus * (elastic_trial - np.trace(elastic_trial) * np.eye(3) / 3)
    q_trial = float(np.sqrt(1.5 * np.sum(trial_deviator**2)))
    delta_p = state.accumulated_strain - previous.accumulated_strain
    active = delta_p > 0
    if active:
        threshold = material.yield_stress(state.accumulated_strain)
        hardening = (material.saturation_increment * material.hardening_rate
                     * np.exp(-material.hardening_rate * state.accumulated_strain))
        viscous_factor = (material.relaxation_time * delta_p / dt)**material.rate_exponent
        # Differentiate q_trial = 3G*delta_p + sigma_y(p)*(1+viscous_factor).
        increment_derivative = 1 / (3 * material.shear_modulus + hardening * (1 + viscous_factor)
                                   + threshold * material.rate_exponent * viscous_factor / delta_p)
        scale = 1 - 3 * material.shear_modulus * delta_p / q_trial
        scale_derivative = -3 * material.shear_modulus * (increment_derivative / q_trial - delta_p / q_trial**2)

    tangent = np.empty((6, 6))
    for column in range(6):
        strain_direction = engineering_strain_basis(column)
        deviator_direction = 2 * material.shear_modulus * (strain_direction - np.trace(strain_direction) * np.eye(3) / 3)
        stress_direction = material.bulk_modulus * np.trace(strain_direction) * np.eye(3)
        if active:
            q_direction = 1.5 * np.sum(trial_deviator * deviator_direction) / q_trial
            stress_direction += scale * deviator_direction + scale_derivative * q_direction * trial_deviator
        else:
            stress_direction += deviator_direction
        tangent[:, column] = stress_vector(stress_direction)
    return stress, state, tangent
