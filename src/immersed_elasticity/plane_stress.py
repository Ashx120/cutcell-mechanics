"""Plane-stress reduction of the 3D J2 material with planar loading history."""

from dataclasses import dataclass

import numpy as np

from .viscoplasticity import MaterialState, update_with_tangent


@dataclass(frozen=True)
class PlaneStressResponse:
    stress: np.ndarray  # xx, yy, xy in MPa
    state: MaterialState
    tangent: np.ndarray  # derivative with respect to exx, eyy, engineering gxy
    strain_zz: float
    stress_zz: float
    iterations: int


def update_plane_stress(in_plane_strain, previous, material, dt):
    """Solve sigma_zz=0, then condense the tangent with a Schur complement.

    Transverse shear strains and their plastic history are zero. Every local
    Newton iteration starts from the same previous committed material state.
    Returned state remains a trial state until the caller commits it.
    """
    values = np.asarray(in_plane_strain, dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("Require finite [epsilon_xx, epsilon_yy, gamma_xy].")
    if np.any(np.abs(previous.plastic_strain[[0, 1], 2]) > 1e-12):
        raise ValueError("Plane stress requires zero transverse plastic shear history.")
    strain = np.zeros((3, 3))
    strain[0, 0], strain[1, 1] = values[:2]
    strain[0, 1] = strain[1, 0] = values[2] / 2
    # An elastic plane-stress estimate supplies the first thickness strain.
    strain[2, 2] = previous.plastic_strain[2, 2] - material.poisson / (1 - material.poisson) * (
        values[0] + values[1] - previous.plastic_strain[0, 0] - previous.plastic_strain[1, 1])

    in_plane = [0, 1, 3]
    for iteration in range(30):
        stress, state, tangent = update_with_tangent(strain, previous, material, dt)
        if abs(stress[2, 2]) <= 1e-8:  # MPa
            reduced = (tangent[np.ix_(in_plane, in_plane)]
                       - np.outer(tangent[in_plane, 2], tangent[2, in_plane]) / tangent[2, 2])
            return PlaneStressResponse(
                np.array([stress[0, 0], stress[1, 1], stress[0, 1]]), state, reduced,
                float(strain[2, 2]), float(stress[2, 2]), iteration)
        strain[2, 2] -= stress[2, 2] / tangent[2, 2]
    raise RuntimeError("Plane-stress material update failed to converge in 30 iterations.")
