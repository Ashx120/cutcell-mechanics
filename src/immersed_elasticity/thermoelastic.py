"""One-way, constant-property thermoelasticity on immersed Cartesian Q4 meshes.

Nodal temperatures are prescribed input from the heat solver. Mechanical units
are mm, N, MPa; temperatures are degrees Celsius and expansion is per kelvin.
"""

from dataclasses import dataclass

import numpy as np

from .quadrature import cut_rule
from .solver import assemble, dofs, element_kinematics, solve


@dataclass(frozen=True)
class ThermoelasticResult:
    displacement: np.ndarray
    residual: np.ndarray
    points: np.ndarray
    temperature: np.ndarray
    strain: np.ndarray  # Total strain: xx, yy, engineering xy.
    stress: np.ndarray  # Plane stress: xx, yy, xy, in MPa.
    physical: np.ndarray


def solve_thermoelastic(mesh, material, nodal_temperature, constraints, *,
                       reference_temperature=20., expansion_coefficient=12e-6,
                       depth=3, alpha=1e-6):
    """Solve K u = f_thermal with zero applied mechanical force.

    Thermal load is integral(B.T D epsilon_thermal) dV. Constraints remove rigid
    motion or impose restraint; residual is the support reaction at fixed DOFs.
    This elastic checkpoint has no history variables or temperature feedback.
    The fictitious region supplies alpha-scaled elastic stiffness with no thermal
    eigenstrain. Reported fictitious stresses include alpha; physical stresses
    do not. Heat and mechanics should use the same mesh, depth, and alpha.
    """
    temperature = np.asarray(nodal_temperature, dtype=float)
    if temperature.shape != (len(mesh.nodes),) or not np.all(np.isfinite(temperature)):
        raise ValueError("Require one finite temperature per mesh node.")
    if not np.all(np.isfinite([reference_temperature, expansion_coefficient])):
        raise ValueError("Require finite reference temperature and expansion coefficient.")
    matrix, _ = assemble(mesh, material, depth=depth, alpha=alpha)
    force = np.zeros(2 * len(mesh.nodes))
    integration = []
    for cell in mesh.cells:
        xy = mesh.nodes[cell]
        points, weights = cut_rule((xy[0, 0], xy[1, 0], xy[0, 1], xy[3, 1]),
                                   mesh.radius, depth, alpha)
        shape, strain_matrix = element_kinematics(xy, points)
        point_temperature = shape @ temperature[cell]
        expansion = expansion_coefficient * (point_temperature - reference_temperature)
        thermal_strain = expansion[:, None] * np.array([1., 1., 0.])
        physical = np.sum(points**2, axis=1) >= mesh.radius**2
        thermal_strain[~physical] = 0.
        load = material.thickness * np.einsum(
            "qia,ij,qj,q->a", strain_matrix, material.D, thermal_strain, weights)
        ids = dofs(cell)
        np.add.at(force, ids, load)
        integration.append((ids, points, strain_matrix, point_temperature, thermal_strain, physical))
    displacement, residual = solve(matrix, force, constraints)
    strains, stresses = [], []
    for ids, _, strain_matrix, _, thermal_strain, physical in integration:
        strain = np.einsum("qia,a->qi", strain_matrix, displacement[ids])
        strains.append(strain)
        stresses.append(((strain - thermal_strain) @ material.D.T) * np.where(physical, 1., alpha)[:, None])
    return ThermoelasticResult(
        displacement, residual, np.concatenate([item[1] for item in integration]),
        np.concatenate([item[3] for item in integration]),
        np.concatenate(strains), np.concatenate(stresses),
        np.concatenate([item[5] for item in integration]))
