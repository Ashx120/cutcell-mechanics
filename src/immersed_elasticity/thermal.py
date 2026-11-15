"""Stationary conduction with constant isotropic conductivity on immersed Q4 meshes.

Units are mm, W, and degrees Celsius (temperature differences in K).
Conductivity is W/(mm K), thickness is mm, and heat flux is W/mm^2.
This first thermal checkpoint has prescribed nodal temperatures, zero volumetric
source, and zero flux on unconstrained boundaries.
"""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix

from .quadrature import cut_rule
from .solver import element_kinematics, solve


@dataclass(frozen=True)
class ThermalResult:
    temperature: np.ndarray
    heat_residual: np.ndarray  # W; reservoir heat into model at constrained nodes.
    points: np.ndarray
    flux: np.ndarray  # Physical Fourier flux at integration points, W/mm^2.
    physical: np.ndarray  # False for the weakly conducting fictitious region.


def solve_stationary_heat(mesh, conductivity, thickness, prescribed_temperatures, *, depth=3, alpha=1e-6):
    """Solve integral(k * grad(N)^T grad(T) * thickness) dA = 0 at free nodes.

    At constrained nodes, the residual reports signed reservoir heat input.
    Hot reservoirs supply positive heat and cold reservoirs remove heat.
    Temperature has one DOF per node. The same Q4 geometry and displacement
    elimination routine as the mechanical solver are reused.
    Cut cells use quadtree integration. The fictitious region has conductivity
    alpha*k, approximating an insulated hole as alpha tends to zero. Flux there
    is reported with that reduced conductivity; use physical to select material.
    """
    if not np.isfinite(alpha) or not 0 < alpha <= 1 or depth < 0 or int(depth) != depth:
        raise ValueError("Require 0 < alpha <= 1 and nonnegative integer depth.")
    if (not np.isfinite(conductivity) or conductivity <= 0
            or not np.isfinite(thickness) or thickness <= 0):
        raise ValueError("Require positive finite conductivity and thickness.")
    if not prescribed_temperatures:
        raise ValueError("At least one prescribed temperature is required.")
    for node, value in prescribed_temperatures.items():
        if (not isinstance(node, (int, np.integer)) or not 0 <= node < len(mesh.nodes)
                or not np.isfinite(value)):
            raise ValueError("Prescribed temperatures must have valid node indices and finite values.")

    rows, columns, values = [], [], []
    integration = []
    for cell in mesh.cells:
        xy = mesh.nodes[cell]
        points, weights = cut_rule((xy[0, 0], xy[1, 0], xy[0, 1], xy[3, 1]),
                                   mesh.radius, depth, alpha)
        _, strain_matrices = element_kinematics(xy, points)
        gradients = np.stack((strain_matrices[:, 0, 0::2], strain_matrices[:, 1, 1::2]), axis=1)
        element_matrix = conductivity * thickness * np.einsum("qia,qib,q->ab", gradients, gradients, weights)
        rows.extend(np.repeat(cell, 4))
        columns.extend(np.tile(cell, 4))
        values.extend(element_matrix.ravel())
        integration.append((cell, points, gradients))

    size = len(mesh.nodes)
    matrix = coo_matrix((values, (rows, columns)), shape=(size, size)).tocsr()
    temperature, residual = solve(matrix, np.zeros(size), prescribed_temperatures)
    flux = [-conductivity * np.einsum("qia,a->qi", gradient, temperature[cell])
            for cell, _, gradient in integration]
    points = np.concatenate([points for _, points, _ in integration])
    physical = np.sum(points**2, axis=1) >= mesh.radius**2
    flux = np.concatenate(flux) * np.where(physical, 1., alpha)[:, None]
    return ThermalResult(temperature, residual, points, flux, physical)
