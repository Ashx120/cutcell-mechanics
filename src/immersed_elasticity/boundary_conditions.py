"""The two supported mechanical boundary layouts for the thermal examples."""

import numpy as np

from .solver import symmetry_constraints


def thermal_constraints(mesh, supports="symmetry"):
    """Return edge rollers or three constraints removing rigid motion only."""
    if supports == "symmetry":
        return symmetry_constraints(mesh)
    if supports != "minimal":
        raise ValueError("Supports must be symmetry or minimal.")
    x, y = mesh.nodes.T
    first = int(np.flatnonzero(np.isclose(x, 0.) & np.isclose(y, mesh.height))[0])
    second = int(np.flatnonzero(np.isclose(x, mesh.length) & np.isclose(y, mesh.height))[0])
    return {2*first: 0., 2*first+1: 0., 2*second+1: 0.}
