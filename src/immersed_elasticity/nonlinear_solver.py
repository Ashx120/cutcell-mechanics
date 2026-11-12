"""Incremental plane-stress FEM with fixed quadtree or compressed quadrature."""

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from .plane_stress import update_plane_stress
from .quadrature import classify, cut_rule
from .moment_fitting import compress_rule
from .solver import Material, dofs, element_kinematics
from .viscoplasticity import MaterialState


@dataclass(frozen=True)
class IncrementState:
    displacement: np.ndarray
    history: tuple


@dataclass(frozen=True)
class IncrementResult:
    state: IncrementState
    residual: np.ndarray  # Internal minus external forces; reactions on fixed DOFs.
    stress: np.ndarray  # One [xx, yy, xy] row per integration point.
    iterations: int
    free_residual_norm: float


class NonlinearModel:
    def __init__(self, mesh, material, thickness=1.0, *, depth=3, alpha=1e-6, method="quadtree"):
        if not np.isfinite(thickness) or thickness <= 0:
            raise ValueError("Require positive finite thickness.")
        if not 0 < alpha <= 1 or depth < 0 or int(depth) != depth:
            raise ValueError("Require 0 < alpha <= 1 and nonnegative integer depth.")
        if method not in ("quadtree", "moment-fitting"):
            raise ValueError("Choose 'quadtree' or 'moment-fitting'.")
        self.material = material
        self.elastic_matrix = Material(material.young, material.poisson).D
        self.ndofs = 2 * len(mesh.nodes)
        self.elements = []
        self.fallback_cells = 0
        self.plane_stress_updates = 0
        self.evaluations = 0
        point_locations = []
        for cell in mesh.cells:
            xy = mesh.nodes[cell]
            bounds = (xy[0, 0], xy[1, 0], xy[0, 1], xy[3, 1])
            points, weights = cut_rule(bounds, mesh.radius, depth=depth, alpha=alpha)
            if method == "moment-fitting" and classify(bounds, mesh.radius) == 0:
                physical = np.sum(points**2, axis=1) >= mesh.radius**2
                regions = []
                try:
                    for mask in (physical, ~physical):
                        if np.any(mask):
                            regions.append(compress_rule(points[mask], weights[mask], bounds, degree=2))
                    points = np.concatenate([region[0] for region in regions])
                    weights = np.concatenate([region[1] for region in regions])
                except RuntimeError:
                    self.fallback_cells += 1
            point_locations.append(points)
            _, matrices = element_kinematics(xy, points)
            self.elements.append((dofs(cell), matrices, weights * thickness))
        self.point_count = sum(len(weights) for _, _, weights in self.elements)
        self.points = np.concatenate(point_locations)
        self.physical = np.sum(self.points**2, axis=1) >= mesh.radius**2
        if self.fallback_cells:
            warnings.warn(f"Used quadtree fallback in {self.fallback_cells} cells.", RuntimeWarning, stacklevel=2)

    def initial_state(self):
        return IncrementState(np.zeros(self.ndofs), tuple(MaterialState() for _ in range(self.point_count)))

    def evaluate(self, displacement, committed_history, dt):
        self.evaluations += 1
        internal = np.zeros(self.ndofs)
        rows, columns, values = [], [], []
        trial_history, stresses = [], []
        point = 0
        for ids, matrices, weights in self.elements:
            stiffness = np.zeros((8, 8))
            element_force = np.zeros(8)
            for matrix, weight in zip(matrices, weights):
                strain = matrix @ displacement[ids]
                if self.physical[point]:
                    self.plane_stress_updates += 1
                    response = update_plane_stress(strain, committed_history[point], self.material, dt)
                    stress, tangent, trial_state = response.stress, response.tangent, response.state
                else:
                    # The void supplies alpha-scaled linear elastic stiffness.
                    # Its integration weights already include alpha.
                    stress = self.elastic_matrix @ strain
                    tangent = self.elastic_matrix
                    trial_state = committed_history[point]
                stiffness += matrix.T @ tangent @ matrix * weight
                element_force += matrix.T @ stress * weight
                trial_history.append(trial_state)
                stresses.append(stress)
                point += 1
            np.add.at(internal, ids, element_force)
            rows.extend(np.repeat(ids, 8))
            columns.extend(np.tile(ids, 8))
            values.extend(stiffness.ravel())
        tangent = coo_matrix((values, (rows, columns)), shape=(self.ndofs, self.ndofs)).tocsr()
        return internal, tangent, tuple(trial_history), np.asarray(stresses)


def solve_increment(model, previous, constraints, dt, force=None, *, max_iterations=25):
    """Solve one time increment without modifying the previous converged state.

    On failure, raise and leave history uncommitted. The caller advances by
    assigning result.state only after this function returns successfully.
    """
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Require dt > 0.")
    if max_iterations < 0 or int(max_iterations) != max_iterations:
        raise ValueError("Require a nonnegative integer iteration limit.")
    if len(previous.history) != model.point_count or previous.displacement.shape != (model.ndofs,):
        raise ValueError("Previous state must belong to this model and integration rule.")
    external = np.zeros(model.ndofs) if force is None else np.asarray(force, dtype=float)
    if external.shape != (model.ndofs,) or not np.all(np.isfinite(external)):
        raise ValueError("Require a finite force vector with one entry per DOF.")
    fixed = np.array(sorted(constraints), dtype=int)
    prescribed = np.array([constraints[i] for i in fixed], dtype=float)
    if np.any(fixed < 0) or np.any(fixed >= model.ndofs) or not np.all(np.isfinite(prescribed)):
        raise ValueError("Invalid prescribed displacement.")
    free = np.setdiff1d(np.arange(model.ndofs), fixed)
    displacement = previous.displacement.copy()
    displacement[fixed] = prescribed
    for iteration in range(max_iterations + 1):
        internal, tangent, trial_history, stress = model.evaluate(displacement, previous.history, dt)
        residual = internal - external
        residual_norm = float(np.linalg.norm(residual[free]))
        tolerance = 1e-8 + 1e-9 * max(np.linalg.norm(internal), np.linalg.norm(external))
        if residual_norm <= tolerance:
            return IncrementResult(IncrementState(displacement, trial_history), residual,
                                   stress, iteration, residual_norm)
        if iteration == max_iterations:
            break
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            correction = spsolve(tangent[free][:, free], -residual[free])
        if not np.all(np.isfinite(correction)):
            raise RuntimeError("Non-finite Newton correction; check supports and tangent.")
        displacement[free] += correction
    raise RuntimeError(f"Equilibrium failed after {max_iterations} corrections; history was not committed.")
