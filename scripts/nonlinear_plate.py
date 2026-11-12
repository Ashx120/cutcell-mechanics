"""Load and hold a quarter plate using fixed quadtree or moment-fitted rules."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from immersed_elasticity.nonlinear_solver import NonlinearModel, solve_increment
from immersed_elasticity.solver import structured_mesh, symmetry_constraints
from immersed_elasticity.viscoplasticity import J2Material


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, default=6)
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--steps", type=int, default=10, help="Increments in each one-second phase.")
    parser.add_argument("--output", type=Path, default=Path("results/walkthrough/v3-plate"))
    parser.add_argument("--method", choices=("quadtree", "moment-fitting"), default="quadtree")
    parser.add_argument("--alpha", type=float, default=1e-6,
                        help="Fictitious-domain stiffness factor (default: 1e-6).")
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("steps must be positive")

    mesh = structured_mesh(args.mesh, args.mesh, radius=1.)
    material = J2Material(relaxation_time=100)
    start = perf_counter()
    model = NonlinearModel(mesh, material, depth=args.depth, method=args.method, alpha=args.alpha)
    state = model.initial_state()
    setup_seconds = perf_counter() - start
    right = np.flatnonzero(np.isclose(mesh.nodes[:, 0], mesh.length))
    times = np.linspace(0, 2, 2 * args.steps + 1)
    imposed = np.interp(times, [0, 1, 2], [0, .012, .012])
    rows = [[0., 0., 0., 0., 0., 0., 0.]]
    print(f"Method: {args.method}; active cells: {len(mesh.cells)}, DOFs: {model.ndofs}, integration points: {model.point_count}", flush=True)
    print(f"Physical history points: {np.count_nonzero(model.physical)}", flush=True)
    print(f"Fictitious stiffness factor: {args.alpha:g}", flush=True)
    start = perf_counter()
    for step in range(1, len(times)):
        constraints = symmetry_constraints(mesh)
        constraints.update({2 * node: imposed[step] for node in right})
        result = solve_increment(model, state, constraints, times[step] - times[step-1])
        state = result.state
        force = result.residual[2 * right].sum()
        accumulated = np.array([point.accumulated_strain for point in state.history])
        support = result.residual[np.array(sorted(constraints))]
        reactions = np.zeros(model.ndofs)
        reactions[np.array(sorted(constraints))] = support
        imbalance = np.linalg.norm(reactions.reshape(-1, 2).sum(axis=0))
        rows.append([times[step], imposed[step], force, accumulated[model.physical].max(),
                     result.iterations, result.free_residual_norm, imbalance])
        if step % args.steps == 0:
            print(f"t={times[step]:g} s: right force={force:.3f} N, max accumulated strain={rows[-1][3]:.6f}", flush=True)

    solve_seconds = perf_counter() - start
    history = np.asarray(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "history.csv", history, delimiter=",", comments="",
               header="time_s,right_displacement_mm,right_force_N,max_accumulated_strain,Newton_corrections,free_residual_N,force_imbalance_N")
    metrics = {
        "method": args.method, "mesh_n": args.mesh, "depth": args.depth,
        "steps_per_phase": args.steps, "alpha": args.alpha, "thickness_mm": 1.,
        "plate_side_mm": 4., "hole_radius_mm": 1., "peak_displacement_mm": .012,
        "material": asdict(material), "integration_points": model.point_count,
        "physical_history_points": int(np.count_nonzero(model.physical)),
        "allocated_history_slots": len(state.history),
        "fallback_cells": model.fallback_cells,
        "plane_stress_updates": model.plane_stress_updates,
        "global_evaluations": model.evaluations, "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds, "python": platform.python_version(),
        "numpy": np.__version__,
        "history_storage_note": "Fictitious points retain empty state placeholders; counts are not measured memory bytes.",
        "timing_note": "Single run; setup includes rules, B matrices, and initial history. Solve includes all increments, Newton assembly and solves, and console reporting; excludes plotting."
    }
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    np.savetxt(args.output / "plastic_strain.csv", np.column_stack((model.points[model.physical], accumulated[model.physical])),
               delimiter=",", comments="", header="x_mm,y_mm,accumulated_strain")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    axes[0].plot(history[:, 0], history[:, 2], "o-", markersize=3)
    axes[0].axvspan(1, 2, alpha=.15, color="#94a3b8")
    axes[0].set(xlabel="Time [s]", ylabel="Right-edge force [N]", title="Loading and displacement hold")
    axes[1].plot(history[:, 1], history[:, 2], "o-", markersize=3)
    axes[1].set(xlabel="Right-edge displacement [mm]", ylabel="Right-edge force [N]", title="Reaction-displacement response")
    fig.savefig(args.output / "response.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
    points = model.points[model.physical]
    sc = ax.scatter(*points.T, c=accumulated[model.physical], s=8, cmap="viridis", vmin=0)
    ax.add_patch(Circle((0, 0), mesh.radius, facecolor="white", edgecolor="#334155"))
    ax.set(xlim=(0, mesh.length), ylim=(0, mesh.height), aspect="equal", xlabel="x [mm]", ylabel="y [mm]",
           title="Accumulated viscoplastic strain after the hold")
    fig.colorbar(sc, ax=ax, label="Accumulated equivalent viscoplastic strain")
    fig.savefig(args.output / "plastic_strain.png", dpi=160)
    plt.close(fig)
    print(f"Maximum free residual: {history[:, 5].max():.3e} N; force imbalance: {history[:, 6].max():.3e} N")
    print(f"Maximum Newton corrections: {int(history[:, 4].max())}")
    print(f"Plane-stress updates: {model.plane_stress_updates}; global evaluations: {model.evaluations}")
    print(f"Setup: {setup_seconds:.3f} s; solve: {solve_seconds:.3f} s (single-run timings)")
    print(f"Saved history, strain samples, plots, and metrics to {args.output}")


if __name__ == "__main__":
    main()
