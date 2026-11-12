"""Load and hold a quarter plate using fixed quadtree or moment-fitted rules."""

import argparse
from copy import deepcopy
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
from immersed_elasticity.plate_settings import load_plate_settings, validate_plate_settings


def main(argv=None, *, initial_settings=None, output=None, overwrite=False):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Editable JSON input file; CLI options override its settings.")
    parser.add_argument("--mesh", type=int, help="Override nx and ny (default: 6).")
    parser.add_argument("--depth", type=int, help="Override quadrature depth (default: 2).")
    parser.add_argument("--steps", type=int, help="Override increments per loading/hold phase (default: 10).")
    parser.add_argument("--output", type=Path, default=output if output is not None else Path("results/walkthrough/v3-plate"))
    parser.add_argument("--overwrite", action="store_true", default=overwrite,
                        help="Replace an existing run in the output folder.")
    parser.add_argument("--method", choices=("quadtree", "moment-fitting"))
    parser.add_argument("--alpha", type=float, help="Override fictitious stiffness factor (default: 1e-6).")
    args = parser.parse_args(argv)
    try:
        settings = (deepcopy(initial_settings) if initial_settings is not None and args.config is None
                    else load_plate_settings(args.config))
        if args.mesh is not None:
            settings["geometry"].update(nx=args.mesh, ny=args.mesh)
        if args.steps is not None:
            settings["loading"]["steps_per_phase"] = args.steps
        for key in ("method", "depth", "alpha"):
            if getattr(args, key) is not None:
                settings["quadrature"][key] = getattr(args, key)
        validate_plate_settings(settings)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    # Stop before the costly solve if this run would overwrite previous results.
    if args.output.exists():
        if not args.output.is_dir():
            parser.error(f"Output path is a file: {args.output}. Choose a directory.")
        if any(args.output.iterdir()) and not args.overwrite:
            parser.error(f"Output folder is nonempty: {args.output}. Choose a new folder or use --overwrite.")
    geometry, loading, quadrature = (settings[key] for key in ("geometry", "loading", "quadrature"))
    args.steps = loading["steps_per_phase"]
    args.method, args.depth, args.alpha = (quadrature[key] for key in ("method", "depth", "alpha"))
    mesh = structured_mesh(**geometry)
    material = J2Material(**settings["material"])
    start = perf_counter()
    model = NonlinearModel(mesh, material, thickness=settings["thickness"], **quadrature)
    state = model.initial_state()
    setup_seconds = perf_counter() - start
    right = np.flatnonzero(np.isclose(mesh.nodes[:, 0], mesh.length))
    load_end = loading["loading_time"]
    end_time = load_end + loading["hold_time"]
    times = np.concatenate((np.linspace(0, load_end, args.steps + 1),
                            np.linspace(load_end, end_time, args.steps + 1)[1:]))
    imposed = np.interp(times, [0, load_end, end_time],
                        [0, loading["peak_displacement"], loading["peak_displacement"]])
    rows = [[0., 0., 0., 0., 0., 0., 0.]]
    print(f"Method: {args.method}; active cells: {len(mesh.cells)}, DOFs: {model.ndofs}, integration points: {model.point_count}", flush=True)
    print(f"Physical history points: {np.count_nonzero(model.physical)}", flush=True)
    print(f"Fictitious stiffness factor: {args.alpha:g}", flush=True)
    print(f"Loading ends at {load_end:g} s; hold ends at {end_time:g} s", flush=True)
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
        if step % max(1, args.steps // 4) == 0 or step == len(times)-1:
            print(f"Step {step}/{len(times)-1}, t={times[step]:g} s: right force={force:.3f} N, max accumulated strain={rows[-1][3]:.6f}", flush=True)

    solve_seconds = perf_counter() - start
    history = np.asarray(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    np.savetxt(args.output / "history.csv", history, delimiter=",", comments="",
               header="time_s,right_displacement_mm,right_force_N,max_accumulated_strain,Newton_corrections,free_residual_N,force_imbalance_N")
    metrics = {
        "method": args.method, "mesh_n": geometry["nx"], "depth": args.depth,
        "steps_per_phase": args.steps, "alpha": args.alpha, "thickness_mm": settings["thickness"],
        "plate_side_mm": geometry["length"], "hole_radius_mm": geometry["radius"],
        "peak_displacement_mm": loading["peak_displacement"], "settings": settings,
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
    axes[0].axvspan(load_end, end_time, alpha=.15, color="#94a3b8")
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
