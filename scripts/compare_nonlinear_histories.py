"""Compare matched quadtree and moment-fitted plate reaction histories."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("results/walkthrough/v3-quadtree"))
    parser.add_argument("--fitted", type=Path, default=Path("results/walkthrough/v3-fitted"))
    parser.add_argument("--output", type=Path, default=Path("results/walkthrough/v3-comparison"))
    args = parser.parse_args()
    reference_metrics = json.loads((args.reference / "metrics.json").read_text())
    fitted_metrics = json.loads((args.fitted / "metrics.json").read_text())
    if reference_metrics["method"] != "quadtree" or fitted_metrics["method"] != "moment-fitting":
        parser.error("Expected quadtree reference and moment-fitting comparison.")
    for key in ("mesh_n", "depth", "steps_per_phase", "alpha", "thickness_mm",
                "plate_side_mm", "hole_radius_mm", "peak_displacement_mm", "material"):
        if reference_metrics[key] != fitted_metrics[key]:
            parser.error(f"Run settings differ for {key}.")

    reference = np.genfromtxt(args.reference / "history.csv", delimiter=",", names=True)
    fitted = np.genfromtxt(args.fitted / "history.csv", delimiter=",", names=True)
    if reference.shape != fitted.shape or reference.size < 2:
        parser.error("Require matching histories with at least two samples.")
    for field in ("time_s", "right_displacement_mm"):
        if not np.allclose(reference[field], fitted[field], rtol=0, atol=1e-12):
            parser.error(f"Histories differ in {field}.")
    time = reference["time_s"]
    original_force = reference["right_force_N"]
    fitted_force = fitted["right_force_N"]
    if not np.all(np.isfinite(original_force)) or not np.all(np.isfinite(fitted_force)):
        parser.error("Require finite reaction forces.")
    peak = float(np.max(np.abs(original_force)))
    if peak == 0:
        parser.error("Reference reactions are all zero; relative comparison is undefined.")
    difference = fitted_force - original_force
    index = int(np.argmax(np.abs(difference)))
    summary = {
        "max_absolute_force_difference_N": float(abs(difference[index])),
        "time_of_max_difference_s": float(time[index]),
        "reference_force_at_max_difference_N": float(original_force[index]),
        "fitted_force_at_max_difference_N": float(fitted_force[index]),
        "max_difference_percent_of_reference_peak": float(100 * abs(difference[index]) / peak),
        "relative_sample_vector_norm": float(np.linalg.norm(difference) / np.linalg.norm(original_force)),
        "norm_note": "Euclidean norm over matching time samples; not a continuous-time integral.",
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    np.savetxt(args.output / "force_difference.csv", np.column_stack((time, difference)),
               delimiter=",", comments="", header="time_s,fitted_minus_quadtree_force_N")
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, layout="constrained")
    axes[0].plot(time, original_force, label="Quadtree")
    axes[0].plot(time, fitted_force, "--", label="Moment fitting")
    axes[0].set(ylabel="Right-edge force [N]", title="Matched nonlinear plate histories")
    axes[0].legend()
    axes[1].plot(time, difference)
    axes[1].axhline(0, color="gray", linewidth=.8)
    axes[1].set(xlabel="Time [s]", ylabel="Fitted minus quadtree [N]")
    for axis in axes:
        axis.grid(alpha=.25)
    fig.savefig(args.output / "reaction_comparison.png", dpi=160)
    plt.close(fig)
    print(f"Largest force difference: {abs(difference[index]):.6f} N at t={time[index]:.3f} s")
    print(f"Reference force there: {original_force[index]:.6f} N; fitted: {fitted_force[index]:.6f} N")
    print(f"Largest difference / reference peak: {summary['max_difference_percent_of_reference_peak']:.6f}%")
    print(f"Relative history sample-vector norm: {summary['relative_sample_vector_norm']:.3e}")
    print(f"Saved comparison plot, differences, and summary to {args.output}")


if __name__ == "__main__":
    main()
