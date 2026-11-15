"""Compare saved thermal runs at matching nodes and times, without solving again."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from immersed_elasticity.thermal_settings import settings_from_saved


def compare(reference, fitted):
    """Reject mismatched physical inputs before computing common-node differences."""
    settings = []
    for folder in (reference, fitted):
        resolved = settings_from_saved(json.loads((folder / "settings.json").read_text(encoding="utf-8")))
        resolved["quadrature"].pop("method")
        settings.append(resolved)
    if settings[0] != settings[1]:
        raise ValueError("Physical inputs differ; compare matched geometry, material, supports, time, depth, and alpha.")
    a = np.genfromtxt(reference / "history.csv", delimiter=",", names=True)
    b = np.genfromtxt(fitted / "history.csv", delimiter=",", names=True)
    if a.shape != b.shape or not np.allclose(a["time_s"], b["time_s"], rtol=0, atol=1e-12):
        raise ValueError("History times differ.")
    n = np.genfromtxt(reference / "nodes.csv", delimiter=",", names=True)
    m = np.genfromtxt(fitted / "nodes.csv", delimiter=",", names=True)
    if n.shape != m.shape:
        raise ValueError("Nodal fields differ in size.")
    for key in ("x_mm", "y_mm", "temperature_C"):
        if not np.allclose(n[key], m[key], rtol=0, atol=1e-12):
            raise ValueError("Nodal coordinates or final temperatures differ.")
    u = np.column_stack((n["ux_mm"], n["uy_mm"]))
    v = np.column_stack((m["ux_mm"], m["uy_mm"]))
    norm = np.linalg.norm(u)
    absolute = float(np.linalg.norm(v-u))
    summary = dict(max_mean_right_displacement_difference_mm=float(np.max(np.abs(a["mean_right_ux_mm"]-b["mean_right_ux_mm"]))),
                   final_nodal_displacement_difference_norm_mm=absolute,
                   final_relative_nodal_displacement_norm=float(absolute/norm) if norm else (0. if absolute == 0 else None),
                   max_sampled_peak_stress_history_difference_MPa=float(np.max(np.abs(a["max_von_mises_MPa"]-b["max_von_mises_MPa"]))),
                   max_sampled_peak_strain_history_difference=float(np.max(np.abs(a["max_accumulated_strain"]-b["max_accumulated_strain"]))),
                   note="Sampled stress/strain maxima are not matched-point field error norms. Nodal displacement norm includes all active nodes.")
    return summary, a, b


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("results/reference/thermal/quadtree"))
    parser.add_argument("--fitted", type=Path, default=Path("results/reference/thermal/fitted"))
    parser.add_argument("--output", type=Path, default=Path("results/runs/comparison"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary, a, b = compare(args.reference, args.fitted)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    if args.output.exists() and (not args.output.is_dir() or (any(args.output.iterdir()) and not args.overwrite)):
        parser.error("Choose a new output directory or use --overwrite.")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "comparison.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), layout="constrained")
    axes[0].plot(a["time_s"], a["mean_right_ux_mm"], label="Quadtree reference")
    axes[0].plot(b["time_s"], b["mean_right_ux_mm"], "--", label="Moment fitting")
    axes[0].set(ylabel="Mean right-edge ux [mm]", title="Same temperature field and physical inputs")
    axes[0].legend()
    axes[1].plot(a["time_s"], b["mean_right_ux_mm"]-a["mean_right_ux_mm"])
    axes[1].set(xlabel="Time [s]", ylabel="Fitted minus reference [mm]")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(args.output / "comparison.png", dpi=160)
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    print(f"Saved comparison to {args.output.resolve()}")


if __name__ == "__main__":
    main()
