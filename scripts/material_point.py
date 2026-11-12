"""Load, hold, and unload one J2 viscoplastic point in simple shear."""

from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from immersed_elasticity.viscoplasticity import J2Material, MaterialState, equivalent_stress, update_material

# Illustrative parameters, not calibrated to a particular alloy.
MATERIAL = J2Material(relaxation_time=100.0)
TIME_STEP = 0.01  # seconds
PEAK_SHEAR_STRAIN = 0.01  # engineering gamma_xy
OUTPUT = Path(__file__).resolve().parents[1] / "results" / "walkthrough" / "v3-material-point"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dt", type=float, default=TIME_STEP,
                        help="Maximum time step in seconds (default: 0.01).")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--loading-time", type=float, default=1.0,
                        help="Time to reach peak strain; hold is 2 s and unloading is 1 s.")
    args = parser.parse_args()
    if not np.isfinite(args.dt) or args.dt <= 0:
        parser.error("dt must be positive and finite")
    if not np.isfinite(args.loading_time) or args.loading_time <= 0:
        parser.error("loading-time must be positive and finite")
    load_end = args.loading_time
    hold_end = load_end + 2.0
    unload_end = hold_end + 1.0
    phase_times = [0.0, load_end, hold_end, unload_end]
    # End each phase exactly at its prescribed change of loading direction.
    times = np.concatenate([
        np.linspace(start, end, int(np.ceil((end-start)/args.dt)) + 1)[:-1]
        for start, end in zip(phase_times[:-1], phase_times[1:])
    ] + [np.array([unload_end])])
    shear_strain = np.interp(times, phase_times, [0, PEAK_SHEAR_STRAIN, PEAK_SHEAR_STRAIN, 0])
    state = MaterialState()
    rows = [[0., 0., 0., 0., MATERIAL.initial_yield, 0.]]
    for index in range(1, len(times)):
        strain = np.zeros((3, 3))
        strain[0, 1] = strain[1, 0] = shear_strain[index] / 2
        stress, state = update_material(strain, state, MATERIAL, times[index] - times[index - 1])
        rows.append([times[index], shear_strain[index], stress[0, 1],
                     equivalent_stress(stress), MATERIAL.yield_stress(state.accumulated_strain),
                     state.accumulated_strain])
    history = np.asarray(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "history.csv", history, delimiter=",", comments="",
               header="time_s,gamma_xy,tau_xy_MPa,equivalent_stress_MPa,yield_stress_MPa,accumulated_strain")
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True, layout="constrained")
    axes[0].plot(times, shear_strain, color="#2563eb")
    axes[0].set_ylabel("Engineering shear strain")
    axes[1].plot(times, history[:, 2], label="Shear stress")
    axes[1].set_ylabel("Shear stress [MPa]")
    axes[2].plot(times, history[:, 5], color="#16a34a")
    axes[2].set(ylabel="Accumulated viscoplastic strain", xlabel="Time [s]")
    for ax in axes:
        ax.axvspan(load_end, hold_end, color="#94a3b8", alpha=.15)
        ax.grid(alpha=.25)
    fig.suptitle("One material point | simple shear | shaded interval: strain hold")
    fig.savefig(args.output / "history.png", dpi=160)
    plt.close(fig)
    print(f"Maximum step: {args.dt:g} s | Increments: {len(times)-1}")
    print(f"Loading shear-strain rate: {PEAK_SHEAR_STRAIN / load_end:g} 1/s")
    for time in (load_end, hold_end, unload_end):
        row = history[np.argmin(abs(times-time))]
        print(f"t={time} s: gamma={row[1]:.4f}, shear stress={row[2]:.3f} MPa, accumulated strain={row[5]:.6f}")
    print(f"Saved history.csv and history.png to {args.output}")


if __name__ == "__main__":
    main()
