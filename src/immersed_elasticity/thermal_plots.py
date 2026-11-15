"""Plots for thermal examples; points are physical integration samples, not contours."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save_thermal_plots(output, history, points, temperature, von_mises, plastic, mesh, loading):
    """Write loading/history plots and unsmoothed final physical-point fields."""
    heating_time = loading["heating_time"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    for ax, column, title, label in zip(
            axes, (1, 3, 4), ("Prescribed temperature", "Stress response", "Plastic strain history"),
            ("Right-edge temperature [deg C]", "Sampled peak von Mises stress [MPa]", "Maximum accumulated strain")):
        ax.plot(history[:, 0], history[:, column], linewidth=2)
        if loading["hold_time"]:
            ax.axvspan(heating_time, history[-1, 0], alpha=.12, color="gray")
        ax.set(xlabel="Time [s]", ylabel=label, title=title)
        ax.grid(alpha=.2)
    fig.savefig(output / "history.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), layout="constrained")
    for ax, values, title, label, cmap in zip(
            axes, (temperature, von_mises, plastic),
            ("Final temperature", "Final von Mises stress", "Final accumulated plastic strain"),
            ("deg C", "MPa", "Dimensionless strain"), ("inferno", "magma", "viridis")):
        field = ax.scatter(*points.T, c=values, s=10, cmap=cmap)
        _outline(ax, mesh)
        ax.set(title=title)
        fig.colorbar(field, ax=ax, label=label, shrink=.8)
    fig.suptitle("Physical integration-point samples | colour scales differ between fields", fontsize=11)
    fig.savefig(output / "fields.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
    field = ax.scatter(*points.T, c=plastic, s=10, cmap="viridis")
    _outline(ax, mesh)
    ax.set(title="Accumulated viscoplastic strain after the hold")
    fig.colorbar(field, ax=ax, label="Dimensionless strain")
    fig.savefig(output / "plastic_strain.png", dpi=160)
    plt.close(fig)


def _outline(ax, mesh):
    angle = np.linspace(0, np.pi/2, 100)
    ax.plot(mesh.radius*np.cos(angle), mesh.radius*np.sin(angle), color="#333333")
    ax.set(xlim=(0, mesh.length), ylim=(0, mesh.height), aspect="equal", xlabel="x [mm]", ylabel="y [mm]")
