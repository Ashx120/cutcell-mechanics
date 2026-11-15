"""Load and validate the editable settings for the nonlinear plate example."""

from copy import deepcopy
import json
from pathlib import Path

import numpy as np

from .viscoplasticity import J2Material


DEFAULT_SETTINGS = {
    "geometry": {"nx": 6, "ny": 6, "length": 4.0, "height": 4.0, "radius": 1.0},
    "material": {"young": 210000.0, "poisson": 0.3, "initial_yield": 250.0,
                 "saturation_increment": 100.0, "hardening_rate": 20.0,
                 "relaxation_time": 100.0, "rate_exponent": 1.0},
    "thickness": 1.0,
    "loading": {"peak_displacement": 0.012, "loading_time": 1.0,
                "hold_time": 1.0, "steps_per_phase": 10},
    "quadrature": {"method": "quadtree", "depth": 2, "alpha": 1e-6},
}


def load_plate_settings(path=None, *, supplied=None):
    """Merge a partial JSON configuration into defaults and reject unknown keys."""
    settings = deepcopy(DEFAULT_SETTINGS)
    if path is not None and supplied is not None:
        raise ValueError("Supply a configuration file or settings object, not both.")
    if path is not None:
        supplied = json.loads(Path(path).read_text(encoding="utf-8"))
    if supplied is not None:
        if not isinstance(supplied, dict):
            raise ValueError("Configuration must be a JSON object.")
        for key, value in supplied.items():
            if key not in settings:
                raise ValueError(f"Unknown configuration setting: {key}")
            if isinstance(settings[key], dict):
                if not isinstance(value, dict):
                    raise ValueError(f"{key} must be a JSON object.")
                unknown = value.keys() - settings[key].keys()
                if unknown:
                    raise ValueError(f"Unknown {key} settings: {', '.join(sorted(unknown))}")
                settings[key].update(value)
            else:
                settings[key] = value
    return settings


def validate_plate_settings(settings):
    """Check settings before mesh creation, material updates, or output writes."""
    geometry, loading, quadrature = (settings[key] for key in ("geometry", "loading", "quadrature"))
    for name, value in (("nx", geometry["nx"]), ("ny", geometry["ny"]),
                        ("steps_per_phase", loading["steps_per_phase"]), ("depth", quadrature["depth"])):
        minimum = 0 if name == "depth" else 1
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}.")
    numeric = [geometry["length"], geometry["height"], geometry["radius"], settings["thickness"],
               loading["peak_displacement"], loading["loading_time"], loading["hold_time"], quadrature["alpha"]]
    if not np.all(np.isfinite(numeric)):
        raise ValueError("Geometry, thickness, loading and alpha must be finite numbers.")
    if (geometry["length"] <= 0 or geometry["height"] <= 0
            or not 0 <= geometry["radius"] < min(geometry["length"], geometry["height"])):
        raise ValueError("Require positive dimensions and 0 <= radius < min(length, height).")
    if settings["thickness"] <= 0 or loading["loading_time"] <= 0 or loading["hold_time"] <= 0:
        raise ValueError("Thickness and phase durations must be positive.")
    if not 0 < quadrature["alpha"] <= 1 or quadrature["method"] not in ("quadtree", "moment-fitting"):
        raise ValueError("Require 0 < alpha <= 1 and method 'quadtree' or 'moment-fitting'.")
    J2Material(**settings["material"])


def settings_from_metrics(metrics):
    """Read both configurable runs and the original fixed-square v3 snapshots."""
    if "settings" in metrics:
        return deepcopy(metrics["settings"])
    settings = deepcopy(DEFAULT_SETTINGS)
    settings["geometry"] = {"nx": metrics["mesh_n"], "ny": metrics["mesh_n"],
                            "length": metrics["plate_side_mm"], "height": metrics["plate_side_mm"],
                            "radius": metrics["hole_radius_mm"]}
    settings["material"] = metrics["material"].copy()
    settings["thickness"] = metrics["thickness_mm"]
    settings["loading"].update(peak_displacement=metrics["peak_displacement_mm"],
                                steps_per_phase=metrics["steps_per_phase"])
    settings["quadrature"] = {key: metrics[key] for key in ("method", "depth", "alpha")}
    return settings
