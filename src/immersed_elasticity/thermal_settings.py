"""Validated inputs shared by the v4 playground and command-line runner."""

from copy import deepcopy

import numpy as np

from .solver import structured_mesh
from .viscoplasticity import J2Material

DEFAULTS = dict(
    geometry=dict(nx=16, ny=16, length=4., height=4., radius=1.),
    material=dict(young=210000., poisson=.3, initial_yield=250., saturation_increment=100.,
                  hardening_rate=20., relaxation_time=100., rate_exponent=1.),
    thickness=1., conductivity=.045, expansion_coefficient=12e-6,
    loading=dict(reference_temperature=20., left_temperature=20., right_temperature=320.,
                 heating_time=1., hold_time=2., steps_per_second=20),
    quadrature=dict(method="quadtree", depth=3, alpha=1e-6),
    properties=dict(elasticity="constant", hardening="constant"),
    supports="symmetry")


def thermal_settings(supplied=None):
    """Merge explicit input with defaults, rejecting misspelled/invalid settings."""
    result = deepcopy(DEFAULTS)
    if supplied is not None:
        if not isinstance(supplied, dict) or set(supplied)-set(result):
            raise ValueError("Unknown thermal setting or invalid settings structure.")
        for key, value in supplied.items():
            if isinstance(result[key], dict):
                if not isinstance(value, dict) or set(value)-set(result[key]):
                    raise ValueError(f"Unknown or invalid thermal settings in {key}.")
                result[key].update(deepcopy(value))
            else:
                result[key] = deepcopy(value)
    g, loading, q = result["geometry"], result["loading"], result["quadrature"]
    for key in ("nx", "ny"):
        if isinstance(g[key], bool) or not isinstance(g[key], (int, np.integer)) or g[key] < 1:
            raise ValueError(f"Require a positive integer geometry.{key}.")
    # Validate dimensions without allocating the requested full mesh.
    structured_mesh(1, 1, length=g["length"], height=g["height"], radius=g["radius"])
    J2Material(**result["material"])
    for key in ("thickness", "conductivity"):
        if not np.isfinite(result[key]) or result[key] <= 0:
            raise ValueError(f"Require positive finite {key}.")
    if not np.isfinite(result["expansion_coefficient"]):
        raise ValueError("Require finite expansion_coefficient.")
    for key in ("reference_temperature", "left_temperature", "right_temperature"):
        if not np.isfinite(loading[key]):
            raise ValueError(f"Require finite loading.{key}.")
    if (not np.isfinite(loading["heating_time"]) or loading["heating_time"] <= 0
            or not np.isfinite(loading["hold_time"]) or loading["hold_time"] < 0):
        raise ValueError("Require positive heating_time and nonnegative hold_time.")
    rate = loading["steps_per_second"]
    if isinstance(rate, bool) or not isinstance(rate, (int, np.integer)) or rate < 1:
        raise ValueError("Require positive integer steps_per_second.")
    if (isinstance(q["depth"], bool) or not isinstance(q["depth"], (int, np.integer))
            or q["depth"] < 0 or not np.isfinite(q["alpha"]) or not 0 < q["alpha"] <= 1):
        raise ValueError("Require integer depth >= 0 and 0 < alpha <= 1.")
    if q["method"] not in ("quadtree", "moment-fitting"):
        raise ValueError("Unknown mechanical quadrature method.")
    for mode in result["properties"].values():
        if mode not in ("constant", "paper-reference", "paper-temperature"):
            raise ValueError("Unknown temperature-property mode.")
    if "paper-temperature" in result["properties"].values():
        if not all(0 <= loading[key] <= 1500 for key in
                   ("reference_temperature", "left_temperature", "right_temperature")):
            raise ValueError("Temperature-property fits require temperatures in 0-1500 C.")
    if result["supports"] not in ("symmetry", "minimal"):
        raise ValueError("Supports must be symmetry or minimal.")
    return result


def settings_from_saved(saved):
    """Read current inputs or the flat settings format from the v4 checkpoints."""
    if not isinstance(saved, dict):
        raise ValueError("Saved settings must be a JSON object.")
    if "inputs" in saved:
        return thermal_settings(saved["inputs"])
    if "conductivity_W_per_mm_K" not in saved:
        return thermal_settings(saved)
    supports = saved.get("supports", "symmetry")
    if supports == "left ux=0; bottom uy=0":
        supports = "symmetry"
    reference = saved["reference_temperature_C"]
    return thermal_settings(dict(
        geometry=saved["geometry"], material=saved["material"],
        thickness=saved["thickness_mm"], conductivity=saved["conductivity_W_per_mm_K"],
        expansion_coefficient=saved["expansion_coefficient_per_K"], supports=supports,
        loading=dict(reference_temperature=reference, left_temperature=saved.get("left_temperature_C", reference),
                     right_temperature=saved["right_temperature_C"], heating_time=saved["heating_time_s"],
                     hold_time=saved["hold_time_s"], steps_per_second=saved["steps_per_second"]),
        quadrature=dict(method=saved.get("mechanical_method", "quadtree"), depth=saved["depth"], alpha=saved["alpha"]),
        properties=dict(elasticity=saved.get("elasticity", "constant"), hardening=saved.get("hardening", "constant"))))
