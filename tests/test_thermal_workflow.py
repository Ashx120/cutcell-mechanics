"""Verify the editable v4 entry point, validation, and saved experiment inputs."""

import json

import numpy as np
import pytest

from immersed_elasticity.thermal_settings import thermal_settings


@pytest.mark.parametrize("supplied", [
    {"unknown": 1}, {"geometry": {"raduis": 1}}, {"geometry": {"nx": 2.5}},
    {"loading": {"heating_time": 0}}, {"loading": {"hold_time": -1}},
    {"properties": {"hardening": "paper-temperature"}, "loading": {"right_temperature": 1600}},
])
def test_thermal_inputs_reject_invalid_settings(supplied):
    with pytest.raises(ValueError):
        thermal_settings(supplied)


def test_thermal_defaults_are_independent_between_calls():
    settings = thermal_settings()
    settings["geometry"]["nx"] = 99
    assert thermal_settings()["geometry"]["nx"] == 16


def test_thermal_mode_uses_editable_inputs_and_protects_results(tmp_path, monkeypatch):
    import config as playground
    import run
    monkeypatch.setattr(playground, "MODE", "thermal")

    monkeypatch.setattr(playground, "GEOMETRY", dict(nx=2, ny=1, length=3., height=2., radius=0.))
    monkeypatch.setattr(playground, "MATERIAL", dict(young=100000., poisson=.2, initial_yield=1e6,
                        saturation_increment=0., hardening_rate=0., relaxation_time=1., rate_exponent=1.))
    monkeypatch.setattr(playground, "THICKNESS", 2.)
    monkeypatch.setattr(playground, "CONDUCTIVITY", .03)
    monkeypatch.setattr(playground, "EXPANSION_COEFFICIENT", 10e-6)
    monkeypatch.setattr(playground, "THERMAL_LOADING", dict(reference_temperature=20., left_temperature=70.,
                        right_temperature=70., heating_time=.25, hold_time=.5, steps_per_second=4))
    monkeypatch.setattr(playground, "PROPERTIES", dict(elasticity="constant", hardening="constant"))
    monkeypatch.setattr(playground, "QUADRATURE", dict(method="quadtree", depth=2, alpha=1e-6))
    monkeypatch.setattr(playground, "SUPPORTS", "symmetry")
    output = tmp_path / "thermal"
    run.main(["--output", str(output)])
    saved = json.loads((output / "settings.json").read_text())
    inputs = saved["inputs"]
    assert inputs["geometry"] == playground.GEOMETRY
    assert inputs["material"] == playground.MATERIAL
    assert inputs["loading"] == playground.THERMAL_LOADING
    assert inputs["conductivity"] == .03 and inputs["thickness"] == 2.
    history = np.genfromtxt(output / "history.csv", delimiter=",", names=True)
    np.testing.assert_allclose(history["time_s"], [0, .25, .5, .75])
    nodes = np.genfromtxt(output / "nodes.csv", delimiter=",", names=True)
    np.testing.assert_allclose(nodes["ux_mm"], .0005*nodes["x_mm"], atol=1e-13)
    np.testing.assert_allclose(nodes["uy_mm"], .0005*nodes["y_mm"], atol=1e-13)
    np.testing.assert_allclose(history["max_von_mises_MPa"], 0, atol=1e-8)
    for name in ("settings.json", "metrics.json", "history.png", "plastic_strain.png", "physical_samples.csv"):
        assert (output / name).stat().st_size > 0
    before = (output / "history.csv").read_bytes()
    with pytest.raises(SystemExit):
        run.main(["--output", str(output)])
    assert (output / "history.csv").read_bytes() == before
    run.main(["--mesh", "2", "--supports", "minimal", "--output", str(tmp_path / "override")])
    overridden = json.loads((tmp_path / "override/settings.json").read_text())["inputs"]
    assert overridden["geometry"]["ny"] == 2 and overridden["supports"] == "minimal"
    assert playground.GEOMETRY["ny"] == 1
