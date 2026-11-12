"""Regressions for editable inputs, saved results, and the public script workflow."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from immersed_elasticity.plate_settings import load_plate_settings, validate_plate_settings, settings_from_metrics


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("supplied", [{"unknown": 1}, {"geometry": {"raduis": 1}}, []])
def test_configuration_rejects_unknown_or_invalid_structure(tmp_path, supplied):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(supplied))
    with pytest.raises(ValueError):
        load_plate_settings(path)


def test_configuration_rejects_invalid_mesh_and_does_not_modify_defaults():
    settings = load_plate_settings()
    settings["geometry"]["nx"] = 2.5
    with pytest.raises(ValueError, match="nx"):
        validate_plate_settings(settings)
    assert load_plate_settings()["geometry"]["nx"] == 6


def test_original_snapshot_remains_readable():
    metrics = json.loads((ROOT / "results/v3/quadtree/metrics.json").read_text())
    settings = settings_from_metrics(metrics)
    validate_plate_settings(settings)
    assert settings["geometry"] == dict(nx=16, ny=16, length=4., height=4., radius=1.)
    assert settings["loading"]["steps_per_phase"] == 40


def test_configurable_cli_workflow_and_overwrite_protection(tmp_path):
    config = tmp_path / "input.json"
    config.write_text(json.dumps({
        "geometry": dict(nx=2, ny=1, length=3, height=2, radius=0),
        "thickness": 2,
        "material": {"initial_yield": 1e9},
        "loading": dict(loading_time=.25, hold_time=.5, peak_displacement=.0001, steps_per_phase=2),
    }))
    folders = {}
    for method in ("quadtree", "moment-fitting"):
        output = tmp_path / method
        command = [sys.executable, str(ROOT / "scripts/nonlinear_plate.py"), "--config", str(config),
                   "--method", method, "--steps", "1", "--output", str(output)]
        subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
        history = np.genfromtxt(output / "history.csv", delimiter=",", names=True)
        np.testing.assert_allclose(history["time_s"], [0, .25, .75])
        np.testing.assert_allclose(history["right_force_N"][1:], 28., atol=1e-7)
        metrics = json.loads((output / "metrics.json").read_text())
        assert metrics["settings"]["geometry"]["ny"] == 1
        assert metrics["settings"]["quadrature"]["method"] == method
        assert metrics["steps_per_phase"] == 1
        for name in ("response.png", "plastic_strain.png", "settings.json"):
            assert (output / name).stat().st_size > 0
        before = (output / "history.csv").read_bytes()
        refused = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        assert refused.returncode != 0 and "nonempty" in refused.stderr
        assert (output / "history.csv").read_bytes() == before
        folders[method] = output
    comparison = tmp_path / "comparison"
    subprocess.run([sys.executable, str(ROOT / "scripts/compare_nonlinear_histories.py"),
                    "--reference", str(folders["quadtree"]), "--fitted", str(folders["moment-fitting"]),
                    "--output", str(comparison)], cwd=ROOT, check=True, capture_output=True, text=True)
    summary = json.loads((comparison / "comparison.json").read_text())
    assert summary["max_absolute_force_difference_N"] < 1e-8
    assert (comparison / "reaction_comparison.png").is_file()


def test_v3_playground_uses_its_editable_settings(tmp_path, monkeypatch):
    import playground

    monkeypatch.setattr(playground, "GEOMETRY", dict(nx=2, ny=1, length=3, height=2, radius=0))
    monkeypatch.setattr(playground, "THICKNESS", 2.)
    monkeypatch.setattr(playground, "MATERIAL", {**playground.MATERIAL, "initial_yield": 1e9})
    monkeypatch.setattr(playground, "LOADING", dict(peak_displacement=.0001, loading_time=.25,
                                                 hold_time=.5, steps_per_phase=1))
    monkeypatch.setattr(playground, "OUTPUT_DIR", tmp_path / "playground")
    playground.main([])
    saved = json.loads((playground.OUTPUT_DIR / "settings.json").read_text())
    assert saved["geometry"] == playground.GEOMETRY
    assert saved["loading"] == playground.LOADING
    history = np.genfromtxt(playground.OUTPUT_DIR / "history.csv", delimiter=",", names=True)
    np.testing.assert_allclose(history["right_force_N"][1:], 28, atol=1e-7)
    with pytest.raises(SystemExit):
        playground.main([])
    playground.main(["--steps", "2", "--output", str(tmp_path / "override")])
    assert playground.LOADING["steps_per_phase"] == 1
    assert json.loads((tmp_path / "override/settings.json").read_text())["loading"]["steps_per_phase"] == 2


def test_relocated_elastic_playground_remains_runnable():
    from scripts.elastic_playground import run_model

    (mesh, material, displacement), energy = run_model(dict(nx=2, ny=1, length=3, height=2, radius=0))
    np.testing.assert_allclose(displacement[0::2], 100 / material.young * mesh.nodes[:, 0], atol=1e-12)
    assert energy > 0
