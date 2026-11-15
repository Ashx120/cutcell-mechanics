"""Public entry-point checks: safe reruns, saved-input replay, and useful outputs."""

import json
from pathlib import Path

import numpy as np
import pytest

from immersed_elasticity.thermal_settings import thermal_settings
from scripts.compare_thermal_runs import compare


def test_public_runner_preserves_runs_and_replays_saved_inputs(tmp_path, monkeypatch):
    import config
    import run

    inputs = thermal_settings(dict(geometry=dict(nx=2, ny=1, radius=0.),
                                  properties=dict(elasticity="constant", hardening="constant"),
                                  loading=dict(left_temperature=70., right_temperature=70.,
                                               heating_time=.25, hold_time=0., steps_per_second=4)))
    monkeypatch.setattr(config, "settings", lambda: inputs)
    monkeypatch.setattr(config, "OUTPUT_ROOT", tmp_path / "runs")
    run.main([])
    first = next(config.OUTPUT_ROOT.iterdir())
    original = (first / "history.csv").read_bytes()
    run.main([])
    assert len(list(config.OUTPUT_ROOT.iterdir())) == 2
    assert (first / "history.csv").read_bytes() == original
    for name in ("summary.txt", "fields.png", "history.png", "settings.json", "metrics.json"):
        assert (first / name).stat().st_size > 0
    saved = json.loads((first / "settings.json").read_text())
    assert saved["inputs"] == inputs
    replay = tmp_path / "replay"
    run.main(["--config", str(first / "settings.json"), "--output", str(replay)])
    a = np.genfromtxt(first / "nodes.csv", delimiter=",", names=True)
    b = np.genfromtxt(replay / "nodes.csv", delimiter=",", names=True)
    for key in a.dtype.names:
        np.testing.assert_allclose(a[key], b[key], atol=1e-12)
    np.testing.assert_allclose(a["ux_mm"], .0006*a["x_mm"], atol=1e-12)
    np.testing.assert_allclose(a["uy_mm"], .0006*a["y_mm"], atol=1e-12)
    with pytest.raises(SystemExit):
        run.main(["--config", str(first / "settings.json"), "--output", str(replay)])


def test_public_runner_reports_unreadable_configuration(tmp_path):
    import run

    with pytest.raises(SystemExit):
        run.main(["--config", str(tmp_path / "missing.json")])


def test_published_thermal_snapshots_remain_reproducibly_comparable():
    base = Path(__file__).resolve().parents[1] / "results/reference/thermal"
    summary, _, _ = compare(base / "quadtree", base / "fitted")
    assert summary["max_mean_right_displacement_difference_mm"] < 3e-9
    assert summary["final_relative_nodal_displacement_norm"] < 2e-5


def test_thermal_comparison_rejects_mismatched_physical_inputs(tmp_path):
    base = Path(__file__).resolve().parents[1] / "results/reference/thermal"
    changed = json.loads((base / "quadtree/settings.json").read_text())
    changed["right_temperature_C"] = 400.
    folder = tmp_path / "different"
    folder.mkdir()
    (folder / "settings.json").write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="Physical inputs differ"):
        compare(folder, base / "fitted")
