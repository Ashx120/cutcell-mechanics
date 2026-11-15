"""Run CutCell Mechanics. Edit config.py for inputs, or use --config saved/settings.json."""

import argparse
from datetime import datetime
import json
from pathlib import Path

import config
from scripts.thermal_viscoplastic_plate import main as run_thermal
from scripts.nonlinear_plate import main as run_displacement
from immersed_elasticity.thermal_settings import settings_from_saved
from immersed_elasticity.plate_settings import load_plate_settings


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--config", type=Path, help="Reuse a saved settings.json (or an inputs-only JSON file).")
    parser.add_argument("--mode", choices=("thermal", "displacement"), help="Override MODE in config.py.")
    options, remaining = parser.parse_known_args(argv)
    mode = options.mode or config.MODE
    try:
        inputs = config.settings(options.mode) if options.mode else config.settings()
    except ValueError as error:
        parser.error(str(error))
    if options.config is not None:
        try:
            saved = json.loads(options.config.read_text(encoding="utf-8"))
            if not isinstance(saved, dict):
                raise ValueError("Configuration must be a JSON object.")
            payload = saved.get("inputs", saved.get("settings", saved))
            if not isinstance(payload, dict):
                raise ValueError("Saved inputs must be a JSON object.")
            loading = payload.get("loading", {})
            detected = ("displacement" if "peak_displacement" in loading else
                        "thermal" if "reference_temperature" in loading or "conductivity_W_per_mm_K" in saved else mode)
            if options.mode and options.mode != detected:
                raise ValueError("--mode conflicts with the saved loading type.")
            mode = detected
            inputs = settings_from_saved(saved) if mode == "thermal" else load_plate_settings(supplied=payload)
        except (OSError, ValueError, KeyError, TypeError) as error:
            parser.error(f"Could not read configuration: {error}")
    folder = Path(config.OUTPUT_ROOT) / (mode+"-"+datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
    if "--help" in remaining or "-h" in remaining:
        print("CutCell Mechanics: edit config.py, then run python run.py.\n"
              "Use --config PATH to replay a saved configuration.\n"
              "Set MODE in config.py or use --mode thermal/displacement.\n"
              "Each normal run has a new output folder; --output PATH selects an explicit folder.\n")
    else:
        print(f"CutCell Mechanics | {mode} loading\n"
              "Inputs: config.py (or --config); new runs preserve previous results.\n", flush=True)
    runner = run_thermal if mode == "thermal" else run_displacement
    result_folder = runner(remaining, initial_settings=inputs, output=folder)
    if mode == "displacement":
        (result_folder / "summary.txt").write_text(
            "CutCell Mechanics - displacement loading\n\n"
            "The right edge moves horizontally, then holds. Left ux=0, bottom uy=0.\n"
            "Open response.png for reaction force and displacement histories.\n"
            "Open plastic_strain.png for final accumulated plastic strain.\n"
            "Settings are in settings.json; metrics and history.csv contain numerical results.\n"
            "This mode is isothermal: thermal settings in config.py do not apply.\n"
            "A small equilibrium residual does not establish mesh convergence.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
