"""Run CutCell Mechanics. Edit config.py for inputs, or use --config saved/settings.json."""

import argparse
from datetime import datetime
import json
from pathlib import Path

import config
from scripts.thermal_viscoplastic_plate import main as run_analysis
from immersed_elasticity.thermal_settings import settings_from_saved


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--config", type=Path, help="Reuse a saved settings.json (or an inputs-only JSON file).")
    options, remaining = parser.parse_known_args(argv)
    inputs = config.settings()
    if options.config is not None:
        try:
            saved = json.loads(options.config.read_text(encoding="utf-8"))
            inputs = settings_from_saved(saved)
        except (OSError, ValueError, KeyError, TypeError) as error:
            parser.error(f"Could not read configuration: {error}")
    folder = Path(config.OUTPUT_ROOT) / datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    if "--help" in remaining or "-h" in remaining:
        print("CutCell Mechanics: edit config.py, then run python run.py.\n"
              "Use --config PATH to replay a saved configuration.\n"
              "Each normal run has a new output folder; --output PATH selects an explicit folder.\n")
    else:
        print("CutCell Mechanics | heating, deformation, and stress\n"
              "A small-strain, one-way thermal/mechanical demonstration.\n"
              "Inputs: config.py (or --config); new runs preserve previous results.\n", flush=True)
    run_analysis(remaining, initial_settings=inputs, output=folder)


if __name__ == "__main__":
    main()
