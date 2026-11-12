"""V3 playground. Edit the settings below, then run: python playground.py."""

from pathlib import Path

from scripts.nonlinear_plate import main as run_analysis


# Geometry in mm. The circular hole stays centred at (0, 0).
# Use radius=0 for a solid rectangle.
GEOMETRY = dict(nx=6, ny=6, length=4.0, height=4.0, radius=1.0)
THICKNESS = 1.0  # mm

# Isothermal J2 material. Stress parameters are in MPa, relaxation time in s.
# These are illustrative parameters, not calibrated to a specific alloy.
MATERIAL = dict(
    young=210_000.0,
    poisson=0.3,
    initial_yield=250.0,
    saturation_increment=100.0,
    hardening_rate=20.0,
    relaxation_time=100.0,
    rate_exponent=1.0,
)

# Increase right-edge ux from zero to peak_displacement, then hold it fixed.
# Left ux=0, bottom uy=0. Right uy is free; top and hole are unloaded.
LOADING = dict(
    peak_displacement=0.012,  # mm
    loading_time=1.0,        # seconds
    hold_time=1.0,           # seconds
    steps_per_phase=10,
)

# Choose "quadtree" or "moment-fitting". Only cut cells are compressed.
QUADRATURE = dict(method="quadtree", depth=2, alpha=1e-6)

# Each run writes settings, CSV histories, plots, and timing diagnostics here.
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "playground" / "v3"
# Set True to replace an earlier run, or choose a different OUTPUT_DIR.
OVERWRITE_OUTPUT = False


def main(argv=None):
    settings = dict(geometry=GEOMETRY, material=MATERIAL, thickness=THICKNESS,
                    loading=LOADING, quadrature=QUADRATURE)
    run_analysis(argv, initial_settings=settings, output=OUTPUT_DIR, overwrite=OVERWRITE_OUTPUT)


if __name__ == "__main__":
    main()
