"""V4 entry point. Edit this settings block, then run: python thermal_playground.py."""

from pathlib import Path

from scripts.thermal_viscoplastic_plate import main as run_analysis

# Dimensions in mm; quarter-circle hole centred at (0, 0). Use radius=0 for solid.
GEOMETRY = dict(nx=20, ny=10, length=4., height=4., radius=1.)
THICKNESS = 1.  # mm
CONDUCTIVITY = .045  # W/(mm K)
EXPANSION_COEFFICIENT = 12e-6  # per K; constant in this milestone

# Stress in MPa and relaxation time in seconds. Initial yield/viscosity stay fixed.
MATERIAL = dict(young=210000., poisson=.3, initial_yield=250.,
                saturation_increment=100., hardening_rate=20.,
                relaxation_time=100., rate_exponent=1.)

# Start uniformly at reference temperature; ramp to left/right reservoir targets,
# then hold. Fields are equilibrium snapshots, not transient heat diffusion.
LOADING = dict(reference_temperature=20., left_temperature=20., right_temperature=320.,
               heating_time=1., hold_time=2., steps_per_second=40)

# constant: use MATERIAL values. paper-reference: freeze the paper fit at 20 C.
# paper-temperature: evaluate the fit at each physical point's current temperature.
# Paper hardening modes replace MATERIAL Q/b; initial yield remains supplied above.
PROPERTIES = dict(elasticity="paper-temperature", hardening="paper-temperature")
QUADRATURE = dict(method="moment-fitting", depth=3, alpha=1e-6)

# symmetry: left ux=0, bottom uy=0 (the verified restrained-heating example).
# minimal: remove only rigid motion; this changes the physical mechanical setup.
SUPPORTS = "symmetry"
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "playground" / "v4"
OVERWRITE_OUTPUT = False


def main(argv=None):
    settings = dict(geometry=GEOMETRY, material=MATERIAL, thickness=THICKNESS,
                    conductivity=CONDUCTIVITY, expansion_coefficient=EXPANSION_COEFFICIENT,
                    loading=LOADING, properties=PROPERTIES, quadrature=QUADRATURE, supports=SUPPORTS)
    run_analysis(argv, initial_settings=settings, output=OUTPUT_DIR, overwrite=OVERWRITE_OUTPUT)


if __name__ == "__main__":
    main()
