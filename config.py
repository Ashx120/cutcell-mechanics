"""Edit this file, then run: python run.py. Units: mm, N, MPa, seconds, degrees C."""

from pathlib import Path

# Choose "thermal" for heating or "displacement" for tensile loading.
# Both modes use the shared geometry, material, thickness, and quadrature below.
MODE = "thermal"

# A coarse, quick demonstration. For the recorded comparison use 16 x 16 and
# 40 steps per second. Increasing resolution improves accuracy and takes longer.
GEOMETRY = dict(nx=8, ny=8, length=4., height=4., radius=1.)
THICKNESS = 1.  # mm
CONDUCTIVITY = .045  # W/(mm K); equivalent to 45 W/(m K)
EXPANSION_COEFFICIENT = 12e-6  # Free thermal strain per degree of temperature rise

# Stress in MPa; relaxation_time in seconds. These are demonstration parameters,
# not a complete calibrated alloy. "young" controls elastic stiffness, "poisson"
# controls lateral contraction; the remaining parameters control plastic flow.
MATERIAL = dict(young=210000., poisson=.3, initial_yield=250.,
                saturation_increment=100., hardening_rate=20.,
                relaxation_time=100., rate_exponent=1.)

# Begin uniformly at reference_temperature. Ramp to the left/right targets,
# then hold them. These are stationary temperature snapshots, not heat diffusion.
THERMAL_LOADING = dict(reference_temperature=20., left_temperature=20., right_temperature=320.,
               heating_time=1., hold_time=2., steps_per_second=20)

# Displacement mode: move the right edge horizontally, then hold it fixed.
# Left ux=0 and bottom uy=0; right uy is free. No temperature is applied.
DISPLACEMENT_LOADING = dict(peak_displacement=.012, loading_time=1.,
                            hold_time=1., steps_per_phase=10)

# Thermal mode only. Displacement mode uses the constant MATERIAL values above.
# constant: use MATERIAL values above.
# paper-reference: use the paper fit frozen at 20 C.
# paper-temperature: evaluate that fit at the local temperature.
# Paper hardening replaces saturation_increment/hardening_rate; initial_yield
# and viscosity stay fixed. See docs/MODEL.md for the scope of these fits.
# Start with constant properties so edits to MATERIAL take effect directly.
# The recorded paper-fit comparison is reproducible through saved settings.
PROPERTIES = dict(elasticity="constant", hardening="constant")

# moment-fitting uses fewer mechanical integration points than quadtree.
# depth resolves the curved hole for integration; it does not refine the mesh.
# alpha gives the fictitious region a tiny conductivity/stiffness for stability.
QUADRATURE = dict(method="moment-fitting", depth=3, alpha=1e-6)

# Thermal mode only. symmetry: rollers on left (no x motion) and bottom (no y motion). This is the
# restrained-heating example. minimal: only prevent whole-body motion; stresses
# will be much lower. Changing supports changes the physical problem.
SUPPORTS = "symmetry"

# Every run gets a new timestamped subfolder here. Old runs are kept automatically.
OUTPUT_ROOT = Path(__file__).resolve().parent / "results" / "runs"


def settings(mode=None):
    selected = MODE if mode is None else mode
    if selected == "displacement":
        return dict(geometry=GEOMETRY, material=MATERIAL, thickness=THICKNESS,
                    loading=DISPLACEMENT_LOADING, quadrature=QUADRATURE)
    if selected != "thermal":
        raise ValueError("MODE must be 'thermal' or 'displacement'.")
    return dict(geometry=GEOMETRY, material=MATERIAL, thickness=THICKNESS,
                conductivity=CONDUCTIVITY, expansion_coefficient=EXPANSION_COEFFICIENT,
                loading=THERMAL_LOADING, properties=PROPERTIES, quadrature=QUADRATURE, supports=SUPPORTS)
