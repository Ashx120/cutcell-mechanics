# V4 walkthrough

V4 is the `v0.4.0` source milestone, developed on `v4-thermomechanics`.
The broader beginner-facing cleanup will follow as a separate final-release stage.
Run commands from the repository root using the installed virtual environment.
The root `playground.py` still runs v3. Use `thermal_playground.py` for v4.

## Edit v4 inputs in one place

```powershell
.\.venv\Scripts\python.exe thermal_playground.py
```

Edit its settings block and rerun the same command:

| Setting | Meaning |
|---|---|
| `GEOMETRY` | Independent nx/ny, length, height, and hole radius in mm |
| `MATERIAL` | E, nu, initial yield, hardening, and viscosity parameters |
| `THICKNESS`, `CONDUCTIVITY` | mm and W/(mm K) |
| `EXPANSION_COEFFICIENT` | Constant thermal expansion per K |
| `LOADING` | Reference/reservoir temperatures, ramp/hold durations, increments per second |
| `PROPERTIES` | Constant, paper-reference, or paper-temperature elastic/hardening modes |
| `QUADRATURE` | Mechanical method, depth, and fictitious factor alpha |
| `SUPPORTS` | Symmetry rollers or minimal rigid-motion constraints |
| `OUTPUT_DIR`, `OVERWRITE_OUTPUT` | Save location and explicit replacement option |

The default matches the verified 16 x 16 fitted thermal case, with both property
fits active and dt=0.025 s. It usually takes around one to two minutes. Outputs
go to `results/playground/v4/`. Paper property modes replace the corresponding
MATERIAL constants, as explained in the comments; initial yield/viscosity remain
fixed. Settings are validated before solving and saved under `inputs` in
settings.json. CLI overrides remain available through `--help`.

The ramp starts uniformly at reference temperature and heats both reservoirs
toward their specified targets. With the default left target equal to reference,
the left edge stays at 20 C. For custom durations, the runner includes the exact
phase endpoints and uses increments no larger than 1/steps_per_second. A zero
hold duration is supported. Boundary layouts are the two documented choices;
arbitrary boundary conditions require the Python solver interface.

## Pipeline and units

The heat solver computes a stationary nodal temperature field. Q4 shape functions
interpolate temperature to mechanical integration points. The constitutive update
subtracts isotropic thermal strain, evaluates the selected material properties,
and updates trial plastic history. Global Newton iterations enforce mechanical
equilibrium; history is committed only when an increment converges.

Lengths are mm, forces N, stress MPa, time s, and temperatures degrees Celsius.
Temperature differences are K. Conductivity is W/(mm K), heat flow W, and flux
W/mm². The constant conductivity 0.045 corresponds to 45 W/(m K).

## Follow the small checks first

Prefix each script below with `.\.venv\Scripts\python.exe` in PowerShell.

| Script and options | What to check |
|---|---|
| `scripts/thermal_rectangle.py` | Exact linear temperature, Fourier flux, reservoir balance |
| `scripts/check_thermal_expansion.py` | Zero stress under free expansion; -180 MPa under in-plane restraint |
| `scripts/check_temperature_properties.py` | Paper K/G fits, with plastic parameters fixed |
| `scripts/check_temperature_hardening.py` | Paper hardening fit, with initial yield fixed |
| `scripts/thermoelastic_rectangle.py` | Uniform heating: right ux=0.006 mm and top uy=0.0012 mm |
| `scripts/thermoelastic_rectangle.py --case restrained` | Zero expansion and -180 MPa biaxial stress in the assembled model |
| `scripts/thermoelastic_gradient.py` | Exact linear temperature and decreasing displacement/stress errors under refinement |
| `scripts/thermal_viscoplastic_rectangle.py --compare-time-steps` | FEM/material-point agreement and independent temporal sensitivity |
| `scripts/thermal_plate.py` | Heat balance around the approximately insulated hole |
| `scripts/thermoelastic_plate.py` | Temperature-driven displacement, stress, and equilibrium |

Output-producing examples protect nonempty folders. Use a new `--output` path or
explicitly add `--overwrite`. Checks that only print verification values write no
results. Generated walkthrough data are ignored by Git.

## Restrained thermal viscoplastic plate

The quarter-shaped domain is 4 x 4 mm with a 1 mm radius hole at the origin and
1 mm thickness. The physical left edge stays at 20 C. The hot edge ramps from
20 to 320 C over one second, followed by a two-second hold. Intermediate fields
are scaled stationary solutions, not transient heat-conduction solutions.

Mechanical rollers enforce left ux=0 and bottom uy=0. Right/top/hole boundaries
have no prescribed mechanical traction. Quadtree weights approximate the hole;
the fictitious region has alpha-scaled conductivity and elastic stiffness, no
thermal eigenstrain, and no plastic flow. Plastic history is retained only at
physical integration points.

Start with the constant-property checkpoint:

```powershell
.\.venv\Scripts\python.exe scripts/thermal_viscoplastic_plate.py --steps-per-second 40 --output results/walkthrough/v4-constant-dt0025
```

This uses the original E=210000 MPa, nu=0.3, initial yield=250 MPa, Q=100 MPa,
b=20, relaxation time=100 s, exponent=1, and expansion coefficient=12e-6 /K.
The plain Newton solver has no line search or automatic step subdivision;
large temperature jumps can fail. Refine increments rather than treating failed
history as committed.

## Isolate property effects

For elastic softening, compare `--elasticity paper-reference` with
`--elasticity paper-temperature`. Leave hardening at its `constant` default.
The reference case freezes the paper's elastic fits at 20 C so both cases have
the same starting stiffness. Use different output folders and matched time steps.

For hardening softening, keep `--elasticity paper-temperature` in both cases,
and compare `--hardening paper-reference` with `--hardening paper-temperature`.
Both use Q(20 C)=206.3425 MPa and b=1/0.02784959. Initial yield stays at 250 MPa;
this is a partial literature-based model, not the complete paper calibration.
See [material notes](V4_MATERIAL_PROPERTIES.md) for numerical reference results.

## Compare mechanical quadrature

These commands use the same quadtree temperature field and physical settings:

```powershell
.\.venv\Scripts\python.exe scripts/thermal_viscoplastic_plate.py --steps-per-second 40 --elasticity paper-temperature --hardening paper-temperature --method quadtree --output results/walkthrough/v4-thermal-quadtree
.\.venv\Scripts\python.exe scripts/thermal_viscoplastic_plate.py --steps-per-second 40 --elasticity paper-temperature --hardening paper-temperature --method moment-fitting --output results/walkthrough/v4-thermal-fitted
```

Each run saves settings.json, metrics.json, history.csv, nodes.csv,
physical_samples.csv, history.png, and plastic_strain.png. The physical samples
contain final temperature, stress, von Mises stress, and accumulated strain.
Metrics include point counts, constitutive updates, residuals, Newton corrections,
and single-run solve time. Earlier runs may predate some saved metadata.

Sampled maxima use different point locations between quadrature rules. Compare
displacements at shared nodes and history times as well; maxima differences are
not matched-point stress-field error norms. Plot colour scales are automatic and
can differ between runs.

## Diagnose support sensitivity

```powershell
.\.venv\Scripts\python.exe scripts/thermoelastic_plate.py --mesh 32 --right-temperature 320 --supports minimal --output results/walkthrough/v4-minimal-supports
```

This fixes top-left ux/uy and top-right uy, removing rigid motion while permitting
free deformation. It represents a different mechanical setup from quarter-plate
symmetry rollers. The comparison diagnoses restraint-driven stress; it does not
make one boundary condition universally preferable. Absolute displacements use
different rigid-motion references and cannot be compared directly as expansion.

See [results and limitations](V4_RESULTS.md) before interpreting local hotspots
or claiming agreement with the full source paper.
