# Technical and development guide

New users should start with [GETTING_STARTED.md](GETTING_STARTED.md).
The public path is `run.py` -> `config.py` -> the selected shared runner.
MODE chooses thermal or displacement loading. The two former root playgrounds
are merged into config.py; their original files remain in the milestone tags.

## Code map

| Location | Responsibility |
|---|---|
| `run.py`, `config.py` | Public execution, saved-input replay, editable inputs |
| `src/immersed_elasticity/thermal_settings.py` | Input validation and v4 saved-run compatibility |
| `solver.py` | Q4 geometry/kinematics, linear assembly, constraints, solve |
| `quadrature.py`, `moment_fitting.py` | Quadtree reference integration and non-negative compression |
| `viscoplasticity.py`, `plane_stress.py` | Constitutive update, history, algorithmic tangents |
| `nonlinear_solver.py` | Incremental global equilibrium and history commit |
| `thermal.py`, `thermal_expansion.py`, `thermoelastic.py` | Stationary conduction and thermal coupling |
| `temperature_properties.py` | Selected literature property fits |
| `boundary_conditions.py` | Shared thermal mechanical support layouts |
| `thermal_plots.py` | Physical-sample field and history visualisation |
| `scripts/` | Shared experiment runners, verification, and comparisons |
| `tests/` | Physical limits, numerical checks, and user-workflow regressions |
| `results/reference/thermal/` | Selected saved thermal comparison with provenance |
| `results/v3/` | Earlier published isothermal verification snapshots |

Source-module names refer to `src/immersed_elasticity/`. The package keeps its
original import name for compatibility. Geometry and loading are example-specific;
the Python solver interfaces support custom constraints but are not a general CAE GUI.

## Run the independent checks

Prefix script names with `.\.venv\Scripts\python.exe` on Windows.

| Script | Purpose |
|---|---|
| `scripts/material_point.py` | Shear loading, relaxation, and unloading |
| `scripts/check_material_tangent.py` | 3D tangent check |
| `scripts/check_plane_stress.py` | Zero thickness stress and condensed tangent |
| `scripts/nonlinear_rectangle.py` | Uniform FEM/material-point agreement |
| `scripts/thermal_rectangle.py` | Analytical temperature/flux and heat balance |
| `scripts/check_thermal_expansion.py` | Free/restrained material-point expansion |
| `scripts/thermoelastic_rectangle.py` | Free thermal FEM expansion; add `--case restrained` |
| `scripts/thermoelastic_gradient.py` | Compatible gradient and mesh convergence |
| `scripts/thermal_viscoplastic_rectangle.py --compare-time-steps` | Thermal history and time sensitivity |
| `scripts/check_temperature_properties.py` | Elastic property-fit checkpoint |
| `scripts/check_temperature_hardening.py` | Hardening property-fit checkpoint |
| `scripts/compare_thermal_runs.py` | Shared-node/shared-time thermal comparison |

Use `--help` for script options. Set MODE in config.py for thermal or isothermal
displacement loading. `scripts/elastic_playground.py` remains an advanced
linear-elastic verification example.

For earlier isothermal comparisons, use `scripts/nonlinear_plate.py` with
`examples/nonlinear_plate.json`, then `scripts/compare_nonlinear_histories.py`.
For elastic quadrature studies use `scripts/run_benchmark.py`, `fit_one_cell.py`,
`compare_quadrature.py`, and `time_quadrature.py` under scripts/.

## Development milestones and further research

The original elastic milestone is retained as `v1.0.0`; the later v3/v4 development
snapshots are `v0.3.0` and `v0.4.0`. The final presentation release uses `v1.1.0`
to preserve those existing tags. Development history is not a requirement for
running the current example.

The [v4 walkthrough](V4_WALKTHROUGH.md) and [material notes](V4_MATERIAL_PROPERTIES.md)
retain the follow-along details. Future numerical research would include higher-
order approximation, adaptive overlays, paper-matched benchmarks, and missing
property laws. None are required to use the declared scope of this release.

See [CONTRIBUTING.md](../CONTRIBUTING.md) for development setup and the
[publishing checklist](PUBLISHING.md) for release verification. Add independent
checks before changing constitutive equations or history management.
