# Contributing

Contributions should improve numerical correctness, reproducibility, or the
experience of running and understanding the examples. For changes to the
constitutive model, element formulation, or supported loading conditions, open
an issue to describe the equations, assumptions, and intended verification first.

## Development setup

Fork the repository, clone your fork, and create a branch from the intended base.
Until the v3 work is merged, its development base is `v3-viscoplasticity`.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pip check
```

On Linux/macOS use `.venv/bin/python`. The editable installation makes source
changes available without reinstalling. After changing package metadata, repeat
the install command. CI tests Windows and Linux with Python 3.11 and 3.12 and
builds an installable wheel.

## Where changes belong

- `playground.py` is the primary editable v3 interface. It calls the shared runner.
- `src/immersed_elasticity/` contains reusable numerical and input-handling code.
- `scripts/` contains thin runnable experiments and comparisons.
- `examples/` contains optional command-line input configurations.
- `tests/` contains independent physical checks and regression tests.
- `docs/` contains detailed verification reports and the roadmap.
- `results/v3/` contains selected published reference snapshots.

Keep temporary runs in `results/walkthrough/` or `results/playground/`, which are
ignored by Git. Use a new output directory to preserve earlier runs.

## Numerical changes

Explain the problem, expected behaviour, and units. Add a focused regression
test for changed numerical behaviour. Prefer analytical limits, finite-difference
derivatives, and independent equilibrium checks over duplicating the same code
in a test. Check elastic behaviour before assessing a nonlinear extension.

Previous material history must remain unchanged during trial iterations and
after failed increments. Physical and fictitious domains must retain their
distinct material behaviour. New meshes or quadrature rules require fresh
histories unless explicit state-transfer support is implemented and verified.

Report timing with the exact settings, environment, measured scope, and number
of trials. Keep hardware-specific performance thresholds out of pass/fail tests.
Use the saved v3 data to reproduce comparisons without expensive fresh solves:

```powershell
.\.venv\Scripts\python.exe scripts/compare_nonlinear_histories.py --reference results/v3/quadtree --fitted results/v3/fitted --output results/walkthrough/review-comparison
```

## Pull requests

Describe what changed, why it matters, and how you verified it. Include model
inputs or a minimal reproduction for numerical bugs. Update the README when
entry points change and the relevant report when published results change.
Record user-facing changes under `Unreleased` in `CHANGELOG.md`.

Contributions are made under the project's MIT license. Keep discussion
respectful, specific, and focused on evidence and the work.
