# Changelog

User-facing changes are grouped by source release and development milestone.
Version 0.3.0 records the verified v3 source baseline.

## Unreleased

## [0.3.0] - 2026-11-12

### Added

- Isothermal J2 Perzyna viscoplasticity with exponential saturation hardening.
- Backward-Euler stress update, analytical tangent, and plane-stress condensation.
- Incremental nonlinear equilibrium with trial and committed material history.
- Separate physical/fictitious moment fitting for history-dependent analysis.
- Material-point, rectangle, and plate experiments, plus full-history comparisons.
- Reference reaction histories and recorded sensitivity studies.
- Editable v3 playground, optional JSON inputs, saved settings, and input validation.
- Contribution guidance, software citation metadata, issue forms, and PR template.
- Cross-platform CI configuration and wheel build/import verification.
- Brief AI-assistance statement in the README.

### Changed

- Public project title changed to CutCell Mechanics. The Python package retains
  its existing import name, `immersed_elasticity`.
- Root `playground.py` now runs v3. The elastic playground moved to
  `scripts/elastic_playground.py`.
- Nonempty nonlinear analysis output folders require explicit overwrite permission.
- README follows installation, editable inputs, execution, output interpretation,
  and verification. Detailed elastic results moved to `docs/LINEAR_RESULTS.md`.

## V2 milestone

- Non-negative polynomial moment fitting alongside quadtree integration.
- Moment, stiffness, displacement, and reaction comparisons on the Kirsch plate.
- Setup and cached-integration timing experiments.

## V1 milestone

- Q4 plane-stress elasticity and immersed circular-hole geometry.
- Quadtree cut-cell integration, sparse assembly, and displacement constraints.
- Patch, uniaxial-tension, area, equilibrium, and Kirsch convergence checks.
- The original milestone is accessible through the `v1.0.0` Git tag.
