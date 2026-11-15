# CutCell Mechanics 1.1.0

This release packages the verified v4 numerical scope as a usable, documented
portfolio project. Start with `python run.py`, edit `config.py`, and inspect the
automatically saved summary, history, and physical-field plots.

The default is a coarse, quick restrained-heating demonstration with constant
properties. Saved configurations replay the finer temperature-dependent cases.
Independent run folders preserve previous results. The README and first-run
guide introduce the engineering problem before the technical details.

Selected thermal reference data, provenance, and a comparison command accompany
the results. In the recorded matched case, moment fitting uses 33.5% fewer total
mechanical points and 17.0% fewer physical history points, with close nodal and
history agreement. These are case-specific observations, not general guarantees.

Verification includes 93 automated checks, successful starter/reference runs,
dependency verification, and package wheel build/import checks. GitHub CI covers
Windows/Linux on Python 3.11/3.12 plus wheel build/import.

No new physics beyond v4 is added. This is a paper-inspired partial implementation:
no higher-order adaptive overlays, matched 3D paper benchmark, transient heat,
plastic heating, or complete temperature-dependent initial-yield calibration.
The plastic hotspot is restraint dependent and is not established as mesh converged.
See docs/MODEL.md and docs/RESULTS.md for assumptions, evidence, and limits.

Earlier milestones remain available at v0.4.0, v0.3.0, and the original elastic
v1.0.0 tag. This release is v1.1.0 to preserve that existing tag history. Code and
documentation remain MIT licensed. The project retains its brief AI-assistance
disclosure. This is a GitHub/source release, not a package-index publication.
