# Repository publishing checklist

The source tree follows a scientific Python layout. The root README and
playground are the user entry points, reusable code lives under `src/`, tests
under `tests/`, and detailed scientific evidence under `docs/` and `results/v3/`.

## Before merging a development milestone

1. Install the current package and run `python -m pytest -q` and `python -m pip check`.
2. Check that the editable playground and published-data comparison work.
3. Review the diff and confirm that local environments, caches, the source-paper
   PDF, and temporary outputs are excluded. Preserve published reference snapshots.
4. Update the changelog and verify that performance claims include their measured
   scope and limitations.
5. Push the branch and wait for all configured CI jobs. The matrix covers Windows
   and Linux with Python 3.11 and 3.12, and a separate job builds/imports the wheel.
   Local verification on one machine does not establish the other CI results.
6. Merge the reviewed branch into the default branch so visitors see the current
   README, citation metadata, contribution guidance, and issue forms.

## Repository presentation

In GitHub's About panel, a suitable description is:

> Verified Q4 immersed elasticity and J2 viscoplasticity with non-negative moment-fitting quadrature.

Relevant topics include `finite-element-method`, `computational-mechanics`,
`numerical-analysis`, `viscoplasticity`, `quadrature`, and `python`.
The repository has no separate deployed website. Link to the README or results
report when sharing the project. Repository settings and topics are separate
from the tracked files and must be configured on GitHub.

## Creating a versioned release

The v3 source baseline is recorded as version `0.3.0` in the package, changelog,
and citation metadata. For later milestones, move changelog entries from
Unreleased to a dated heading and update the version and citation metadata before
tagging the reviewed commit. A GitHub release page is a separate publishing step.
Retain earlier tags, including the original `v1.0.0` milestone.

Release notes should state the new numerical behaviour, reproducible commands,
verification results, migration of the elastic playground, and known limitations.
Never infer a package-index release from a GitHub release. Publishing to a package
index would be a separate action.

## Guidance used

- [GitHub README guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes)
- [GitHub community files](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/creating-a-default-community-health-file)
- [GitHub citation files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
- [Scientific Python Development Guide](https://learn.scientific-python.org/development/)
- [Python package metadata](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)
