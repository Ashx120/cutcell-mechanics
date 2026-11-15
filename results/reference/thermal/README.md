# Selected thermal reference snapshots

These files are project-generated numerical data, not experimental measurements
or figures copied from the paper. They preserve the verified v4 comparison.
See provenance.json and each settings.json for the solver baseline and inputs.

- `quadtree/`: saved settings, nodal fields, and history for the reference rule.
- `fitted/`: corresponding moment-fitted case, final physical samples, and
  metrics from the selected verification run.
- `comparison.json`: common-node/history differences computed from these files.

The old flat settings format remains readable by run.py and compare_thermal_runs.py.
The selected timing in metrics.json is a single observation and may differ from
other runs recorded in the v4 report. No speedup benchmark follows from it alone.
Plots in docs/images/reference/ are generated from these data for presentation.

Recompute the comparison with:

    python scripts/compare_thermal_runs.py

Rerun either numerical experiment with:

    python run.py --config results/reference/thermal/fitted/settings.json

Change fitted to quadtree for the reference. New numerical experiments should go
under results/runs/; do not overwrite the published snapshots.
