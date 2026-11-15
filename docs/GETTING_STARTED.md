# Your first simulation

CutCell Mechanics asks: **how does a plate with a hole deform and develop stress
when one side is heated, and can we do the numerical integration with fewer points?**
You can run the example without knowing the finite element method (FEM). The
technical notes explain the equations when you want to look deeper.

## Install and run

Install Python 3.11 or newer and Git. Open a terminal in the project directory.
On Windows PowerShell:

```powershell
git clone https://github.com/Ashx120/cutcell-mechanics.git
cd cutcell-mechanics
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe run.py
```

On Linux/macOS, use `.venv/bin/python` instead of `.\.venv\Scripts\python.exe`.
No environment activation is required. If you already installed the project,
start with the last command. The coarse default took about 12 seconds on the
development machine; speed varies with hardware and settings.

The terminal prints progress and the absolute output-folder path. Every normal
run creates a separate timestamped folder under `results/runs/`; previous results
are preserved. The first example is for learning, not a converged engineering
prediction. Finer meshes and smaller increments take longer.

## What happens in this example?

The plate is 4 x 4 mm, with a quarter-circle hole at the lower-left and thickness
1 mm. Initially the temperature is 20 C everywhere. Over one second, the right
edge rises to 320 C while the left stays at 20 C. Those temperatures then hold
for two seconds.

Rollers prevent horizontal motion on the left and vertical motion on the bottom.
The hot material wants to expand, but neighbouring regions and supports resist
that deformation. If stress exceeds the material's current yield level, some
plastic strain accumulates. During the hold, that time-dependent plastic flow
can reduce the stress. This is called **stress relaxation**.

The temperature fields are stationary equilibrium snapshots. The time axis
describes the imposed schedule and mechanical response; it does not simulate
the time required for heat to diffuse through the plate.

## Open these files first

| File | What it tells you |
|---|---|
| `summary.txt` | Run settings, headline results, and interpretation reminders |
| `history.png` | Heating schedule, peak stress, and accumulated plastic strain over time |
| `fields.png` | Final temperature, stress, and plastic-strain locations |
| `settings.json` | The complete resolved inputs, including command-line overrides |
| `metrics.json` | Integration-point counts, equilibrium residuals, iterations, and measured solve time |
| `history.csv`, `nodes.csv`, `physical_samples.csv` | Numerical data for your own plots or comparisons |

In the plots, each dot is a physical **integration point**, where the code
evaluates material behaviour. Denser dots near the hole come from integration
subdivision; they are not necessarily a finer displacement mesh. Plot colour
scales adjust separately, so compare numbers as well as colours.

**Stress** measures internal force per area, here in MPa. **Von Mises stress**
combines the stress components into a scalar used by the yielding model.
**Accumulated plastic strain** measures irreversible deformation history; it is
dimensionless and does not decrease when loading reverses. It is not the same
as displacement or total strain.

The hotspot may change with resolution. In the finer recorded example it is
near the cold supported corner, not the hottest edge. Supports and compatibility
matter as well as temperature. Read [the model scope](MODEL.md) before treating
a coloured hotspot as a reliable physical prediction.

## Change one input at a time

Open `config.py`, edit a setting, save, and run `run.py` again.

1. Change `right_temperature` from 320 to 120. Expect less thermal loading and
   potentially no plastic accumulation.
2. Return to 320 and change `SUPPORTS` to `"minimal"`. This releases edge
   restraints while preventing whole-body motion. It changes the physical setup;
   stresses are much lower in the recorded fine-mesh comparison.
3. Increase `nx` and `ny` from 8 to 16. This refines the displacement mesh.
   Small residuals alone do not show that the mesh is sufficiently fine.
4. Change `QUADRATURE["method"]` between `"quadtree"` and `"moment-fitting"`
   with all other settings identical. Compare common-node displacements as well
   as point counts; the integration-point locations differ.

The starter uses constant properties so changes to MATERIAL take effect directly.
Units are mm, N, MPa, seconds, and degrees C. Conductivity is W/(mm K), so
0.045 represents 45 W/(m K). The `paper-*` property modes replace the relevant
MATERIAL elastic/hardening constants with the documented literature fits;
initial yield and viscosity remain supplied by MATERIAL. See [MODEL.md](MODEL.md).

To reproduce a saved configuration:

```powershell
.\.venv\Scripts\python.exe run.py --config results/runs/YOUR-RUN/settings.json
```

Command-line overrides remain available, for example:

```powershell
.\.venv\Scripts\python.exe run.py --mesh 16 --steps-per-second 40
.\.venv\Scripts\python.exe run.py --help
```

An explicit `--output PATH` uses that exact folder. If it is nonempty, choose a
new folder or add `--overwrite` deliberately. No automatic deletion is performed.

## If something goes wrong

- **Import error:** use the environment's Python and repeat the installation command.
- **Configuration error:** check the named setting and units; unknown keys are rejected.
- **Newton equilibrium fails:** reduce the heating increment by increasing
  `steps_per_second`, or reduce the loading. Failed trial history is not committed.
- **Runs take longer:** use the coarse default to learn, and reserve the recorded
  16 x 16 / 40-steps-per-second case for detailed comparisons.
- **Unexpected peak stress:** check supports, mesh/time sensitivity, and whether
  you are comparing different integration-point samples or colour scales.

For recorded evidence, see [RESULTS.md](RESULTS.md). For equations and individual
verification examples, see [DEVELOPMENT.md](DEVELOPMENT.md).
