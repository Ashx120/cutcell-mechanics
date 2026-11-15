# V4 results and limitations

These results describe checkpoints preserved in the `v0.4.0` source milestone,
not experimental validation. See the [walkthrough](V4_WALKTHROUGH.md) for commands and
[material notes](V4_MATERIAL_PROPERTIES.md) for equations and controlled property
comparisons. The source paper is arXiv:2604.15921v2; this project implements a
selected, simplified subset of its methods.

## Checkpoint status

| Checkpoint | Evidence and scope |
|---|---|
| Solid stationary conduction | Exact linear temperature/flux and heat balance |
| Thermal expansion | Free and restrained material-point and FEM checks |
| Spatial temperature coupling | Linear-temperature rectangle with mechanical refinement |
| Thermal plastic history | Restrained rectangle agrees with material-point history |
| Perforated conduction | Heat balance plus mesh, depth, and alpha sensitivity |
| Perforated thermoelasticity | Equilibrium, refinement, and support comparison |
| Thermal viscoplastic plate | Heating/hold history and time-step sensitivity with constant properties |
| Temperature-dependent K/G and hardening | Matched frozen-property controls; initial yield and viscosity fixed |
| Mechanical moment fitting | Shared thermal field, fewer physical history points, close nodal/history agreement |
| Initial-yield temperature fit | Unresolved coefficients; no invented curve |
| Unified configuration | Both playgrounds merged into config.py; run.py selects thermal/displacement loading |

## Perforated stationary conduction

Domain 4 x 4 mm, hole radius 1 mm, thickness 1 mm, k=0.045 W/(mm K), physical
left edge 20 C, right edge 120 C, top/bottom insulated. Depth=3, alpha=1e-6:

| Mesh | Hot-reservoir input [W] |
|---|---:|
| 16 x 16 | 4.080917191 |
| 32 x 32 | 4.079494971 |
| 64 x 64 | 4.079043968 |

Successive changes are about 0.0349% and 0.0111%. Cold-edge removal matches hot
input to numerical precision. On 16 x 16, alpha=1e-8 gives 4.080917014 W. On
32 x 32, depth=4 gives 4.079432027 W. This supports low sensitivity of total heat
transfer in these tests; it does not establish local heat-flux convergence.

## Elastic mechanical response

Use the same thermal case, E=210000 MPa, nu=0.3, alpha_thermal=12e-6 /K,
reference temperature 20 C, and left-x/bottom-y rollers:

| Mesh | Mean right ux [mm] | Sampled peak von Mises [MPa] |
|---|---:|---:|
| 16 x 16 | 0.002552529 | 103.570280 |
| 32 x 32 | 0.002559813 | 109.699101 |
| 64 x 64 | 0.002563755 | 113.672490 |

At 32 x 32, depth=4 gives 0.002559829 mm and 109.695274 MPa. Integration-depth
sensitivity is smaller than mesh sensitivity for these quantities. Displacement
stabilises more quickly than the sampled stress peak; the latter is not declared
converged. Mechanical residuals and net reaction forces are near numerical zero.

## Support-driven hotspot

At 32 x 32 and hot-edge temperature 320 C, the original rollers give a sampled
elastic peak of 329.097303 MPa near (0.098584, 3.973584) mm. Replacing edge rollers
with three constraints removing only rigid motion gives 26.195310 MPa near
(0.496698, 0.878302) mm. The temperature field is unchanged and the minimal-support
reactions are near zero.

Thus the demonstrated plastic hotspot near the top-left corner is strongly
restraint driven. It is not evidence of plasticity concentrated at the hole in a
freely expanding body. The support alternatives represent different physical
problems. A possible corner singularity has not been established or ruled out;
global equilibrium alone does not validate a local stress peak.

## Constant-property thermal viscoplastic time study

Mesh 16 x 16, depth=3, alpha=1e-6, original constant properties. Hot edge ramps
20 -> 320 C over 1 s; cold edge remains at 20 C; temperatures hold to 3 s.

| dt [s] | Peak stress at 1 s [MPa] | Max p at 1 s | Peak stress at 3 s [MPa] | Max p at 3 s |
|---|---:|---:|---:|---:|
| 0.05 | 284.389393 | 0.00018538 | 250.929009 | 0.00046017 |
| 0.025 | 285.592429 | 0.00017634 | 250.924476 | 0.00045973 |
| 0.0125 | 286.225165 | 0.00017145 | 250.922255 | 0.00045939 |

The shrinking endpoint differences are consistent with first-order temporal
convergence. Plastic strain during heating remains more sensitive than the
end-of-hold response. dt=0.025 s is a development compromise; its 1 s maximum
plastic strain differs from the finer result by about 2.8%. This study does not
establish time convergence for every property mode or local field, and nonlinear
mesh convergence remains outstanding.

## Mechanical quadrature under temperature-dependent properties

Both runs use temperature-dependent elasticity/hardening, initial yield=250 MPa,
fixed viscosity/expansion, mesh 16 x 16, depth=3, alpha=1e-6, dt=0.025 s.
The stationary heat problem always uses quadtree; nodal temperature fields match.

| Mechanical rule | Total points | Physical history points | Peak stress at 3 s [MPa] | Max p at 3 s |
|---|---:|---:|---:|---:|
| Quadtree | 1628 | 1227 | 252.658292 | 0.00036872 |
| Moment fitting | 1082 | 1019 | 252.658301 | 0.00036872 |

Moment fitting uses 33.5% fewer total points and 17.0% fewer physical history
points. The largest mean-right-ux history difference is 2.152691e-9 mm. The final
nodal displacement difference has relative Euclidean norm 1.663814e-5 across
active nodes. Sampled-maximum stress histories differ by at most 0.000350191 MPa;
this is not a matched-point stress-field norm.

The recorded user runs took 98.993 s for quadtree and 80.794 s for fitting, an
18.4% shorter individual solve. These are single-run observations, not a repeated
performance benchmark. The fitted run used 286339 plane-stress updates, 281 global
evaluations, and at most three Newton corrections per increment.

## Limits of the current model

- Heat conduction is stationary with constant conductivity. Heating histories
  prescribe scaled equilibrium fields; no thermal inertia or transient equation.
- Coupling is one-way: no plastic heating, mechanical feedback into temperature,
  or full thermodynamic energy-balance implementation.
- Thermal expansion is constant. Initial yield and viscosity remain fixed even
  when elastic/hardening fits are active. Missing yield-fit coefficients remain
  unresolved. The model is not the complete paper calibration.
- The hole is approximated by cut integration and weak fictitious conductivity
  and stiffness; its boundary conditions are approached as alpha becomes small.
- Q4, Cartesian, small-strain plane stress only. No hp adaptivity, remeshing,
  moving integration rules, or history transfer between rules.
- Plain global Newton has no line search or adaptive time stepping. A large
  single temperature jump can fail without committing history.
- Nonlinear mesh convergence, independent local-field validation, and repeated
  timing benchmarks remain open. Agreement between quadrature rules is not
  proof of agreement with an exact continuum solution or experiment.

All 89 automated tests pass at this checkpoint. Generated walkthrough results
remain ignored by Git; these numerical summaries document the recorded runs.
Broader beginner-facing UX and repository cleanup follow the v4 source milestone
as a separate final-release stage.
