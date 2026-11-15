# V4 material-property checkpoints

The local source paper, arXiv:2604.15921v2, equations (18)-(20) and Appendix D,
Table 2, supplies the elastic and hardening fits used here. Temperatures are in
degrees Celsius; the reference temperature is 20 C (293.15 K). Stress is in MPa.
The implementation permits 0-1500 C, the appendix's plotted range; this does not
establish calibration accuracy over that entire interval.

`check_temperature_properties.py` changes only bulk and shear moduli, converting
them to Young's modulus and Poisson's ratio. All plastic parameters stay fixed.

`check_temperature_hardening.py` changes only the saturation increment and
hardening rate in the existing law:

    sigma_y(p,T) = sigma_y0 + Q(T) * (1 - exp(-b*p))
    Q(T) = 206.3425 * Phi_chi(T) MPa
    b = 1 / 0.02784959

The arctangent factor Phi_chi follows equation (19), with omega=-39.22882,
chi=-0.009042314 /K, and phi=-1.461198. The demonstration fixes sigma_y0 at
302.3605 MPa to isolate hardening. At p=0.02 the yield stress is 408.078680,
367.155028, and 347.812611 MPa for 20, 400, and 700 C, respectively.

These are equivalent yield stresses; a flowing material can carry additional
viscous overstress. They are not simple-shear stress values.

The appendix names a temperature-dependent initial-yield law but does not list
its separate fit coefficients. We have not substituted the hardening coefficients
for them or invented a curve. The cited underlying model is
[Oppermann, Denzer and Menzel (2022)](https://doi.org/10.1007/s00466-021-02103-4).
Initial-yield softening remains a separate, unresolved property checkpoint.
Neither demonstration is a full replication of the paper's thermal material model.

## Elastic softening in the restrained plate

The nonlinear thermal plate provides a controlled comparison at 16 x 16,
depth 3, alpha=1e-6, dt=0.025 s, with hot-edge heating from 20 to 320 C over
one second and a two-second hold. `--elasticity paper-reference` freezes the
paper's moduli at 20 C; `--elasticity paper-temperature` updates them at the
current interpolated physical-point temperature. Plastic history is retained.

| Elastic mode | Peak stress at 1 s [MPa] | Max p at 1 s | Peak stress at 3 s [MPa] | Max p at 3 s |
|---|---:|---:|---:|---:|
| paper-reference | 284.569082 | 0.00016523 | 250.884441 | 0.00043972 |
| paper-temperature | 281.458613 | 0.00013819 | 250.769617 | 0.00038436 |

The observed maximum plastic strain after the hold is about 12.6% lower with
temperature-dependent elasticity for this setup. Yield, hardening, viscosity,
expansion coefficient, conductivity, and mechanical supports are identical.
Fictitious stiffness uses the paper's 20 C moduli in both modes. These pointwise
maxima are not mesh-converged predictions. The original constant mode retains
the previous E=210000 MPa, nu=0.3 parameters for reproducibility.

## Hardening softening in the same plate

Keep `--elasticity paper-temperature --steps-per-second 40` active in both runs,
then compare `--hardening paper-reference` against `--hardening paper-temperature`.
Both use Q(20 C)=206.3425 MPa and b=1/0.02784959, so this is a new matched pair,
not a direct comparison with the old Q=100 MPa, b=20 runs. Initial yield is still
250 MPa and viscosity is unchanged. At fixed T the previous algorithmic tangent
continues to apply; plastic history is retained across temperature changes.

| Hardening mode | Peak stress at 1 s [MPa] | Max p at 1 s | Peak stress at 3 s [MPa] | Max p at 3 s |
|---|---:|---:|---:|---:|
| paper-reference | 281.717918 | 0.00013635 | 252.713806 | 0.00036834 |
| paper-temperature | 281.709909 | 0.00013640 | 252.658292 | 0.00036872 |

The effect is small: maximum plastic strain after the hold increases by about
0.10%. Its sample is at (0.197169, 3.947169) mm with a local temperature of
35.595 C, far below the 320 C hot reservoir. Reservoir temperature alone does
not indicate the property change at the plastic hotspot. This comparison
demonstrates the implementation without claiming strong hardening softening
in this particular restraint-driven example.

## Mechanical moment fitting with a shared thermal solution

For the same temperature-dependent elasticity/hardening comparison above,
`--method moment-fitting` compresses mechanical quadrature while leaving thermal
quadrature as quadtree. Saved nodal coordinates and temperatures were identical,
and all physical/time settings matched the quadtree reference.

Mechanical points decrease from 1628 to 1082; physical history points decrease
from 1227 to 1019. Across shared history times, the largest mean-right-ux
difference is 2.152691e-9 mm, the difference between sampled peak stress histories
is at most 0.000350191 MPa, and the maximum-plastic-strain histories differ by at
most 1.149979e-9. The final nodal displacement difference has relative Euclidean
norm 1.663814e-5. Stress/history maxima are diagnostics across different point
sets, not matched-point field norms. The fitted run needed 286339 plane-stress
updates and 281 global evaluations, with no more than three Newton corrections
per increment. Single-run timing is diagnostic and is not a speedup benchmark.
