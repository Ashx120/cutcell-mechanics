# Development status and roadmap

V1 is an educational, bilinear immersed elasticity solver. It does not reproduce
the complete high-order thermo-viscoplastic method in the source paper.

## Implemented in V2

- Fit scaled tensor-product polynomial moments using non-negative least squares.
- Obtain reference moments from fixed-depth quadtree integration. The default
  benchmark depth is four; the physical-area demonstration uses six.
- Use the reference rule's positive-weight points as fixed candidates.
- Fall back to the quadtree rule when fitting fails or its relative residual
  exceeds `1e-10`, and emit a warning.
- Check positivity, polynomial integrals, stiffness, full-plate displacement
  and reaction agreement, and fallback behavior.
- Compare total assembly timings and cached element-integration costs.
- Report estimated break-even counts for the measured integration workload.

Candidate enrichment and automatic reference-integration error control remain
potential extensions. The current method inherits the quadtree geometry error.

## Implemented in V3

1. Verify a material-point model under loading, unloading, and different rates.
2. Add load increments, nonlinear equilibrium iterations, and trial/committed history.
3. Compare both quadrature methods using the same plate and loading history.
4. Measure reaction curves, plastic strain, iteration counts, history storage,
   and total runtime.

Completed material-point relaxation, rate, and time-step studies; analytical
and numerical tangent checks; plane-stress condensation; uniform rectangle
verification; nonlinear plate comparisons; and fictitious-stiffness sensitivity.
History counts are recorded, while actual memory bytes remain unmeasured.
See [V3 results and reproduction commands](V3_RESULTS.md).

## Further research features

1. High-order integrated-Legendre basis with p-convergence verification.
2. Hierarchical overlay refinement and gradient-jump error indicators.
3. Linear thermoelasticity with a verified stationary conduction problem.
4. Robust nonlinear globalization and automatic time-step reduction.
5. Temperature-dependent properties and calibrated material parameters.

Each feature needs its own independent verification before combined benchmarks.
Further features will follow independent verification of the mechanical model.
