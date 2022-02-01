Analysis of scheme when spatial and/or temporal scales are
under-resolved.

1D cases:

- r1: Uniform flow with perturbed field; 256x256 grid. CFL=0.5,
  limiters and 1D variation in X only
- r2: Uniform flow with perturbed field; 256x256 grid.  CFL=0.5,
  limiters and 1D variation in Y only
- r3: Same as r2, cfl=0.1

2D cases:

- r4: Same as r1, except full 2D variation in perturbation. "zero" limiter (first-order)
- r5: Same as r4, except 512x512 grid
- r6: Same as r1, except CFL = 0.95
- r7: Same as r4, except CFL = 1.0
