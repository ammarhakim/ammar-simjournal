Tests for an iterative solver for Poisson solver
================================================

Iterative solver for Poisson equations with Super-Time-Stepping (STS)
scheme.

STS RKL2 Tests
--------------

Each directory contains two runs, with and without extrapolation.

- s1: 32x32, RKL2 tests
- s2: 64x64, RKL2 tests
- s3: 128x128, RKL2 tests
- s4: Same as s1, p=2
- s5: Same as s2, p=2
- s6: Same as s3, p=2

STS RKL1 Tests
--------------

This is the 'm' series of tests corresponding to the 's' series of
simulations. Only difference is use of RKL1. Note that extrapolation
can't be applied each step for RKL1, for reasons I do not fully
understand.

Convergence Tests
-----------------

These are based on
http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html#convergence-of-second-order-solver-with-periodic-boundary-conditions

- c1: 8x8 grid, p=1
- c2: 16x16 grid, p=1
- c3: 32x32 grid, p=1.
- c4: 8x8 grid, p=2
- c5: 16x16 grid, p=2
- c6: 32x32 grid, p=2.

See 'd' series for RKL1 convergence.
