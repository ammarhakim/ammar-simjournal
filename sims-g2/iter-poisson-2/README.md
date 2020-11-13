Tests for an iterative solvers for Poisson equations
====================================================

Iterative solvers for Poisson equations. Uses the in-built Gkeyll
updater. Most (not all) of these driectories contain three Lua files:
one for RKL1, one for RKL1-extrapolated and one for Richarson
Second-order method.

Convergence tests 1D
--------------------

- a1: 8 grid, p=1
- a2: 16 grid, p=1
- a3: 32 grid, p=1
- a4: 64 grid, p=1

- b1: 8 grid, p=2
- b2: 16 grid, p=2
- b3: 32 grid, p=2
- b4: 64 grid, p=2

Convergence tests 2D
--------------------

http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html#convergence-of-second-order-solver-with-periodic-boundary-conditions

- c1: 8x8 grid, p=1
- c2: 16x16 grid, p=1
- c3: 32x32 grid, p=1
- c4: 64x64 grid, p=1
- c5: 128x128 grid, p=1
- c6: 256x256 grid, p=1

- d1: 8x8 grid, p=2
- d2: 16x16 grid, p=2
- d3: 32x32 grid, p=2
- d4: 64x64 grid, p=2
- d5: 128x128 grid, p=2

Convergence tests 3D
--------------------

- e1: 8x8x8 grid, p=1
- e2: 16x16x16 grid, p=1
- e3: 32x32x32 grid, p=1
- e4: 64x64x64 grid, p=1

- f1: 8x8x8 grid, p=2
- f2: 16x16x16 grid, p=2

Gaussian source in 2D
---------------------

- m1: 16x16 grid, p=1
- m2: 32x32 grid, p=1
- m3: 64x64 grid, p=1
- m4: 128x128 grid, p=1
- m5: 128x1 grid, p=1
- m6: 128x128 grid, p=1, parallel
- m7: 64x64x64 grid, p=1, parallel 2 cores
- m8: 64x64x64 grid, p=1, parallel 4 cores

- n1: 16x16 grid, p=2
- n2: 32x32 grid, p=2
- n3: 64x64 grid, p=2
- n4: 128x128 grid, p=2
