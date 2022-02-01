Benchmark problems for Vlasov-Maxwell solver.

Free Streaming
--------------

- f1: Free-streaming, 1x1v, Serendipity, 32x16. p=1
- f2: Free-streaming, 1x1v, Serendipity, 32x16. p=2

Free Streaming on open domain
-----------------------------

- o1: Free-streaming, 1x1v, Serendipity, 32x16. p=1
- o1: Free-streaming, 1x1v, Serendipity, 32x16. p=2

Stand-alone Maxwell
-------------------

- m1: Plane-wave in 1D. Serendipity p=1. NX=16
- m2: Plane-wave in 1D. Serendipity p=2. NX=8
- m3: Plane-wave in 2D. Serendipity p=1. 16x16 grid.
- m4: Plane-wave in 2D. Max-order p=1. 16x16 grid.
- m5: Plane-wave in 2D. Serendipity p=2. 8x8 grid.
- m6: Plane-wave in 2D. Max-order p=2. 8x8 grid.
- m7: Plane-wave in 3D. Serendipity p=1. 16x16x16 grid.
- m8: Plane-wave in 3D. Max-order p=1. 16x16x16 grid.
- m9: Plane-wave in 3D. Serendipity p=2. 8x8x8 grid.
- m10: Plane-wave in 3D. Max-order p=2. 8x8x8 grid.
- m11: EM pulse in a box. Serendipity p=1, 32x32 grid

Potential Well
--------------

- p1: Potential well, Serendipity 1x1v, p=1
- p2: Potential well, Serendipity 1x1v, p=2
- p3: Potential well, Maximal-order 1x1v, p=1

Two-stream instability
----------------------

- t1: Two-stream instability. Serendipity 1x1v, 64x16 p=2
- t2: Two-stream instability. Serendipity 1x1v, 64x32 p=2
- t3: Two-stream instability. Serendipity 1x1v, 32x24 p=3
- t4: Two-stream instability. Max-order 1x1v, 32x24 p=3

Electrostatic Sod-shock problem
-------------------------------

- esss1: Electrostatic Sod-shock. mi=25, Te/T1=1

Advection in a constant magnetic field
--------------------------------------

- c1: Constant magnetic field, Serendipity 1x2v, 4x16^2 p=1
- c2: Same as c1, except p=2
- c3: Constant magnetic field, Max-order 1x2v, 4x16^2 p=1
- c4: Same as c3, except p=2

Weibel instability
------------------

- w1: Serendipity, 64x16^2, p=1
- w2: Max-order, 64x16^2, p=1
- w3: Serendipity, 64x16^2, p=2
- w4: Max-order, 64x16^2, p=2

Orsag-Tang Problem
------------------

- ot1: 64^2x12^3, p=2