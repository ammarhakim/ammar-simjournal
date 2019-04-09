A series of tests to look at the peculiar behavior of our FV Maxwell
solver.

1D periodic domain, propagation of an EM wave:

- s1: 32 cells, cfl = 0.25
- s2: 32 cells, cfl = 0.5
- s3: 32 cells, cfl = 0.75
- s4: 48 cells, cfl = 0.75

- d1: same as s1, except 8 cells DG scheme; p=1
- d2: same as d1, except 16 cells DG scheme; p=1
- d3: same as d2, except p=2.
- d4: same as d2 (p=1, 16 cells) but central fluxes

- f1: Same as s1, except FDTD method (G1 input file)
- f2: Same as s2, except FDTD method (G1 input file)
- f3: Same as s3, except FDTD method (G1 input file)
- f3: Same as s4, except FDTD method (G1 input file)