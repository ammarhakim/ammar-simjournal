Input files for benchmarking the FOP (LBO version). Used in the FPO
paper.

Relaxation tests (Problem 1)

- r1: Relaxation of step function to Maxwellian. p=1 case
- r2: Relaxation of step function to Maxwellian. p=2 case
- r3: Same as r1, except CFL factor of 2x smaller
- r4: Same as r1, except 2x more cells in V (entropy convergence comparison)

Sod-shock tests (Problem 2)

- s1: Sod-shock problem. MFP = 0.1*Lx
- s2: Sod-shock problem. MFP = 0.01*Lx
- s3: Sod-shock problem. MFP = 0.001*Lx
- s4: Exact solution to Euler Sod-shock
