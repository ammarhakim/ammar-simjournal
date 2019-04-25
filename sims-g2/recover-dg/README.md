Contains code and input files for testing recovery based discontinuous
Glaerkin/Differences based schemes.

Code directory has:

- advection.lua: Containing code to update 1D advection equation with
  periodic BCs

Advection equation tests:

- a1: p=1, N=4 (for convergence) 
- a2: p=1, N=8
- a3: p=1, N=16
- a4: p=1, N=32
- a101: p=1 N=128

- a5: p=2, N=4 (for convergence) 
- a6: p=2, N=8
- a7: p=2, N=16
- a8: p=2, N=32
- a105: p=2, N=128

- a9: p=3, N=2 (for convergence) 
- a10: p=3, N=4
- a11: p=3, N=8
- a12: p=3, N=16
- a109: p=3, N=128
