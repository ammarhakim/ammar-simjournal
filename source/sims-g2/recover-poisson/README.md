# Maxima code for Poisson equations with recovery appraoch

## 1D simulations

p=1 simulations:

- p1-s0: p=1, testing stencil
- p1-s1: p=1, N=4
- p1-s2: p=1, N=8
- p1-s3: p=1, N=16
- p1-s4: p=1, N=4. Cubic solution (to show it is exactly captured)
- p1-s5: p=1, N=4, periodic BCs
- p1-s6: p=1, N=8, periodic BCs

- p1-s20: p=1, testing stencil, FD approach. Taylor series
- p1-s21: p=1, N=4, FD approach (using calcRecov1CellGen for all cells)

p=2 simulations:

- p2-s0: p=2, testing stencil
- p2-s1: p=2, N=4
- p2-s2: p=2, N=8
- p2-s3: p=2, N=16
- p2-s4: p=2, N=4. k=5 solution (to show it is exactly captured)
- p2-s5: p=2, N=4. Neumann on left, Dirichlet of right BCs (sol is quartic, so exactly captured)
- p2-s6: p=2, N=4. periodic BCs
- p2-s7: p=2, N=8. periodic BCs
- p2-s8: p=2, Taylor series analysis
- p2-s9: p=2, N=4. Neumann on left, Dirichlet of right BCs (sol is quintic, so exactly captured)
- p2-s10: p=2, N=4. Neumann on left, Dirichlet of right BCs (sol is k=8)
- p2-s11: p=2, N=4. Attempt at Taylor series

- p2-s21: p=2, N=4, FD approach (using calcRecov1CellGen for all cells)

p=3 simulations:

- p3-s1: p=3, N=2
- p3-s2: p=3, N=4
- p3-s3: p=3, N=8
- p3-s4: p=3, N=16
- p3-s5: p=3, N=4. k=6 solution (to show it is exactly captured)
- p3-s6: p=3, N=4. periodic BCs
- p3-s7: p=3, N=8. periodic BCs

p=4 simulations:

- p4-s1: p=4, N=2
- p4-s2: p=4, N=4
- p4-s3: p=4, N=8
- p4-s4: p=4, N=2. k=7 solution (to show it is exactly captured)

## 2D simulations

p=1 simulations:

- p1-2d-s0: p=1, testing stencil

p=2 simulations:

- p2-2d-s0: p=2, Serendipity basis, testing stencil for bi-linear (all terms match)
- p2-2d-s1: p=2, Serendipity basis, testing stencil for bi-cubic (only bilinear terms match)
- p2-2d-s2: p=2, Tensor product basis, testing stencil for bi-cubic. All terms match
- p2-2d-s3: p=2, Tensor product basis, attempt at Taylor analysis

p=3 simulations:

- p3-2d-s1: p=2, Tensor product basis, testing stencil for 5th order source. All terms match
