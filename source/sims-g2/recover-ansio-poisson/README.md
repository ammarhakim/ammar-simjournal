Test cases for debugging/understanding behaviour of the RDG solver for
anisotropic Poisson equations.

# Constant, inclined field

- bconst-1: 4x4, 15 deg field inclination, bi-cubic. Solution is exact
- bconst-2: 4x4, 15 deg field inclination, fifth order.
- bconst-3: 8x8, 15 deg field inclination, fifth order.
- bconst-4: 16x16, 15 deg field inclination, fifth order.
- bconst-5: 32x32, 15 deg field inclination, fifth order.
- bconst-6: 16x16, 15 deg field inclination, Gaussian blob src

# O-point simulations

- opoint-1: 4x4 bi-quadratic
- opoint-2: 8x8 bi-quadratic
- opoint-3: 16x16 bi-quadratic
- opoint-4: 16x16 Gaussian blob src
- opoint-5: 17x17 Gaussian blob src

# O-point outside domain simulations

- out-opoint-1: 4x4 bi-cubic
- out-opoint-2: 8x8 bi-cubic
- out-opoint-3: 16x16 bi-cubic

- out-opoint-5: 8x8 bi-cubic

# X-point simulations


- xpoint-1: 4x4 bi-quadratic
- xpoint-2: 8x8 bi-quadratic
- xpoint-3: 16x16 bi-quadratic

- xpoint-4: 4x4 bi-cubic
- xpoint-5: 8x8 bi-cubic
- xpoint-6: 16x16 bi-cubic

- xpoint-7: 16x16 Gaussian blob src

# Projection of O-point Dij

Just projecting Dij to check if recovery is correct.

- Dij-proj-1: 4x4 grid; zet=1e9
- Dij-proj-2: 16x16 grid; zet=1e1
- Dij-proj-3: 16x16 grid; zet=1e9
