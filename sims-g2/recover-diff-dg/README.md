Recovery scheme for diffusion equation with cross-derivatives. Uses
recovery to compute the needed numerical fluxes, including for
cross-terms.

See code directory for the Maxima code needed to generate the kernels.

Note: the time-steps in these tests are tiny to reduce the error from
dt. In actual simulations one can set cflFrac to 1.0.

- s1: kx=ky=1.0 with Dxx = Dyy = 1.0; Dxy = Dyx = 0.9 p=1 8x8 grid
- s2: Same as s1, 16x16 grid
- s3: Same as s1, 32x32 grid
- s4: Dxx = Dyy = 0. Dxy = Dyx = 1. This shows eventually instability
  as D is not positive definite.
- s5: ky = -1 test. This mode damps very slowly compared to kx=ky=1 mode
- s6: Decay of Gaussian with Dxx = Dyy = 1.0; Dxy = Dyx = 0.9 p=1 12x12 grid

- t1: Same as s1, except p=2, 4x4 grid
- t2: Same as t1, except p=2, 8x8 grid
- t3: Same as t1, except p=2, 16x16 grid
- t4: Same as t2, run for longer with cflFrac = 0.9

- w1: with Dxx = Dyy = 1.0; Dxy = Dyx = 0.9 p=1 8x8 grid; source term
- w2: Anisotropic diffusion problem. p=1, 4x4 grid
- w3: Anisotropic diffusion problem. p=1, 8x8 grid
- w4: Anisotropic diffusion problem. p=1, 16x16 grid
- w5: Anisotropic diffusion problem. p=1, 32x32 grid

- x2: Same as w2, except p=2

The following series for "bad" kernels, in which the stencil is
5-point and not 9-point

- b1: kx=ky=1.0 with Dxx = Dyy = 1.0. Dxy = Dyx = 1/6. p=1 8x8 grid
- b2: kx=ky=1.0 with Dxx = Dyy = 1.0. Dxy = Dyx = 0.9 p=1 8x8 grid

Tests for diffusion terms in FPO. g = (Dxx*x^2+Dyy*y^2)/2 + Dxy*x*y to
mimic the constant diffusion case.

- f1: kx=ky=1.0 with Dxx = Dyy = 1.0; Dxy = Dyx = 0.9 p=1 8x8 grid
