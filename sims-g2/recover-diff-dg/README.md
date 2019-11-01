Recovery scheme for diffusion equation with cross-derivatives. Uses
recovery to compute the needed numerical fluxes, including for
cross-terms.

See code directory for the Maxima code needed to generate the kernels.

Note: the time-steps in these tests are tiny to reduce the error from
dt. In actual simulations one can set cflFrac to 1.0.

- s1: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 8x8 grid
- s2: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 16x16 grid
- s3: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 32x32 grid
- s4: kx=ky=1.0 with Dxx = Dyy = 0.0; Dxy = Dyx = 1. p=1 8x8 grid

- t1: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=2 4x4 grid
- t2: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 8x8 grid
- t3: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 16x16 grid

The following series for "bad" kernels, in which the stencil is
5-point and not 9-point

- b1: kx=ky=1.0 with Dxx = Dyy = Dxy = Dyx = 1. p=1 8x8 grid