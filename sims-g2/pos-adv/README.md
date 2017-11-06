Tests for a enhancer based positivity scheme for 2D advection
equation. All tests use polyOrder=1 and are run on periodic domain,
using a RK3 stepper.

"s" series of simulations are "vanilla" DG and "m" series are various
enhancer based DG schemes.

- s1: Gaussian initial condition; no sub-cell diffusion.
- s2: Gaussian initial condition; with sub-cell diffusion.
- s3: Cylinder initial condition; no sub-cell diffusion.
- s4: Cylinder initial condition; with sub-cell diffusion.
- s5: Square-hat initial condition; no sub-cell diffusion.
- s6: Square-hat initial condition; with sub-cell diffusion.
- s9: 5 cell sim with single non-zero cell

- m1: Gaussian initial condition; no sub-cell diffusion.
- m2: Gaussian initial condition; with sub-cell diffusion.
- m3: Cylinder initial condition; no sub-cell diffusion.
- m4: Cylinder initial condition; with sub-cell diffusion.
- m5: Square-hat initial condition; no sub-cell diffusion.
- m6: Square-hat initial condition; with sub-cell diffusion.
- m7: Same as m6, except dt increased to 1/6
- m8: IGNORE
- m9: 5 cell sim with single non-zero cell
- m10: Advection of cosine

Convergence tests

- c-s1: Gaussian initial condition; no sub-cell diffusion. 8 \times 8
- c-s2: Gaussian initial condition; no sub-cell diffusion. 16 \times 16
- c-s3: Gaussian initial condition; no sub-cell diffusion. 32 \times 32

- c-m1: Gaussian initial condition; no sub-cell diffusion. 8 \times 8
- c-m2: Gaussian initial condition; no sub-cell diffusion. 16 \times 16
- c-m3: Gaussian initial condition; no sub-cell diffusion. 32 \times 32