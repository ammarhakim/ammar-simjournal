

Simulation Index
================

Following is a list of simulation numbers with one-line descriptions.

.. list-table::
  :header-rows: 1
  :widths: 10,90

  * - Number
    - Description
  * - :doc:`s1 <s1/s1-periodic-poisson>` 
    - Poisson solve on a 2D periodic grid, with Gaussian source
  * - :doc:`s2 <s2/s2-periodic-poisson>` 
    - Poisson solve on a 2D periodic grid, with anisotropic Gaussian source
  * - :doc:`s3 <s3/s3-periodic-poisson>` 
    - Poisson solve on a 2D periodic grid, with source a sum of two Gaussians
  * - :doc:`s4 <s4/s4-periodic-poisson>` 
    - Same as s3, but with different grid spacing in X and Y.
  * - :doc:`s5 <s5/s5-euler-shock-wave>` 
    - 1D Euler Sod-shock with sonic point in rarefaction. Wave-propagation algorithm.
  * - :doc:`s6 <>` 
    - Exact solution to s5.
  * - :doc:`s7 <s7/s7-euler-shock-muscl>` 
    - Same as s5, except using the MUSCL-Hancock scheme.
  * - :doc:`s8 <s8/s8-euler-shock-wave>` 
    - 1D Euler shock with low density/pressure region. Wave-propagation scheme.
  * - :doc:`s9 <>` 
    - Exact solution to s8.
  * - :doc:`s10 <s10/s10-euler-shock-muscl>` 
    - Same as s8, except using the MUSCL-Hancock scheme.
  * - :doc:`s11 <s11/s11-euler-shock-wave>` 
    - 1D Noh problem. Wave-propagation algorithm.
  * - :doc:`s12 <>` 
    - Exact solution to s11.
  * - :doc:`s13 <s13/s13-euler-shock-muscl>` 
    - Same as s11, except using the MUSCL-Hancock scheme.
  * - :doc:`s14 <s14/s14-euler-shock-wave>` 
    - 1D Euler shock with a stationary contact discontinuity. Wave-propagation scheme.
  * - :doc:`s15 <>` 
    - Exact solution to s14
  * - :doc:`s16 <s16/s16-euler-shock-muscl>` 
    - Same as s14, except using the MUSCL-Hancock scheme.
  * - :doc:`s17 <s17/s17-euler-shock-wave>` 
    - 1D Euler shock with two strong shocks. Wave-propagation scheme.
  * - :doc:`s18 <>` 
    - Exact solution to s17
  * - :doc:`s19 <s19/s19-euler-shock-muscl>` 
    - Same as s17, except using the MUSCL-Hancock scheme.
  * - :doc:`s20 <s20/s20-euler-shock-wave>` 
    - 1D Euler with stationary contact discontinuity. Wave-propagation scheme.
  * - :doc:`s21 <>` 
    - Exact solution to s20
  * - :doc:`s22 <s22/s22-euler-shock-muscl>` 
    - Same as s20, except using the MUSCL-Hancock scheme.
  * - :doc:`s23 <s23/s23-euler-shock-wave>` 
    - 1D Euler with slowly moving contact discontinuity. Wave-propagation scheme.
  * - :doc:`s24 <>` 
    - Exact solution to s23
  * - :doc:`s25 <s25/s25-euler-shock-muscl>` 
    - Same as s23, except using the MUSCL-Hancock scheme.
  * - :doc:`s26 <s26/s26-euler-shock-wave>` 
    - 1D Euler with sharp spike in density. Wave-propagation scheme.
  * - :doc:`s27 <>` 
    - Exact solution to s26
  * - :doc:`s28 <s28/s28-euler-shock-muscl>` 
    - Same as s26, except using the MUSCL-Hancock scheme.
  * - :doc:`s29 <s29/s29-euler-blastwave-wave>` 
    - 1D Euler Woodward-Collela blast wave problem. Wave-propagation scheme.
  * - :doc:`s30 <s30/s30-euler-blastwave-wave>` 
    - Same as s29 run with higher-resolution to serve as an "exact" solution.
  * - :doc:`s31 <s31/s31-euler-blastwave-muscl>` 
    - Same as s29, except using the MUSCL-Hancock scheme.
  * - :doc:`s32 <s32/s32-rte-slab>` 
    - Slab RTE with Mie scattering with :math:`L=8`.
  * - :doc:`s33 <s33/s33-rte-slab>` 
    - Slab RTE with Haze-L phase function with :math:`L=82`.
  * - :doc:`s34 <s34/s34-rte-slab>` 
    - Same as s33 but with :math:`\mu_0=0.5`.
  * - :doc:`s35 <s35/s35-rte-slab>` 
    - Same as s33 but with :math:`\varpi=1.0`.
  * - :doc:`s36 <s36/s36-twofluid-shock>` 
    - Two-fluid shock problem with :math:`q_i/m_i = 1.0`.
  * - :doc:`s37 <s37/s37-twofluid-shock>` 
    - Two-fluid shock problem with :math:`q_i/m_i = 10.0`.
  * - :doc:`s38 <s38/s38-twofluid-shock>` 
    - Two-fluid shock problem with :math:`q_i/m_i = 100.0`
  * - :doc:`s39 <s39/s39-twofluid-shock>` 
    - Two-fluid shock problem with :math:`q_i/m_i = 1000.0`
  * - :doc:`s40 <s40/s40-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 10` and 100 cells.
  * - :doc:`s41 <s41/s41-sqpulse-exact>` 
    - Exact solution of dispersive Euler equations.
  * - :doc:`s42 <s42/s42-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 10` and 200 cells.
  * - :doc:`s43 <s43/s43-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 10` and 300 cells.
  * - :doc:`s44 <s44/s44-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 10` and 400 cells.
  * - :doc:`s45 <s45/s45-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 100` and 200 cells.
  * - :doc:`s46 <s46/s46-dispersive-euler>` 
    - Dispersive Euler equations with :math:`\omega_c = 100` and 400 cells.
  * - :doc:`s47 <s47/s47-dispersive-euler>` 
    - Sod-shock for dispersive Euler equations.
  * - :doc:`s48 <s48/s48-dispersive-euler>` 
    - Same as s45 but larger time-step.
  * - :doc:`s49 <s49/s49-tm-maxwell-wave>` 
    - Wave-propagation scheme for Maxwell equation TM (8,5) mode. :math:`80\times 40` grid.
