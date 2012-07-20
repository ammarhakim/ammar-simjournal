

Simulation Index
================

Following is a list of simulation numbers with one-line
descriptions. Click on the simulation number to see the complete Lua
program for that simulation.

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
    - Wave-propagation scheme for Maxwell equation TM (8,5) mode, :math:`80\times 40` grid.
  * - :doc:`s50 <s50/s50-tm-maxwell-wave>` 
    - Same as s49 but with :math:`160\times 80` grid.
  * - :doc:`s51 <s51/s51-tm-maxwell-wave>` 
    - Same as s49 but with :math:`240\times 120` grid.
  * - :doc:`s52 <s52/s52-tm-maxwell-wave>` 
    - Same as s49 but with :math:`320\times 160` grid.
  * - :doc:`s53 <s53/s53-tm-maxwell-fdtd>` 
    - FDTD scheme for Maxwell equation TM (8,5) mode, :math:`80\times 40` grid.
  * - :doc:`s54 <s54/s54-tm-maxwell-fdtd>` 
    - Same as s53 but with :math:`160\times 80` grid.
  * - :doc:`s55 <s55/s55-tm-maxwell-fdtd>` 
    - Same as s53 but with :math:`240\times 120` grid.
  * - :doc:`s56 <s56/s56-tm-maxwell-fdtd>` 
    - Same as s53 but with :math:`320\times 160` grid.
  * - :doc:`s57 <s57/s57-pulsebox-wave>` 
    - Pulse in metal box with wave-propagation scheme.
  * - :doc:`s58 <s58/s58-pulsebox-fdtd>` 
    - Pulse in metal box with FDTD scheme.
  * - :doc:`s59 <s59/s59-pulsebox-wave>` 
    - Same as s58 but using a :math:`400 \times 400` grid.
  * - :doc:`s60 <s60/s60-pulsebox-fdtd>` 
    - Same as s58 but using a :math:`400 \times 400` grid.
  * - :doc:`s61 <s61/s61-riem-wave>` 
    - 1D Riemann problem with wave-propagation scheme.
  * - :doc:`s62 <s62/s62-riem-fdtd>` 
    - 1D Riemann problem with FDTD scheme.
  * - :doc:`s63 <s63/s63-tm-maxwell-fdtd-dual>` 
    - FDTD scheme on dual Yee-cell for Maxwell equation TM (8,5) mode,
      :math:`80\times 40` grid.
  * - :doc:`s64 <s64/s64-pulsebox-fdtd-dual>` 
    - FDTD scheme on dual Yee-cell for Maxwell equation. Pulse problem.
  * - :doc:`s65 <s65/s65-plasmabeach-maxwell>` 
    - Plasma wave beach problem, without the plasma. Using 100 cells.
  * - :doc:`s66 <s66/s66-plasmabeach-maxwell>` 
    - Plasma wave beach problem, without the plasma. Using 200 cells.
  * - :doc:`s67 <s67/s67-plasmabeach>` 
    - Plasma wave beach problem using 100 cells.
  * - :doc:`s68 <s68/s68-plasmabeach>` 
    - Plasma wave beach problem using 200 cells.
  * - :doc:`s69 <s69/s69-plasmabeach>` 
    - Plasma wave beach problem using 400 cells.
  * - :doc:`s70 <s70/s70-plasmabeach>` 
    - Plasma wave beach problem using 800 cells.
  * - :doc:`s71 <s71/s71-plasmabeach>` 
    - Plasma wave beach problem using 1600 cells.
  * - :doc:`s72 <s72/s72-cyclotron-cutoff>` 
    - Tunneling through an electron-cyclotron cutoff layer, 200 cells
  * - :doc:`s73 <s73/s73-cyclotron-cutoff>` 
    - Tunneling through an electron-cyclotron cutoff layer, 400 cells
  * - :doc:`s74 <s74/s74-icw>` 
    - Ion-cyclotron wave propagation in a 1D tokamak-like configuration, 400 cells.
  * - :doc:`s75 <s75/s75-icw>` 
    - Ion-cyclotron wave propagation in a 1D tokamak-like configuration, 200 cells.
      This simulation does not work very well as the resolution is too low.
  * - :doc:`s76 <s76/s76-icw>` 
    - Ion-cyclotron wave propagation in a 1D tokamak-like configuration, 
      low-field incidence.
  * - :doc:`s77 <s77/s77-poisson-1d>` 
    - 1D Poisson convergence test with 8 elements.
  * - :doc:`s78 <s78/s78-poisson-1d>` 
    - 1D Poisson convergence test with 16 elements.
  * - :doc:`s79 <s79/s79-poisson-1d>` 
    - 1D Poisson convergence test with 32 elements.
  * - :doc:`s80 <s80/s80-poisson-1d>` 
    - 1D Poisson convergence test with 64 elements.
  * - :doc:`s81 <s81/s81-poisson-1d>` 
    - 1D Poisson convergence test with 8 elements and Dirichlet/Neumann Bcs.
  * - :doc:`s82 <s82/s82-poisson-1d>` 
    - 1D Poisson convergence test with 16 elements and Dirichlet/Neumann Bcs.
  * - :doc:`s83 <s83/s83-poisson-1d>` 
    - 1D Poisson convergence test with 32 elements and Dirichlet/Neumann Bcs.
  * - :doc:`s84 <s84/s84-poisson-1d>` 
    - 1D Poisson convergence test with 64 elements and Dirichlet/Neumann Bcs.
  * - :doc:`s85 <s85/s85-poisson-2d>` 
    - 2D Poisson convergence test with :math:`8\times 8` elements.
  * - :doc:`s86 <s86/s86-poisson-2d>` 
    - Same as s86, except with :math:`16 \times 16` elements.
  * - :doc:`s87 <s87/s87-poisson-2d>` 
    - Same as s86, except with :math:`32 \times 32` elements.
  * - :doc:`s88 <s88/s88-poisson-2d>` 
    - Same as s86, except with :math:`64 \times 64` elements.
  * - :doc:`s89 <s89/s89-poisson-o3-1d>` 
    - 1D Poisson convergence test, 3rd order scheme, with 4 elements.
  * - :doc:`s90 <s90/s90-poisson-o3-1d>` 
    - 1D Poisson convergence test, 3rd order scheme, with 8 elements.
  * - :doc:`s91 <s91/s91-poisson-o3-1d>` 
    - 1D Poisson convergence test, 3rd order scheme, with 16 elements.
  * - :doc:`s92 <s92/s92-poisson-o3-1d>` 
    - 1D Poisson convergence test, 3rd order scheme, with 32 elements.
  * - :doc:`s93 <s93/s93-poisson-o4-1d>` 
    - 1D Poisson convergence test, 4th order scheme, with 2 elements.
  * - :doc:`s94 <s94/s94-poisson-o4-1d>` 
    - 1D Poisson convergence test, 4th order scheme, with 4 elements.
  * - :doc:`s95 <s95/s95-poisson-o4-1d>` 
    - 1D Poisson convergence test, 4th order scheme, with 8 elements.
  * - :doc:`s96 <s96/s96-poisson-o4-1d>` 
    - 1D Poisson convergence test, 4th order scheme, with 16 elements.
  * - :doc:`s97 <s97/s97-poisson-o3-2d>` 
    - 2D Poisson convergence test, 3rd order scheme, with :math:`4\times 4` elements.
  * - :doc:`s98 <s98/s98-poisson-o3-2d>` 
    - 2D Poisson convergence test, 3rd order scheme, with :math:`8\times 8` elements.
  * - :doc:`s99 <s99/s99-poisson-o3-2d>` 
    - 2D Poisson convergence test, 3rd order scheme, with :math:`16\times 16` elements.
  * - :doc:`s100 <s100/s100-poisson-o3-2d>` 
    - 2D Poisson convergence test, 3rd order scheme, with :math:`32\times 32` elements.
  * - :doc:`s101 <s101/s101-periodic-poisson-2d>` 
    - 2D Poisson convergence test with periodic BCs and :math:`32\times 32` elements
  * - :doc:`s102 <s102/s102-periodic-poisson-2d>` 
    - 2D Poisson convergence test with periodic BCs and :math:`64\times 64` elements
  * - :doc:`s103 <s103/s103-periodic-poisson-2d>` 
    - 2D Poisson convergence test with periodic BCs and :math:`128\times 128` elements
  * - :doc:`s104 <s104/s104-pb-advection-1d>` 
    - Temporal convergence study of Poisson bracket with Gaussian pulse. CFL 0.2.
  * - :doc:`s105 <s105/s105-pb-advection-1d>` 
    - Temporal convergence study of Poisson bracket with Gaussian pulse. CFL 0.1.
  * - :doc:`s106 <s106/s106-pb-advection-1d>` 
    - Temporal convergence study of Poisson bracket with Gaussian pulse. CFL 0.05.
  * - :doc:`s107 <s107/s107-pb-advection-1d>` 
    - Temporal convergence study of Poisson bracket with Gaussian pulse. CFL 0.025.
  * - :doc:`s108 <s108/s108-pb-advection-1d>` 
    - Same as s104, except with RK3 time-stepping.
  * - :doc:`s109 <s109/s109-pb-advection-1d>` 
    - Same as s108. CFL 0.1.
  * - :doc:`s110 <s110/s110-pb-advection-1d>` 
    - Same as s108. CFL 0.05.
  * - :doc:`s111 <s111/s111-pb-advection-1d>` 
    - Same as s108. CFL 0.025.
  * - :doc:`s112 <s112/s112-pb-advection-2d>` 
    - Convergence of Poisson bracket algorithm with 2nd order scheme, :math:`32\times 32` grid.
  * - :doc:`s113 <s113/s113-pb-advection-2d>` 
    - Convergence of Poisson bracket algorithm with 2nd order scheme, :math:`64\times 64` grid.
  * - :doc:`s114 <s114/s114-pb-advection-2d>` 
    - Convergence of Poisson bracket algorithm with 2nd order scheme, :math:`128\times 128` grid.
  * - :doc:`s115 <s115/s115-pb-advection-2d>` 
    - Convergence of Poisson bracket algorithm with 3rd order scheme, :math:`8\times 8` grid.
  * - :doc:`s116 <s116/s116-pb-advection-2d>` 
    - Same as s115, except with :math:`16\times 16` grid.
  * - :doc:`s117 <s117/s117-pb-advection-2d>` 
    - Same as s115, except with :math:`32\times 32` grid.
  * - :doc:`s118 <s118/s118-pb-advection-rb>` 
    - Rigid-body rotation problem for Poisson bracket, with :math:`32\times 32` grid.
  * - :doc:`s119 <s119/s119-pb-advection-rb>` 
    - Same as s118, except on :math:`64\times 64` grid.
  * - :doc:`s120 <s120/s120-pb-advection-rb>` 
    - Same as s118, except with 3-order spatial scheme.
  * - :doc:`s121 <s121/s121-pb-advection-sf>` 
    - Swirling flow problem with 3-order spatial scheme on :math:`32\times 32` grid.
  * - :doc:`s122 <s122/s122-pb-advection-2d>` 
    - Same as s115 for testing enstrophy convergence. CFL of 0.2, central flux.
  * - :doc:`s123 <s123/s123-pb-advection-2d>` 
    - Same as s122 for testing enstrophy convergence. CFL of 0.1, central flux.
  * - :doc:`s124 <s124/s124-pb-advection-2d>` 
    - Same as s122 for testing enstrophy convergence. CFL of 0.05, central flux.
  * - :doc:`s125 <s125/s125-double-shear>` 
    - Double shear problem, on :math:`64\times 64` grid, DG 2, upwind fluxes.
  * - :doc:`s126 <s126/s126-double-shear>` 
    - Same as s125, except with :math:`128\times 128` grid points.
  * - :doc:`s127 <s127/s127-double-shear>` 
    - Same as s125, except with DG3 scheme.
  * - :doc:`s128 <s128/s128-double-shear>` 
    - Same as s126, except with DG3 scheme.
  * - :doc:`s129 <s129/s129-double-shear>` 
    - Same as s125, except with CFL of 0.1.
  * - :doc:`s130 <s130/s130-double-shear>` 
    - Same as s130, except with CFL of 0.05.
  * - :doc:`s131 <s131/s131-double-shear>` 
    - Same as s125, except with a central flux.
  * - :doc:`s132 <s132/s132-double-shear>` 
    - Same as s125, except with CFL of 0.1
  * - :doc:`s133 <s133/s133-double-shear>` 
    - Same as s125, except with CFL of 0.05
  * - :doc:`s134 <s134/s134-vortex-waltz>` 
    - Vortex waltx problem, :math:`64 \times 64` grid.
  * - :doc:`s135 <s135/s135-vortex-waltz>` 
    - Vortex waltx problem, :math:`128 \times 128` grid.
  * - :doc:`s136 <s136/s136-vortex-waltz>` 
    - Vortex waltx problem, :math:`256 \times 256` grid.
  * - :doc:`s137 <s137/s137-vortex-waltz>` 
    - Vortex waltx problem, DG3, :math:`32 \times 32` grid.
  * - :doc:`s138 <s138/s138-vortex-waltz>` 
    - Vortex waltx problem, DG3, :math:`64 \times 64` grid.
  * - :doc:`s139 <s139/s139-vortex-waltz>` 
    - Vortex waltx problem, DG3, :math:`128 \times 128` grid.
  * - :doc:`s140 <s140/s140-vortex-waltz>` 
    - Vortex waltz problem, :math:`64 \times 64` grid with central flux, CFL 0.2
  * - :doc:`s141 <s141/s141-vortex-waltz>` 
    - Vortex waltz problem, :math:`64 \times 64` grid with central flux, CFL 0.1
  * - :doc:`s142 <s142/s142-vortex-waltz>` 
    - Vortex waltz problem, :math:`64 \times 64` grid with central flux, CFL 0.05
  * - :doc:`s143 <s143/s143-vlasov-free-stream>` 
    - Vlasov free-streaming operator test. DG2, upwind flux on :math:`64 \times 64` grid.
  * - :doc:`s144 <s144/s144-vlasov-free-stream>` 
    - Vlasov free-streaming operator test. DG2, upwind flux on :math:`32 \times 8` grid.
  * - :doc:`s145 <s145/s145-vlasov-free-stream>` 
    - Vlasov free-streaming operator test. DG2, upwind flux on :math:`32 \times 16` grid.
  * - :doc:`s146 <s146/s146-vlasov-free-stream>` 
    - Vlasov free-streaming operator test. DG3, upwind flux on :math:`32 \times 8` grid.
  * - :doc:`s147 <s147/s147-vlasov-free-stream>` 
    - Vlasov free-streaming operator test. DG3, upwind flux on :math:`32 \times 16` grid.
  * - :doc:`s148 <s148/s148-vlasov-fp>` 
    - Vlasov in potential well. DG2, upwind flux on :math:`32 \times 64` grid.
  * - :doc:`s149 <s149/s149-vlasov-fp>` 
    - Vlasov in potential well. DG2, upwind flux on :math:`64 \times 128` grid.
  * - :doc:`s150 <s150/s150-vlasov-fp>` 
    - Vlasov in quadratic potential well. DG2, upwind flux on :math:`64 \times 128` grid.
  * - :doc:`s151 <s151/s151-landau-damping-vp>` 
    - Vlasov-Poisson, Landau damping :math:`k=0.5`, :math:`T_e = 1.0`.
  * - :doc:`s152 <s152/s152-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 0.5`.
  * - :doc:`s153 <s153/s153-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 0.6`.
  * - :doc:`s154 <s154/s154-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 0.75`.
  * - :doc:`s155 <s155/s155-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 1.25`.
  * - :doc:`s156 <s156/s156-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 1.5`.
  * - :doc:`s157 <s157/s157-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 1.75`.
  * - :doc:`s158 <s158/s158-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 1.8`.
  * - :doc:`s159 <s159/s159-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 1.9`.
  * - :doc:`s160 <s160/s160-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 2.0`.
  * - :doc:`s161 <s161/s161-landau-damping-vp>` 
    - Same as s151, except :math:`T_e = 2.2`.
  * - :doc:`s162 <s162/s162-landau-damping-vp>` 
    - Same as s151, except :math:`\alpha=0.5`. (Nonlinear Landau damping)


