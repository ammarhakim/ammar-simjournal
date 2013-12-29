


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
  * - :doc:`s163 <s163/s163-landau-damping-vp>` 
    - Vlasov-Poisson, Landau damping :math:`k=\sqrt{0.08}`, :math:`T_e = 1.0` on a large domain.
  * - :doc:`s164 <s164/s164-landau-damping-vp>` 
    - Same as s163, except :math:`k=\sqrt{0.1}`.
  * - :doc:`s165 <s165/s165-landau-damping-vp>` 
    - Same as s163, except :math:`k=\sqrt{0.12}`.
  * - :doc:`s166 <s166/s166-landau-damping-vp>` 
    - Same as s163, except :math:`k=\sqrt{0.06}`.
  * - :doc:`s167 <s167/s167-landau-damping-vp>` 
    - Same as s163, except :math:`k=0.25`.
  * - :doc:`s168 <s168/s168-landau-damping-vp>` 
    - Same as s163, except :math:`k=0.5`.
  * - :doc:`s169 <s169/s169-landau-damping-vp>` 
    - Same as s163, except :math:`k=0.75`.
  * - :doc:`s170 <s170/s170-landau-damping-vp>` 
    - Same as s163, except :math:`k=1.0`.
  * - :doc:`s171 <s171/s171-landau-damping-vp>` 
    - Same as s163, except :math:`k=2.5`.
  * - :doc:`s172 <s172/s172-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 0.1`.
  * - :doc:`s173 <s173/s173-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 0.5`.
  * - :doc:`s174 <s174/s174-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 0.75`.
  * - :doc:`s175 <s175/s175-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 1.0`.
  * - :doc:`s176 <s176/s176-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 1.5`.
  * - :doc:`s177 <s177/s177-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 2.0`.
  * - :doc:`s178 <s178/s178-ion-sound>` 
    - Ion-acoustic waves with quasi-neutrality condition and :math:`T_i/T_e = 0.3`.
  * - :doc:`s179 <s179/s179-landau-damping-vp>` 
    - Vlasov-Poisson simulation to test momentum conservation, :math:`8\times 128` cells.
  * - :doc:`s180 <s180/s180-landau-damping-vp>` 
    - Same as s179, except :math:`8\times 32` cells.
  * - :doc:`s181 <s181/s181-landau-damping-vp>` 
    - Same as s179, except :math:`16\times 128` cells.
  * - :doc:`s182 <s182/s182-landau-damping-vp>` 
    - Same as s179, except :math:`32\times 128` cells.
  * - :doc:`s183 <s183/s183-landau-damping-vp>` 
    - Same as s179 except, :math:`64\times 128` cells.
  * - :doc:`s184 <s184/s184-landau-damping-vp>` 
    - Vlasov-Poisson simulation to test energy conservation, :math:`16\times 32` cells, CFL 0.3.
  * - :doc:`s185 <s185/s185-landau-damping-vp>` 
    - Same as s184, except CFL :math:`0.3/2`
  * - :doc:`s186 <s186/s186-landau-damping-vp>` 
    - Same as s184, except CFL :math:`0.3/4`
  * - :doc:`s187 <s187/s187-landau-damping-vp>` 
    - Same as s184, except CFL :math:`0.3/8`
  * - :doc:`s188 <s188/s188-landau-damping-vp>` 
    - Same as s179, except :math:`16\times 32` cells.
  * - :doc:`s189 <s189/s189-landau-damping-vp>` 
    - Same as s179, except :math:`16\times 64` cells.
  * - :doc:`s190 <s190/s190-landau-damping-vp>` 
    - Same as s179, except :math:`16\times 128` cells.
  * - :doc:`s191 <s191/s191-landau-damping-vp>` 
    - Same as s179 except, :math:`128\times 128` cells.
  * - :doc:`s192 <s192/s192-landau-damping-vp>` 
    - Vlasov-Poisson simulation to test energy conservation, :math:`16\times 32` cells, CFL 0.2, DG3.
  * - :doc:`s193 <s193/s193-landau-damping-vp>` 
    - Same as s192, except CFL :math:`0.1`.
  * - :doc:`s194 <s194/s194-landau-damping-vp>` 
    - Same as s192, except CFL :math:`0.05`.
  * - :doc:`s195 <s195/s195-landau-damping-vp>` 
    - Same as s192, except CFL :math:`0.025`.
  * - :doc:`s196 <s196/s196-landau-damping-vp>` 
    - Vlasov-Poisson simulation to test momentum conservation, :math:`8\times 128` cells, DG3
  * - :doc:`s197 <s197/s197-landau-damping-vp>` 
    - Same as s196, except :math:`16\times 128` cells.
  * - :doc:`s198 <s198/s198-landau-damping-vp>` 
    - Same as s196, except :math:`32\times 128` cells.
  * - :doc:`s199 <s199/s199-landau-damping-vp>` 
    - Same as s196, except :math:`64\times 128` cells.
  * - :doc:`s200 <s200/s200-landau-damping-vp>` 
    - Same as s196, except :math:`128\times 128` cells.
  * - :doc:`s201 <s201/s201-aux-dg-advection-rb>` 
    - Rigid-body flow with nodal DG updater, :math:`16\times 16` cells.
  * - :doc:`s202 <s202/s202-aux-dg-advection-rb>` 
    - Rigid-body flow with nodal DG updater, :math:`32\times 32` cells.
  * - :doc:`s203 <s203/s203-aux-dg-advection-rb>` 
    - Rigid-body flow with nodal DG updater, :math:`64\times 64` cells.
  * - :doc:`s204 <s204/s204-aux-dg-advection-rb>` 
    - Rigid-body flow with nodal DG updater polyOrder 2, :math:`32\times 32` cells.
  * - :doc:`s205 <s205/s205-aux-dg-advection-swirl>` 
    - Swirling flow with nodal DG updater polyOrder 2, :math:`32\times 32` cells.
  * - :doc:`s206 <s206/s206-aux-dg-advection-helix-3d>` 
    - Helical flow with nodal DG updater polyOrder 2, :math:`16\times 16\times 16` cells.
  * - :doc:`s207 <s207/s207-advect-diffuse>` 
    - Advection-diffusion test, 8 cells
  * - :doc:`s208 <s208/s208-advect-diffuse>` 
    - Advection-diffusion test, 16 cells
  * - :doc:`s209 <s209/s209-advect-diffuse>` 
    - Advection-diffusion test, 32 cells
  * - :doc:`s210 <s210/s210-advect-diffuse>` 
    - Advection-diffusion test, 64 cells
  * - :doc:`s211 <s211/s211-advect-diffuse>` 
    - Advection-diffusion test, 8 cells, 3-point flux
  * - :doc:`s212 <s212/s212-advect-diffuse>` 
    - Advection-diffusion test, 16 cells, 3-point flux
  * - :doc:`s213 <s213/s213-advect-diffuse>` 
    - Advection-diffusion test, 32 cells, 3-point flux
  * - :doc:`s214 <s214/s214-advect-diffuse>` 
    - Advection-diffusion test, 64 cells, 3-point flux
  * - :doc:`s215 <s215/s215-hw>` 
    - Hasegawa-Wakatani, with adiabacity parameter 0.1
  * - :doc:`s216 <s216/s216-hw>` 
    - Hasegawa-Wakatani, with adiabacity parameter 0.01
  * - :doc:`s217 <s217/s217-hw>` 
    - Hasegawa-Wakatani, with adiabacity parameter 0.3
  * - :doc:`s218 <s218/s218-hw>` 
    - Hasegawa-Wakatani, with adiabacity parameter 1.0
  * - :doc:`s219 <s219/s219-hw>` 
    - Hasegawa-Wakatani, with adiabacity parameter 2.0
  * - :doc:`s220 <s220/s220-euler-shock-wave>` 
    - 1D Euler shock with low density/pressure
      region. Wave-propagation scheme with pressure/density fix. (See s8)
  * - :doc:`s221 <s221/s221-mhw>` 
    - Modified Hasegawa-Wakatani, with adiabacity parameter 0.5
  * - :doc:`s222 <s222/s222-mhw>` 
    - Modified Hasegawa-Wakatani, with adiabacity parameter 1.0
  * - :doc:`s223 <s223/s223-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with zero guide field, :math:`256x128` domain
  * - :doc:`s224 <s224/s224-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.25, :math:`256x128` domain
  * - :doc:`s225 <s225/s225-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.5, :math:`256x128` domain
  * - :doc:`s226 <s226/s226-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.75, :math:`256x128` domain
  * - :doc:`s227 <s227/s227-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 1.0, :math:`256x128` domain
  * - :doc:`s228 <s228/s228-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 2.0, :math:`256x128` domain
  * - :doc:`s229 <s229/s229-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 5.0, :math:`256x128` domain
  * - :doc:`s230 <s230/s230-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 10.0, :math:`256x128` domain
  * - :doc:`s231 <s231/s231-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with zero guide field, :math:`512x256` domain
  * - :doc:`s232 <s232/s232-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.25, :math:`512x256` domain
  * - :doc:`s233 <s233/s233-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.5, :math:`512x256` domain
  * - :doc:`s234 <s234/s234-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 0.75, :math:`512x256` domain
  * - :doc:`s235 <s235/s235-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge with guide field 1.0, :math:`512x256` domain
  * - :doc:`s236 <s236/s236-gemguide-5m>` 
    - Same as 235, except using open BCs instead of walls
  * - :doc:`s237 <s237/s237-gemguide-5m>` 
    - Same as 235 (guide field 1.0), except with 5% perturbation instead of 10%.
  * - :doc:`s238 <s238/s238-gemguide-5m>` 
    - Two-fluid 5-moment GEM challenge, :math:`768\times768` on :math:`25\times 25` domain and open BCs
  * - :doc:`s239 <s239/s239-gemguide-5m>` 
    - Same as s238, except on a :math:`252\times 252` grid.
  * - :doc:`s240 <s240/s240-gemguide-5m>` 
    - Same as s238, except on a :math:`512\times 512` grid.
  * - :doc:`s241 <s241/s241-gemguide-5m>` 
    - Open domain reconnection on :math:`50\times 25` domain, :math:`1024\times 512` grid.
  * - :doc:`s242 <s242/s242-dg-diffuse>` 
    - Diffusion equation with RDG scheme, using 8 cells, polyOrder 1.
  * - :doc:`s243 <s243/s243-dg-diffuse>` 
    - Same as 242, except with 16 cells.
  * - :doc:`s244 <s244/s244-dg-diffuse>` 
    - Same as 242, except with 32 cells.
  * - :doc:`s245 <s245/s245-dg-diffuse>` 
    - Same as 242, except with 64 cells.
  * - :doc:`s246 <s246/s246-dg-diffuse>` 
    - Diffusion equation with RDG scheme, using 4 cells, polyOrder 2.
  * - :doc:`s247 <s247/s247-dg-diffuse>` 
    - Same as 246, except with 8 cells.
  * - :doc:`s248 <s248/s248-dg-diffuse>` 
    - Same as 246, except with 16 cells.
  * - :doc:`s249 <s249/s249-dg-diffuse>` 
    - Same as 249, except with 32 cells.
  * - :doc:`s250 <s250/s250-dg-diffuse-2d>` 
    - Diffusion equation with RDG scheme in 2D, using :math:`8\times 8` cells, polyOrder 1.
  * - :doc:`s251 <s251/s251-dg-diffuse-2d>` 
    - Same as s250, except on :math:`16\times 16` grid.
  * - :doc:`s252 <s252/s252-dg-diffuse-2d>` 
    - Same as s250, except on :math:`32\times 32` grid.
  * - :doc:`s253 <s253/s253-dg-diffuse-2d>` 
    - Same as s250, except on :math:`64\times 64` grid.
  * - :doc:`s254 <s254/s254-free-stream-bounded>` 
    - Free streaming of particles on a bounded domain.
  * - :doc:`s255 <s255/s255-vlasov-fp-bounded>` 
    - Vlasov equation with fixed potential on a bounded domain.
  * - :doc:`s256 <s256/s256-gemguide-5m>` 
    - Open domain reconnection on :math:`50\times 25` domain, :math:`1536\times 768` grid.
  * - :doc:`s257 <s257/s257-blobs>` 
    - Blob simulation on :math:`30\times 20` domain and :math:`96\times 64` grid, no diffusion.
  * - :doc:`s258 <s258/s258-blobs>` 
    - Blob simulation on :math:`60\times 40` domain and :math:`192\times 128` grid, no diffusion.
  * - :doc:`s259 <s259/s259-blobs>` 
    - Same as s257, except with Rayleigh number :math:`1\times10^{4}`
  * - :doc:`s260 <s260/s260-blobs>` 
    - Same as s259, except on :math:`192\times 128` grid.
  * - :doc:`s261 <s261/s261-blobs>` 
    - Same as s260, except with Rayleigh number :math:`1\times10^{6}`
  * - :doc:`s262 <s262/s262-blobs>` 
    - Same as s260, except with Rayleigh number :math:`1\times10^{8}`
  * - :doc:`s263 <s263/s263-blobs>` 
    - Same as s260, except on a :math:`60\times 40` domain.
  * - :doc:`s264 <s264/s264-blobs>` 
    - Same as s260, except on a :math:`60\times 40` domain, :math:`384\times 256` grid.
  * - :doc:`s265 <s265/s265-blobs>` 
    - Same as s260, except on a :math:`60\times 40` domain and fixed BCs in X and periodic in Y.
  * - :doc:`s266 <s266/s266-blobs>` 
    - Same as s265, except on a :math:`384\times 256` grid.
  * - :doc:`s267 <s267/s267-blobs>` 
    - Same as s266, except on a :math:`450\times 300` grid.
  * - :doc:`s268 <s268/s268-blobs>` 
    - Same as s266, except with central fluxes.
  * - :doc:`s269 <s269/s269-free-stream-bounded>` 
    - Free streaming of particles on a bounded domain, with partially reflecting walls. See s254.      
  * - :doc:`s270 <s270/s270-gemguide-5m>` 
    - Open domain reconnection on :math:`50\times 25` domain, :math:`250\times 125` grid.
  * - :doc:`s271 <s271/s271-gemguide-5m>` 
    - Open domain reconnection on :math:`50\times 25` domain, :math:`500\times 250` grid.
  * - :doc:`s272 <s272/s272-gem-tenmom>` 
    - GEM challenge problem with 10-momnent equations, :math:`256\times 128` grid.
  * - :doc:`s273 <s273/s273-gem-tenmom>` 
    - GEM challenge problem with 10-momnent equations, :math:`512\times 256` grid.
  * - :doc:`s274 <s274/s274-gemguide-5m>` 
    - GEM challenge problem with 5-moment equations, :math:`512\times 256` grid and open domain (see s231)
  * - :doc:`s275 <s275/s275-gem-tenmom>` 
    - Open domain reconnection with 10-moment equations, :math:`500\times 250` grid. (See s271)
  * - :doc:`s276 <s276/s276-gem-tenmom>` 
    - Open domain reconnection with 10-moment equations, :math:`1000\times 500` grid. (See s275)
  * - :doc:`s277 <s277/s277-gem-tenmom>` 
    - Same as s275 expect collision frequency of 1.0.
  * - :doc:`s278 <s278/s278-gem-tenmom>` 
    - Same as s275 expect collision frequency of 0.0.
  * - :doc:`s279 <s279/s279-gem-tenmom>` 
    - Same as s277 expect on :math:`1000\times 500` grid.
  * - :doc:`s280 <s280/s280-gemguide-5m>` 
    - Same as s271 (5-moment) except on :math:`1000\times 500` grid.
  * - :doc:`s281 <s281/s281-gem-tenmom>` 
    - Open domain reconnection with 10-moment equations, :math:`500\times 250` grid with local relaxation heat-flux. (See s271)
  * - :doc:`s282 <s282/s282-gem-tenmom>` 
    - Same as s281, except on :math:`1000\times 500` grid.
  * - :doc:`s283 <s283/s283-gem-tenmom>` 
    - Same as s281, except on :math:`250\times 125` grid
  * - :doc:`s284 <s284/s284-pulsebox-wave>` 
    - EM pulse in box (see s57). This is to test flux and electric field relations.
  * - :doc:`s285 <s285/s285-pulsebox-wave>` 
    - EM pulse in box (see s57 and s284). With in-plane electric field.
  * - :doc:`s286 <s286/s286-harris-tenmom>` 
    - Ten-moment doubly periodic Harris sheet simulation to compare with PSC. On :math:`250\times 125` grid.
  * - :doc:`s287 <s287/s287-harris-tenmom>` 
    - Same as s286, except on a :math:`1000\times 500` grid.
  * - :doc:`s288 <s288/s288-ecell-fdtd>` 
    - FDTD Maxwell solver, in preparation for FD/FV two-fluid scheme. (See s63).
  * - :doc:`s289 <s289/s289-pulse-box-5m>` 
    - A plamsa in a box, perturbed with a EM pulse. Another attempt to nail down div(E) problems.
  * - :doc:`s290 <s290/s290-harris-tenmom>` 
    - Same as s286, except on a :math:`500\times 250` grid.
  * - :doc:`s291 <s291/s291-pulsebox-wave>` 
    - EM pulse in box (see s284). Using dimensional splitting.
  * - :doc:`s292 <s292/s292-harris-tenmom>` 
    - Same as s286, except using dimensional splitting.
  * - :doc:`s293 <s293/s293-harris-tenmom>` 
    - Same as s286, except using explicit RK4 source update (this just for comparsion with s286).
  * - :doc:`s294 <s294/s294-harris-tenmom>` 
    - Same as s292, except on a :math:`500\times 250` grid.
  * - :doc:`s295 <s295/s295-harris-tenmom>` 
    - Same as s292, except on a :math:`500\times 250` grid and :math:`100\times 50` domain.
  * - :doc:`s296 <s296/s296-harris-tenmom>` 
    - Same as s295, except on a :math:`1000\times 500` grid and :math:`100\times 50` domain.
  * - :doc:`s297 <s297/s297-gem-tenmom>` 
    - Same as s282 (open domain reconnection), except using a dimensional splitting.
  * - :doc:`s298 <s298/s298-harris-tenmom>` 
    - Same as s296 (double periodic), except with collisional relaxation (elcCollisionFreq 1.0).
  * - :doc:`s299 <s299/s299-harris-tenmom>` 
    - Same as s296 (double periodic), except on :math:`2000\times 1000`.
  * - :doc:`s300 <s300/s300-harris-tenmom>` 
    - Same as s296 (double periodic), except with CFL of 0.95.
  * - :doc:`s301 <s301/s301-5m-double-periodic>` 
    - Doubly periodic reconnection with 5-moment model. :math:`1000\times 500` grid.
  * - :doc:`s302 <s302/s302-5m-double-periodic>` 
    - Doubly periodic reconnection with 5-moment model. :math:`1000\times 500` grid, on :math`200\times 100` grid.
  * - :doc:`s303 <s303/s303-5m-double-periodic>` 
    - Same as s302, except on :math:`2000\times 1000`.
  * - :doc:`s304 <s304/s304-5m-gem>` 
    - Two-fluid 5-moment GEM challenge with zero guide field, :math:`512x256` grid, dimensional splitting
  * - :doc:`s305 <s305/s305-5m-gem>` 
    - Same as s304, except on double the domain and :math:`1024x512` grid.
  * - :doc:`s306 <s306/s306-5m-gem>` 
    - Same as s304, except :math:`1024x512` grid.
  * - :doc:`s307 <s307/s307-harris-tenmom>` 
    - Same as s296 (double periodic 10M), except using :math:`\eta J` term and on :math:`500\times 250` grid.
  * - :doc:`s308 <s308/s308-harris-tenmom>` 
    - Same as s296 (double periodic 10M), except using a perturbation of 0.5.
  * - :doc:`s309 <s309/s309-5m-double-periodic>` 
    - Same as 301 (double periodic 5M), except enforcing quasi-neutrality
  * - :doc:`s310 <s310/s310-5m-double-periodic>` 
    - Same as 309, with div(E) also enforced using hyperbolic cleaning.
  * - :doc:`s311 <s311/s311-5m-gem>` 
    - Same as 304, shorter wih more output and odd number of cells
  * - :doc:`s312 <s312/s312-5m-gem>` 
    - Same as 311, except with Ey initialized to equilibrium value.
  * - :doc:`s313 <s313/s313-5m-gem>` 
    - Same as 312, except with :math:`1601\times 801` mesh.
  * - :doc:`s314 <s314/s314-5m-gem>` 
    - Same as s311, except both electrons and ions carry current, on a :math:`1001\times 501` grid.
  * - :doc:`s315 <s315/s315-harris-tenmom>` 
    - Same as s296 (double tearing), except both electrons and ions carry current, and using different perturbation.
  * - :doc:`s316 <s316/s316-harris-tenmom>` 
    - Same as s315, except using :math:`k=1/\delta_i`.
  * - :doc:`s317 <s317/s317-5m-gem>` 
    - Same as s314, except fewer output files, and run longer
  * - :doc:`s318 <s318/s318-gem-tenmom>` 
    - Same s297 (open domain reconnection, 10M), except electrons and ions carry current
  * - :doc:`s319 <s319/s319-gemguide-5m>` 
    - Same as s280 (open domain reconnection, 5M), dimensional splitting, both electrons and ions carry current.
  * - :doc:`s320 <s320/s320-gem-tenmom>` 
    - Standard GEM reconnection, with 10M model and periodic in X and walls in Y
  * - :doc:`s321 <s321/s321-gem-tenmom>` 
    - Same as s320, except with wave number of ionSkinDepth.
  * - :doc:`s322 <s322/s322-harris-tenmom>` 
    - Same as s316, except using :math:`1002\times 502` mesh.
  * - :doc:`s323 <s323/s323-harris-tenmom>` 
    - Same as s316, using no ion relaxation.
  * - :doc:`s324 <s324/s324-harris-tenmom>` 
    - Same as s316, with new equilibrium for double tearing (:math:`k=1/\delta_i`)
  * - :doc:`s325 <s325/s325-harris-tenmom>` 
    - Same as s324, expect with :math:`k=1/\delta_e`
  * - :doc:`s326 <s326/s326-iso-gem>` 
    - Two-fluid, isothermal standard GEM reconnection.
  * - :doc:`s327 <s327/s327-gem-tenmom>` 
    - 10-moment with anisotropic heat relaxation closure.
  * - :doc:`s328 <s328/s328-lgd-diffuse>` 
    - Diffusion with Local DG scheme, same parameters as s245.
  * - :doc:`s329 <s329/s329-modal-dg-diffuse>` 
    - Diffusion with modal DG scheme. Run to steady-state with sin(x) source. LDG-L scheme.
  * - :doc:`s330 <s330/s330-modal-dg-diffuse>` 
    - Same as s329, except with LDG-R scheme.
  * - :doc:`s331 <s331/s331-modal-dg-diffuse>` 
    - Same as s329, except with LDG-S scheme.
  * - :doc:`s332 <s332/s332-modal-dg-diffuse>` 
    - Same as s329, except with RDG scheme.
  * - :doc:`s333 <s333/s333-modal-dg-diffuse>` 
    - Decaying temperature problem, LDG-L scheme, 8 cells.
  * - :doc:`s334 <s334/s334-modal-dg-diffuse>` 
    - Same as s333, except with LDG-R scheme, 8 cells.
  * - :doc:`s335 <s335/s335-modal-dg-diffuse>` 
    - Same as s333, except with LDG-S scheme, 8 cells.
  * - :doc:`s336 <s336/s336-modal-dg-diffuse>` 
    - Same as s333, except with RDG scheme, 8 cells.
  * - :doc:`s337 <s337/s337-modal-dg-diffuse>` 
    - Same as s333 (LDG-L), except 16 cells.
  * - :doc:`s338 <s338/s338-modal-dg-diffuse>` 
    - Same as s334 (LDG-R), except 16 cells.
  * - :doc:`s339 <s339/s339-modal-dg-diffuse>` 
    - Same as s335 (LDG-S), except 16 cells.
  * - :doc:`s340 <s340/s340-modal-dg-diffuse>` 
    - Same as s340 (RDG), except 16 cells.
  * - :doc:`s341 <s341/s341-modal-dg-diffuse>` 
    - Same as s333 (LDG-L), except 32 cells.
  * - :doc:`s342 <s342/s342-modal-dg-diffuse>` 
    - Same as s334 (LDG-R), except 32 cells.
  * - :doc:`s343 <s343/s343-modal-dg-diffuse>` 
    - Same as s335 (LDG-S), except 32 cells.
  * - :doc:`s344 <s344/s344-modal-dg-diffuse>` 
    - Same as s340 (RDG), except 32 cells.
  * - :doc:`s345 <s345/s345-modal-dg-diffuse>` 
    - Same as s333 (LDG-L), except 4 cells. (This is out of whack, but needed for proposal)
  * - :doc:`s346 <s346/s346-modal-dg-diffuse>` 
    - Same as s340 (RDG), except 4 cells. (This is out of whack, but needed for proposal)

