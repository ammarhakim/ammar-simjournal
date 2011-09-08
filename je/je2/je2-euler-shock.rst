JE2: Benchmarking Two Finite-Volume Schemes for 1D Euler Equations
==================================================================

:Author: Ammar Hakim
:Date: September 6th 2011

.. contents::

Overview of problems and schemes
--------------------------------

In this entry two algorithms, the wave-propagation algorithm and the
`MUSCL-Hancock scheme <http://ammar-hakim.org/hancock-muscl.html>`_,
are benchmarked against a series of 1D shock-tube problems. For most
of these problems exact solutions can be computed using an exact
Riemann solver, described, for example, in [Kulikovskii2001]_.

The wave-propagation scheme is implemented in the
``WavePropagationUpdater`` in the ``slvrs`` directory, while the
MUSCL-Hancock scheme is implemented in the ``MusclHancock1DUpdater``
in the `proto` directory. Note that as this point the MUSCL-Hancock
scheme is *not* production quality but is hard-coded to solve just the
1D Euler equations. The code to compute the exact solution is located
in the ``sims/code/exactrp`` directory.

The problems solved in this entry are quite difficult and really tax
the ability of schemes to capture discontinuities accurately. In
particular, schemes have a lot of trouble in problems in which vacuum
states can form or exist in parts of the domain. Also, ensuring
positivity of density and pressure is a difficult problem.

Shock-tube problems
-------------------

Eight shock-tube problems are solved. Each problem is described below
and results of comparing both schemes with the exact solution are
shown. Internal energy in the following is computed as
:math:`p/(\gamma-1)\rho`.

Problem 1
+++++++++

Sod-shock with sonic point in rarefaction. Domain is :math:`x \in
[0,1]`, discretized with 100 cells. Gas adiabatic constant of 1.4 was
used. Simulation is initialized with a shock at :math:`x=0.3`, with
left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1 \\
      0.75 \\
      1.0
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      0.125 \\
      0.0 \\
      0.1
    \end{matrix}
  \right].

and is run to :math:`t=0.2`.

.. figure:: s5-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s5] with exact
  solution (red) [s6] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

.. figure:: s7-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of MUSCL-Hancock solution (black) [s7] with exact
  solution (red) [s6] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

Problem 2
+++++++++

This problem has a near-vaccum near the location of the
discontinuity. Domain is :math:`x \in [0,1]`, discretized with 100
cells. Gas adiabatic constant of 1.4 is used. Simulation is
initialized with a shock at :math:`x=0.5`, with left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      -2.0 \\
      0.4
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      2.0 \\
      0.4
    \end{matrix}
  \right].

and is run to :math:`t=0.15`.

Both wave-propagation and MUSCL-Hancock **fail** on this problem. The
solution quickly develops negative pressure and density. A positivity
fix is required for both schemes (not implemented as of September 6
2011). First-order MUSCL-Hancock, however, works and results are shown
below.

.. figure:: s10-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of 1st-order MUSCL-Hancock solution (black) [s10] with
  exact solution (red) [s9] for density (top left), velocity (top
  right), pressure (bottom left) and internal energy (bottom
  right).

Problem 3
+++++++++

The 1D Noh problem. Domain is :math:`x \in [0,1]`, discretized with
100 cells. Gas adiabatic constant of :math:`5/3` is used. Simulation
is initialized with a shock at :math:`x=0.5`, with left and right
states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      1.0 \\
      10^{-6}
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      -1.0 \\
      10^{-6}
    \end{matrix}
  \right].

and is run to :math:`t=1.0`.

.. figure:: s11-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s11] with exact
  solution (red) [s12] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

The MUSCL-Hancock scheme **fails** on this problem. A positivity fix
needs to be implemented. However, the 1st-order MUSCL-Hancock scheme
works and results are shown below.

.. figure:: s13-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of 1st-order MUSCL-Hancock solution (black) [s13] with
  exact solution (red) [s12] for density (top left), velocity (top
  right), pressure (bottom left) and internal energy (bottom right).

Problem 4
+++++++++

1D Euler shock with a stationary contact discontinuity at
:math:`x=0.8`. Domain is :math:`x \in [0,1]`, discretized with 100
cells. Gas adiabatic constant of :math:`1.4` is used. Simulation is
initialized with a shock at :math:`x=0.8`, with left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      -19.59745 \\
      1000
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      -19.59745 \\
      0.01
    \end{matrix}
  \right].

and is run to :math:`t=0.012`.

.. figure:: s14-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s14] with exact
  solution (red) [s15] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

The MUSCL-Hancock scheme **fails** on this problem. Results with the
1st-order MUSCL-Hancock method is shown below.

.. figure:: s16-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of 1st-order MUSCL-Hancock solution (black) [s16] with
  exact solution (red) [s15] for density (top left), velocity (top
  right), pressure (bottom left) and internal energy (bottom right).

Problem 5
+++++++++

1D Euler shock with two strong shocks. Domain is :math:`x \in [0,1]`,
discretized with 100 cells. Gas adiabatic constant of :math:`1.4` is
used. Simulation is initialized with a shock at :math:`x=0.4`, with
left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      5.99924 \\
      19.5975 \\
      460.894
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      5.99242 \\
      -6.19633 \\
      46.0895
    \end{matrix}
  \right].

and is run to :math:`t=0.035`.

.. figure:: s17-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s17] with exact
  solution (red) [s18] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

.. figure:: s19-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of MUSCL-Hancock solution (black) [s19] with exact
  solution (red) [s18] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

Problem 6
+++++++++

1D Euler with a stationary contact discontinuity. Domain is :math:`x
\in [0,1]`, discretized with 100 cells. Gas adiabatic constant of
:math:`1.4` is used. Simulation is initialized with a shock at
:math:`x=0.5`, with left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.4 \\
      0.0 \\
      1.0
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      0.0 \\
      1.0
    \end{matrix}
  \right].

and is run to :math:`t=2.0`.

.. figure:: s20-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s20] with exact
  solution (red) [s21] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

.. figure:: s22-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of MUSCL-Hancock solution (black) [s22] with exact
  solution (red) [s21] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

Problem 7
+++++++++

1D Euler with a slowly moving contact discontinuity. Domain is
:math:`x \in [0,1]`, discretized with 100 cells. Gas adiabatic
constant of :math:`1.4` is used. Simulation is initialized with a
shock at :math:`x=0.5`, with left and right states

.. math::

  \left[
    \begin{matrix}
      \rho_l \\
      u_l \\
      p_l
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.4 \\
      0.1 \\
      1.0
    \end{matrix}
  \right],
  \qquad
  \left[
    \begin{matrix}
      \rho_r \\
      u_r \\
      p_r
    \end{matrix}
  \right]
  = 
  \left[
    \begin{matrix}
      1.0 \\
      0.1 \\
      1.0
    \end{matrix}
  \right].

and is run to :math:`t=2.0`.

.. figure:: s23-euler-shock-wave_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of wave-propagation solution (black) [s23] with exact
  solution (red) [s24] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right).

.. figure:: s25-euler-shock-muscl_exact_cmp.png
  :width: 100%
  :align: center

  Comparison of MUSCL-Hancock solution (black) [s25] with exact
  solution (red) [s24] for density (top left), velocity (top right),
  pressure (bottom left) and internal energy (bottom right). The
  oscillations in the velocity and pressure appear exaggerated by the
  plot scale but are actually very small.

Woodward-Collela blast wave problem
-----------------------------------

XXX

Conclusions
-----------

One of the aims of this note was to determine what modifications are
needed to the wave-propagation scheme and the MUSCL-Hancock scheme to
make them more robust and accurate. Note that the MUSCL-Hancock scheme
tested here is only a prototype version and the tests in conducted in
this entry will allow a better production quality solver.

The lessons learned are:

- The wave-propagation scheme needs a positivity fix. For this, a
  density and pressure floor should be added. More importantly, if the
  Roe averages lead to a NaN or negative pressure, the Roe fluxes
  should be replaced (automatically) with a diffusive, but positivity
  preserving, Rusanov (Lax) flux.

- The MUSCL-Hancock scheme needs a positivity fix also: essentially,
  if the predicted edge values are negative the slope in the cell
  should be simply set to zero.

- More accurate (than Rusanov flux) numerical flux needs to be
  implemented. An HLLC flux will help reduce the diffusion as compared
  to the wave-propagation scheme.

References
----------

.. [Kulikovskii2001] Andrei G. Kulikoviskii and Nikolai V. Pogorelov
   and Andrei Yu. Semenov, *Mathematical Aspects of Numerical
   Solutions of Hyperbolic Systems*, Chapman and Hall/CRC, 2001.