:Author: Ammar Hakim
:Date: 29th December 2013
:Completed: 
:Last Updated:

JE21: Testing a solver for linearized electromagnetic GK equations
==================================================================

.. contents::

.. warning::

  This note is incomplete.

  - I need to put full dimensional EM equations.
  - Put dispersion relations and link to Python code.
  - I may need to redo some plots, as perhaps I have some
    normalizations confused.

  There is also some notational confusion between
  :math:`f(z,p_\parallel,t)` and :math:`f(z,v_\parallel,t)`.

In this note I test a solver for the linearized electromagentic
gyrokinetic equations. This system is described by the Vlasov equation
for the electrons

.. math::

  \frac{\partial f}{\partial t} + \{f,H\} = 0

where :math:`f(z,p_\parallel,t)` is the electron distribution function
[#dist-function]_ and the Hamiltonian is given by

.. math::
  :label: nonlin-H

  H = \frac{1}{2m_e}(p_\parallel-q_e A_\parallel)^2 + q_e \phi

where :math:`A_\parallel` is the parallel magnetic potential and
:math:`\phi` is the electrostatic potential, :math:`m_e` is the
electron mass and :math:`q_e = -|e|` is the electron charge. The ions
are assumed to be stationary.

Normalized linear system
------------------------

To linearize the system we first linearize the Hamiltonian around a
field-free equilibirum to write

.. math::

  H &= H_0 + H_1.

Here :math:`H` is the linearized Hamiltonian. The zeroth-order
Hamiltonian describes free-streaming

.. math::

  H_0 = \frac{1}{2m_e} p_\parallel^2

and the first-order Hamiltonian describes motion of the perturbations

.. math::

  H_1 = -\frac{q_e}{m_e}p_\parallel A_\parallel + q_e\phi.

Writing :math:`f(z,p_\parallel,t) = f_0(z,p_\parallel) +
f_1(z,p_\parallel,t)` the linearized Vlasov equation is given by

.. math::

  \frac{\partial f_1}{\partial t} + \{f_1,H_0\} + \{f_0,H_0+H_1\} = 0.

Note that only the zeroth-order Hamiltonian acts on the perturbed
distribution, while the full linearized Hamiltonian acts on the
equilibirum distribution. The linearized electromagnetic field
equations become: Ampere's law

.. math::

  k_\perp^2 A_\parallel = \mu_0 q_e \int v_\parallel f_1\thinspace dv_\parallel

and the GK Poisson equation

.. math::

  \frac{q_i n_{0i}}{T_{0i}}
  k_\perp^2\rho_i^2 \phi
  =
  q_e \int f_1\thinspace dv_\parallel

where :math:`q_i = Z_i |e|` is ion charge, and :math:`v_\parallel =
(p_\parallel-q_e A_\parallel)/m_e` is the parallel velocity.

We normalize the equations as follows. We pick :math:`q_e=-1`,
:math:`m_e=1` and :math:`v_{te}=1`. To simpify the notation we now use
:math:`w` as the momentum space coordinate. The form of normalized
Vlasov equation remains unchanged, with the Hamiltonians given by
:math:`H_0 = w^2/2` and

.. math::

  H_1 = w A_\parallel - \phi.

Note that in this normalization, effectively :math:`k_\perp^2` is
normalized to :math:`Z_i\rho_i^2` and :math:`\hat{\beta} \equiv (\beta_e/2)
m_i/m_e`, where :math:`\beta_e` is the electron plasma-beta. In these
units, the Alfven speed is :math:`v_A=1/\sqrt{\hat{\beta}}`. The normalized
equilibrium distribution is :math:`f_0 = e^{-w^2/2}/\sqrt{2\pi}`.

The EM equations become

.. math::

 k_\perp^2 \phi &= -\int f_1\thinspace dw\\
 (k_\perp^2+\hat{\beta}) A_\parallel &= -\hat{\beta} \int w f_1\thinspace dw.

Note that the normalized system is identical to the one described in
the appendix of Belli and Hammett 2005 [#belli-hammett-2005]_.

To solve this system with Gkeyll we use two Poisson Bracket updaters,
one acting on the perturbed distribution and the other acting on the
equilibrium distribution. The perturbed Hamiltonian is computed from
the EM fields, after they are projected onto continious basis
functions.

.. note::

  In the implementation used here, the Hamiltonian is computed by
  evaluating it at the nodes of the nodal basis functions. This
  ensures continuity, however, is not an :math:`L_2` fit onto the
  continuous basis functions. This could results in a small loss in
  accuracy, although total energy will still be conserved.

Electron acoustic waves
-----------------------

We first look at the electrostatic limit, i.e. set :math:`\hat{\beta}=0`. In
this limit the system supports electron acoustic waves, which are
increasingly Landau damped as :math:`k_\perp^2` increases. 

The first set of simulations were performed with a piece-wise
first-order polynomial basis functions on a :math:`16\times 32`
grid. The second set of simulations were performed on the same grid,
except using piece-wise second-order polynomials. 

.. figure:: freq-damp-ion-sound.png
  :width: 100%
  :align: center

  Frequency (magenta, left axis) and damping (green, right axis) for
  electron acoustic waves. Solid dots are simulation results on a
  :math:`16\times 32` grid with piece-wise first-order polynomial
  basis functions. Note that the resolution is rather coarse, and the
  discretization errors are particularly noticeable in the damping
  rates. In these simulations :math:`\hat{\beta}=0.0`, and
  :math:`k_\perp^2=0.01,\ldots,1.0`. See simulations [:doc:`347
  <../../sims/s347/s347-lin-em-gke>`] to [:doc:`s354
  <../../sims/s354/s354-lin-em-gke>`] for details.

.. figure:: p2-freq-damp-ion-sound.png
  :width: 100%
  :align: center

  Same as the previous figure, except using piece-wise second order
  polynomial basis functions. The damping rates are much better
  predicted than in the first-order polynomial case, however, the
  simulations take twice as long. See simulations [:doc:`s355
  <../../sims/s355/s355-lin-em-gke>`] to [:doc:`s362
  <../../sims/s362/s362-lin-em-gke>`] for details.

Shear Kinetic Alfven waves
--------------------------

In the next set of calculations, we look at the case in which EM terms
are included, i.e. :math:`\hat{\beta}>0.0`. In this case the system supports
shear kinetic Alfven waves (KAWs), which asymptote to undamped waves
as :math:`k_\perp\rightarrow 0`. 

Case when :math:`k_\perp^2=0.1`
+++++++++++++++++++++++++++++++

For first set of tests I hold :math:`k_\perp^2=0.1` and vary
:math:`\hat{\beta}=0.1,\ldots,10.0`. For all simulations piece-wise
second-order basis functions on a grid of :math:`16\times 32` were
used. One of the reasons to use a second-order polynomial basis
functions is that it delays recurrence issues, rather severe in the
lower-order basis function case. Eventually, we will implement a
hyper-collision term to damp out the recurrence, but this has not been
tested yet.

The results are shown in the following figure.

.. figure:: freq-damp-shear-alf-kp-0p1-beta-scan.png
  :width: 100%
  :align: center

  Frequency (magenta, left axis) and damping (green, right axis) for
  shear kinetic Alfven waves (KAWs). Solid dots are simulation results
  on a :math:`16\times 32` grid with piece-wise second-order
  polynomial basis functions. Note that the damping rates do not agree
  very well for the :math:`\hat{\beta}=10.0` case. See simulations [:doc:`s363
  <../../sims/s363/s363-lin-em-gke>`] to [:doc:`s369
  <../../sims/s369/s369-lin-em-gke>`] for details.

Case when :math:`k_\perp^2=0.05`
+++++++++++++++++++++++++++++++

In this set of tests I hold :math:`k_\perp^2=0.05` and vary
:math:`\hat{\beta}=0.1,\ldots,10.0`. All other parameters are the same as
for the :math:`k_\perp^2=0.1` case.

The results are shown in the following figure.

.. figure:: freq-damp-shear-alf-kp-0p05-beta-scan.png
  :width: 100%
  :align: center

  Frequency (magenta, left axis) and damping (green, right axis) for
  shear kinetic Alfven waves (KAWs). Solid dots are simulation results
  on a :math:`16\times 32` grid with piece-wise second-order
  polynomial basis functions. See simulations [:doc:`s370
  <../../sims/s370/s370-lin-em-gke>`] to [:doc:`s370
  <../../sims/s370/s370-lin-em-gke>`] for details.

Case when :math:`k_\perp^2=0.01`
+++++++++++++++++++++++++++++++

In this set of tests I hold :math:`k_\perp^2=0.01` and vary
:math:`\hat{\beta}=0.1,\ldots,10.0`. All other parameters are the same as
for the :math:`k_\perp^2=0.1` case.

The results are shown in the following figure.

.. figure:: freq-damp-shear-alf-kp-0p01-beta-scan.png
  :width: 100%
  :align: center

  Frequency (magenta, left axis) and damping (green, right axis) for
  shear kinetic Alfven waves (KAWs). Solid dots are simulation results
  on a :math:`16\times 32` grid with piece-wise second-order
  polynomial basis functions. See simulations [:doc:`s370
  <../../sims/s370/s370-lin-em-gke>`] to [:doc:`s370
  <../../sims/s370/s370-lin-em-gke>`] for details.

Case when :math:`\hat{\beta}=10`
+++++++++++++++++++++++++++++++

In the next set of tests :math:`\hat{\beta}=10`, while
:math:`k_\perp^2=10^{-4},\ldots,1.0`. A grid of :math:`16\times 32`
cells, with piece-wise second-order polynomial basis functions were
used. For the case :math:`k_\perp^2=10^{-4}`, a fixed time-step
:math:`\Delta t = 10^{-3}` was used to avoid continuous adjustments
from the Poisson Bracket algorithm. For all other tests the largest
time-step allowed by the CFL condition was used.

.. figure:: freq-damp-shear-alf-beta-10-kp-scan.png
  :width: 100%
  :align: center

  Frequency for shear kinetic Alfven waves (KAWs), with fixed
  :math:`\hat{\beta}=10`. Solid dots are simulation results on a
  :math:`16\times 32` grid with piece-wise second-order polynomial
  basis functions. Gkeyll predicts the frequency to at least two
  significant figures in each case.

Comparison with nonlinear solver
--------------------------------

In this section, I test a solver for the full nonlinear
electromagnetic gyrokinetic equations. Only the electron dynamics is
retained, i.e. the ions are assumed to be stationary. The evolution
equation is

.. math::

  \frac{\partial f_1}{\partial t} + \{f_1+f_0, H\} = 0

where, now the Hamiltonian is given by :eq:`nonlin-H`. The field
equations remain unchanged. The current implementation in Gkeyll
evolves only the perturbed distribution function, however the Poisson
Bracket operator is now applied to the total distribution and the
nonlinear Hamiltonian.

The follow table compares the values obtained for the frequencies for
the :math:`\beta=10.0` case and various values of :math:`k_\perp`.

.. list-table::
  :header-rows: 1
  :widths: 10,30,30,30
   
  * - :math:`k_\perp^2`
    - Linear Code
    - Nonlinear Code
    - Exact
  * - :math:`1 \times 10^{-2}`
    - 0.22465
    - 0.22468
    - 0.22466
  * - :math:`5\times 10^{-2}`
    - 0.22856
    - 0.22856
    - 0.22879

This test shows that both the nonlinear and linear codes give the same
answer to two significant figures, showing that there is no
fundamental problem with the DG algorithm for nonlinear EM/GKE.

-----

.. [#dist-function] The distribution function is for the guiding
   centers. However, in this note the zero gyro-radius approximation
   is used for the electrons, and the particle and guiding center
   distributions coincide.

.. [#belli-hammett-2005] Belli, E. A., & Hammett, G. W. "A numerical
   instability in an ADI algorithm for gyrokinetics", *Computer
   Physics Communications*, **172** (2),
   119–132, 2005. doi:10.1016/j.cpc.2005.06.007
