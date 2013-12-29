:Author: Ammar Hakim
:Date: 29th December 2013
:Completed: 
:Last Updated:

JE21: Testing a solver for linearized electromagnetic GK equations
==================================================================

In this note I test a solver for the linearized electromagentic
gyrokinetic equations. This system is described by the Vlasov equation
for the electrons

.. math::

  \frac{\partial f}{\partial t} + \{f,H\} = 0

where :math:`f(z,p_\parallel,t)` is the electron distribution function
[#dist-function]_ and the Hamiltonian is given by

.. math::

  H = \frac{1}{2m_e}(p_\parallel-q_e A_\parallel)^2 + q_e \phi

where :math:`A_\parallel` is the parallel magnetic potential and
:math:`\phi` is the electrostatic potential, :math:`m_e` is the
electron mass and :math:`q_e = -|e|` is the electron charge. The ions
are assumed to be stationary.

.. warning::

  Put full EM equations here. Will require some re-wording of the
  following section.

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
normalized to :math:`Z_i\rho_i^2` and :math:`\beta \equiv (\beta_e/2)
m_i/m_e`, where :math:`\beta_e` is the electron plasma-beta. Further,
the equilibrium distribution is :math:`f_0 = e^{-w^2/2}/\sqrt{2\pi}`.

The EM equations become

.. math::

 k_\perp^2 \phi &= -\int f_1\thinspace dw\\
 (k_\perp^2-\beta) A_\parallel &= -\beta \int w f_1\thinspace dw.

Note that the normalized system is identical to the one described in
the appendix of Belli and Hammett 2005 [#belli-hammett-2005]_.

To solve this system with Gkeyll we use two Poisson Bracket updaters,
one acting on the perturbed distribution and the other acting on the
equilibrium distribution. The perturbed Hamiltonian is computed from
the EM fields, after they are projected onto continious basis
functions.

-----

.. [#dist-function] The distribution function is for the guiding
   centers. However, in this note the zero gyro-radius approximation
   is used for the electrons, and the particle and guiding center
   distributions coincide.

.. [#belli-hammett-2005] Belli, E. A., & Hammett, G. W. "A numerical
   instability in an ADI algorithm for gyrokinetics", *Computer
   Physics Communications*, **172** (2),
   119–132, 2005. doi:10.1016/j.cpc.2005.06.007
