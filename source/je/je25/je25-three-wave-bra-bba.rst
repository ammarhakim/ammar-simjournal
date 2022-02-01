:Author: Ammar Hakim
:Date: August 4th 2014
:Completed: 
:Last Updated:

JE25: Three-wave models for Backward Raman and Brillouin Amplification
======================================================================

.. note::

  This note is written in conjunction with Updesh Verma, a C.V. Raman
  Fellow visiting Princeton University (2013-2014).

In this note I test a solver for the three-wave resonant interaction
(TWRI) model, and apply it to Backward Raman Amplification (BRA) and
Backward Brillouin Amplification (BBA) problems. The TWRI equations
are written in symmetric form

.. math::

  \left(
    \frac{\partial}{\partial t} + v_1 \frac{\partial}{\partial x}
  \right) E_1 & = c_1 E_2^* E_3^* \\
  \left(
    \frac{\partial}{\partial t} + v_2 \frac{\partial}{\partial x}
  \right) E_2 & = c_2 E_1^* E_3^* \\
  \left(
    \frac{\partial}{\partial t} + v_3 \frac{\partial}{\partial x}
  \right) E_3 & = c_3 E_1^* E_2^*

where, for :math:`i=1,2,3`, :math:`E_i` are (complex) wave amplitudes,
:math:`v_i` are wave velocities and :math:`c_i` are real coupling
constants.

It should be mentioned that for the BBA in the *strong coupling*
regime, the final equation above is replaced instead by

.. math::

  \left(
    \frac{\partial^2}{\partial t^2} - v_3^2 \frac{\partial^2}{\partial x^2}
  \right) E_3 & = c_3 E_1^* E_2^*  

I.e, a *wave-equation* rather than a advection equation. See, for
example, [#lehman-2013]_. To solve this second order equation, an
auxiliary variable is introduced

.. math::

  w = \frac{\partial E_3}{\partial t} - v_3 \frac{\partial E_3}{\partial x}.

Using this, we can rewrite the second-order wave equation as a
first-order equation

.. math::

 \frac{\partial w}{\partial t} + v_3 \frac{\partial w}{\partial x}
   = c_3 E_1^* E_2^*

These two coupled equations are then solved (instead of the single
second-order equation) using the standard advection equation solver.

The coupled system of equations are solved in Gkeyll using a Strang
operator splitting method. In each time-step, the source coupling is
updated first by half time-step, then the advection terms are updated
by a full time-step, and finally the source terms are updated again by
half time-step. To solve the source coupling ODEs two solvers are
implemented: an Crank-Nicholson implicit solver and an explicit
Runge-Kutta 4th order scheme. Both of these ODE solvers, when used in
the full TWRI solver, given nearly identical results, and hence can be
used interchangeably.

As is well known (see for example [#kaup-1979]_) the nature of the
solution, for :math:`v_1<v_2<v_3` depends on the *signs* of the
coefficients :math:`c_1,c_2,c_3`. Let
:math:`\eta_i=\mathrm{sign}(c_i)`, then the cases can be distinguished
by

.. math::

  \eta_1 \eta_2 \eta_3 = \eta, \quad \eta = \pm 1

The three cases identified are the *explosive*, *simulated
backscatter* and *soliton exchange*. For :math:`\eta=1`, we have

.. math::

  (\eta_1,\eta_2,\eta_3) &= (+,+,+) &\quad\mathrm{explosive} \\
  (\eta_1,\eta_2,\eta_3) &= (+,-,-) &\quad\mathrm{backscatter} \\
  (\eta_1,\eta_2,\eta_3) &= (-,-,+) &\quad\mathrm{backscatter} \\
  (\eta_1,\eta_2,\eta_3) &= (-,+,-) &\quad\mathrm{soliton\ exchange}

The case of interest here is the simulated backscatter, which can be
used for laser pulse amplification using the BRA and BBA processes.


Backward Raman Amplification
----------------------------

For BRA, the 

.. warning::

  I hope to complete this note in the near future.

References
----------

.. [#lehman-2013] Lehmann, G., & Spatschek, K. H.. "Nonlinear
   Brillouin amplification of finite-duration seeds in the strong
   coupling regime". *Physics of Plasmas*, **20**
   (7), 073112. doi:10.1063/1.4816030

.. [#kaup-1979] Kaup, D. J., Reiman, A., & Bers, A. "Space-time
   evolution of nonlinear three-wave interactions. I. Interactions in
   a homogeneous medium". *Reviews of Modern Physics*, **51** (2), 275.
