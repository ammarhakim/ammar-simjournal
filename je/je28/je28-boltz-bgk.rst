:Author: Ammar Hakim
:Date: 30nd Jan 2016
:Completed: 
:Last Updated:

JE28: Some tests for Boltzmann-BGK equations
============================================

.. note::

   This note is written with Petr Cagas, a graduate student at
   Virginia Tech. Petr wrote the code to compute a Maxwellian from
   moments. This updater at present does not conserve the
   moments. I.e. moments computed *from* the Maxwellian won't exactly
   match the moments used to *compute* the Maxwellian. We hope to fix
   this soon.

.. contents::

Introduction
------------

Gkeyll includes the ability to incorporate (as of this writing) both a
BGK and a Lenard-Bernstein (LB) collision operator. The LB operator is
modified to ensure momentum and energy conservation. In this note, I
test the BGK operator in the context of kinetic simulations of neutral
gas dynamics.

The neutral gas dynamics system, with the BGK operator, is given by

.. math::

  \frac{\partial f}{\partial t} + v \frac{\partial f}{\partial x} =
  \nu (f_M - f)

where :math:`f_M` is the Maxwellian distribution computed from the
moments of :math:`f(x,v,t)`:

.. math::

   f_M(x,v,t) = \frac{n}{\sqrt{2\pi v_{th}^2}} \exp(-(v-u)^2/2 v_{th}^2)

where

.. math::

   n \equiv \int f(x,v,t) \thinspace dv \\
   nu \equiv \int v f(x,v,t) \thinspace dv \\
   nu^2 + n v_{th}^2 \equiv \int v^2 f(x,v,t) \thinspace dv

The above equations are written in 1X/1V, however, the code works in
higher dimensions.

Test Problems
-------------

Problem 1: Relaxation of step function to Maxwellian
++++++++++++++++++++++++++++++++++++++++++++++++++++

In this test the relaxation of an initial non-Maxwellian distribution
function to Maxwellian (due to collisions) is studied. The initial
distribution function is a step-function in velocity space:

.. math::

   f(x,v,t) &= \frac{1}{2v_0} \quad &|v| < v_0 \\
            &= 0 \quad &|v| > v_0

where :math:`v_0 = 3 v_{th}/2`. The simulation is run to steady-state,
and the resulting Maxwellian compared with the exact solution. Note
that as the BGK operator (as all collision operators) conserves
density, momentum and energy, we can easily calculate the expected
solution (for comparison) from the parameters of the selected initial
distribution function. The results are shown in the figure below. As
can be seen, the steady-state solution matches the exact solution very
well.

.. figure:: s1-max-relax-cmp.png
  :width: 100%
  :align: center

  Relaxation of an initial step-function distribution function
  (red-line) to a Maxwellian. Black line is the numerical solution,
  while blue dots are the exact solution computed from the moments of
  the initial condition. See :doc:`s1
  <../../sims-2/boltz-bgk/s1/s1-bgk-boltz>` for input file.

