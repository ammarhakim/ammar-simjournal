:Author: Ammar Hakim
:Date: 30nd Jan 2016
:Completed: 
:Last Updated:

JE28: Some tests for Boltzmann-BGK equations
============================================

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
higher dimensions also.


