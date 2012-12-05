:Author: Ammar Hakim
:Date: December 4th
:Completed: 
:Last Updated:

JE16: Auxiliary equations and tests of local DG scheme for advection-diffusion equations
========================================================================================

.. contents::

The local DG scheme can be used to solve equations with diffusion
terms. This is achieved by replacing the second order terms in the
original system with first order terms. This leads to a system of
first order equations which can then be updated with the standard DG
algorithm.

Gkeyll has an arbitrary dimensional nodal DG updater that can be used
to solve such equations. The basic idea is to pass extra fields as
auxiliary variables. This allows a very general scheme in which
several different types of equations can be solved. In this document I
test the nodal DG updater and use the auxiliary variables framework to
solve a number of problems, including those with diffusion.

Auxiliary variables for solving 2D advection equations
------------------------------------------------------

In the journal entry :doc:`JE12 <../je12/je12-poisson-bracket>` the
Poisson bracket updater is tested with passive advection in a 2D flow
field given by :math:`u_x = \partial \psi/ \partial y` and :math:`u_y
= -\partial \psi/ \partial x`, where :math:`\psi` is a specified
potential. The same problems can also be solved with the nodal DG
updater by specifying the flow velocity as an auxiliary variable. This
allows tracing passive advection in a general flow field (i.e. not
constrained to be incompressible or 2D). In this set of tests this
ability is exercised by solving the 2D scalar advection equation

.. math::

  \frac{\partial f}{\partial t} + \nabla\cdot (\mathbf{u}f) = 0

with specified flow field :math:`\mathbf{u}(x,y,t)`.

Rigid-body rotating flow
++++++++++++++++++++++++

In this test a rigid body rotating flow is initialized by initializing
the flow field as

.. math::

  u_x(x,y) &= -y+1/2 \\
  u_y(x,y) &= x-1/2

which represents a counter-clockwise rigid body rotation about
:math:`(x_c,y_c)=(1/2,1/2)` with period :math:`2\pi`. Hence,
structures in :math:`f` will perform a circular motion about
:math:`(x_c,y_c)`, returning to their original position at
:math:`t=2\pi`.

The simulation was performed with  with an initial cosine hump of the
form

.. math::

  f(x,y,0) = \frac{1}{4}
  \left[
    1 + \cos(\pi r)
  \right]

where

.. math::

  r(x,y) = \min(\sqrt{(x-x_0)^2 + (y-y_0)^2}, r_0)/r_0

For this problem, :math:`r_0=0.2` and :math:`(x_0,y_0) = (1/4,
1/2)`. To test convergence, the simulation was run to :math:`t=2\pi`
with grid sizes :math:`16\times 16`, :math:`32\times 32` and
:math:`64\times 64` grid (keeping time-step fixed). The results were
compared to the initial condition and errors computed and are shown in
the following table.

Next, a third order spatial scheme was used to compute
the solution to :math:`t=4\pi` at which point the cosine hump has
advected twice about the origin. The figure below shows the solution
at four different times, indicating that the algorithm essentially
advects the initial hump without any significant distortion.
