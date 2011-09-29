:Author: Ammar Hakim
:Date: September 28st 2011

JE6: Solving Maxwell equations with wave-propagation and FDTD schemes
=====================================================================

.. contents::

In this note I compare solutions for Maxwell equations obtained with
the wave-propagation scheme and Finite-Difference Time Domain (FDTD)
scheme. These schemes actually solve different forms of the equations
and with different layout of the electromagnetic fields on the
grid. The wave-propagation scheme collocates all field components
while the FDTD scheme staggers the field components on the cell faces
and edges. The advantage of the collocated fields is that limiters can
be applied that allow the solver to be used when the fields have
discontinuities or sharp gradients. The staggered field formulation,
on the other hand, satisfies basic vector identities for the discrete
fields. This is specially important for the Maxwell equations in which
the divergence relations on electric and magnetic fields are not
explicitly used in the update equations. However, the staggered field
formulation is harder to extend to non-rectangular grids.

One of the intended applications of these Maxwell solvers is to the
solution of two-fluid equations. It is straightforward to apply the
wave-propagation Maxwell scheme to solve the two-fluid
equations. However, some work is required to make the FDTD scheme work
with the cell-centered fluid solvers, in particular, the staggered
fields need to be interpolated back to cell centers before computing
the Lorentz force.

Problem 1: Transverse Magnetic modes in a box
---------------------------------------------

In this problem the domain is a rectangular metal box :math:`[0, X]
\times [0, Y]`. The electric field is initialized with

.. math::

  E_z = E_0 \sin(ax) \sin(by)

where :math:`a = m\pi/X` and :math:`b = n\pi/Y`. All other field
components are set to zero. Also, :math:`E_0 = 1`, :math:`m=8`,
:math:`n=5`, :math:`X = 80` and :math:`Y=40`. All quantities are in SI
units. The simulation is evolved to 150 ns. The exact solution for
:math:`E_z` is

.. math::

  E_z = E_0 \sin(ax) \sin(by) \cos(\omega t)

where the frequency :math:`\omega = c \sqrt{a^2 + b^2}` and :math:`c`
is the speed of light.

The first set of simulations were performed on a grid of :math:`80
\times 40` cells. Results for the computed :math:`E_z` are compared to
the exact solution at 75 ns and 150 ns.