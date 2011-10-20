:Author: Ammar Hakim
:Date: October 18th 2011
:Completed: 

JE8: Propagation into a plasma wave beach
=========================================

A plasma wave beach is a slab configuration in which the density
increases monotonically. An electromagnetic wave is pumped into the
beach (hence the name wave beach). The wave frequency :math:`\omega`
is arranged such that at some location :math:`x=x_c`, :math:`\omega =
\omega_p(x_c)`. At this location the electromagnetic wave suffers a
cutoff and is reflected back towards the drive plane, creating a
standing wave pattern.

In this entry this problem is simulated with Lucee. The plasma profile
is selected as

.. math::

  \omega_p(x,t) \delta t = \left(\frac{1-x}{L}\right)^5

where :math:`0<x<L` is the domain and :math:`\delta t` is a constant
with units of time. The electromagnetic wave is driven by a current
applied at the center of the last cell, i.e.,

.. math::

  J_y(x,t) = \delta(x-x_e) J_0\sin(\omega t)

where :math:`x_e = L-\Delta x /2`, where :math:`\Delta x` is the cell
size.

We pick :math:`\omega \delta t = \pi /10`. With :math:`L=1` the cutoff
location is :math:`x_c \approx 0.58`. We pick :math:`\delta t =
L/100c`, where :math:`c` is the speed of light. Further, the electron
temperature is set to 1 eV and the ion fluid is not evolved, i.e., on
the time-scale of the problem the ions are assumed to be immobile.

This problem is interesting as with conventional explicit FDTD methods
a numerical instability occurs at cutoff and resonance layers. This
means that implicit algorithms are needed to simulate such
problems. See the paper by David Smithe [Smithe2007]_ for an algorithm
in which the plasma is treated as a cold linear medium that is evolved
implicitly to avoid the instability.

EM wave propagation in vacuum
-----------------------------

In the first test the problem is simulated without the plasma. This is
essentially a test to ensure the EM wave propagation works correctly
with the current source. Simulations are performed by setting
:math:`J_0=1\times 10^{-12}` and run to :math:`t=5` ns. Grids with 100
and 200 cells are used. Results are shown below. Note that although
the solutions are smooth limiters need to be applied as the sudden
appearance of :math:`E_y` due to the current source causes a shock
that propagates into the ghost cell, which in turns spoils the
interior solution. The limiters are what cause the maxima in the
solutions shown below to get slightly flattened.

.. figure:: plasmabeach-maxwell-cmp.png
  :width: 100%
  :align: center

  Electromagnetic wave propagation in vacuum driven by a current
  source in the last cell. Shown here is the electric field
  :math:`E_y` at :math:`t=2.5` ns (top) and :math:`t=5.0` ns (bottom)
  for 100 cells [:doc:`s65 <../../sims/s65/s65-plasmabeach-maxwell>`]
  (red line) and 200 cells [:doc:`s66
  <../../sims/s66/s66-plasmabeach-maxwell>`] (black line). In the
  upper panel the electromagnetic wave has not yet propagated through
  the domain.


References
----------

.. [Smithe2007] David N Smithe, "Finite-difference time-domain
   simulation of fusion plasmas at radiofrequency time scales",
   *Physics of Plasmas*, **14**, Pg. 056104 (2007).




