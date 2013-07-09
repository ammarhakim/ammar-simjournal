:Author: Ammar Hakim
:Date: 7th July 2013
:Completed: 
:Last Updated:

JE20: Vlasov equation on bounded domain
=======================================

In many practical problems of interest the Vlasov equation needs to be
solved on a domain bounded by walls. In this document I describe the
boundary conditions for such bounded domain simulations and present
some example calculations with specified electric potential. For
periodic domain calculations, see :doc:`JE14
<../je14/je14-vlasov-fixed-pot.rst>` and :doc:`JE15
<../je14/je15-vlasov-poisson.rst>`.

Free-streaming
--------------

Consider the Vlasov equation with zero potential (free streaming of
particles)

.. math::

  \frac{\partial f}{\partial t} + v\frac{\partial f}{\partial x} = 0

on :math:`x\in[0,L]`, where :math:`f(x,v,t)` is the distribution
function. Let :math:`f(x,v,0)=f_0(v)` be the initial condition, and
assume that the particles streaming into the wall are completely
absorbed. The boundary condition can be hence written as

.. math::

  f(0,v,t) &= 0 \\
  f(L,-v,t) &= 0

for :math:`v>0`. In Gkeyll these boundary conditions are impemented by
simply setting the distribution function to zero in the ghost
cells. Proper upwinding insures that all the particles flow out of the
domain and none enter.

The total number of particles in the domain, :math:`n(t)`, can be
computed analytically as follows. There are :math:`f_0(v)dv` particles
with velocity :math:`v` and spread :math:`dv`. Hence, after time
:math:`t` the number of particles with that velocity will be
:math:`(L-vt)f_0(v)dv`. Also, after time :math:`t>L/v` all particles
with velocity :math:`v` will have been absorbed by the wall. Hence, we
can write

.. math::

  n(t) = \int_0^{L/t} (L-vt) f_0(v) dv 
       + \int_0^{L/t} (L-vt) f_0(-v) dv.

Here the first integral accounts for the particles flowing out of the
right wall, and the second the particles flowing out of the left wall.

To test the boundary conditions a simulation was run on
:math:`x\in[0,2\pi]` with :math:`f_0(v)` initialized to a Maxwellian

.. math::

  f_0(v) = \frac{n_0}{\sqrt{2\pi v_t^2}} e^{-v^2/2v_t^2}

with :math:`n_0=v_t=1`. A grid of :math:`16\times 32` was used with
piece-wise quadratic Serendipity basis functions. The number of
particles in the domain was compared to the exact solution. The
results are shown below.


.. Footnotes
.. ---------

.. .. [BC] For a more general situation, some particles
..   may be reflected back into the domain, often with a
..   time-delay. Assuming a linear response of the wall a more general
..   set of boundary conditions on the left wall can be written as

..   .. math::

..     f(0,v,t) &= \int_0^\infty dv' \int_0^t dt'\thinspace K(v,v') T(t-t')f(0,-v',t')

..   for :math:`v>0`, where :math:`K(v,v')` is a particle reflection
..   kernel and :math:`T(t-t')` is a time-delay kernel. Similar
..   expression can be written for the right wall.
