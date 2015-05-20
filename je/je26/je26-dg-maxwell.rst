:Author: Ammar Hakim
:Date: May 19th 2015
:Completed: 
:Last Updated:

JE26: Benchmarking a discontinuous Galerkin algorithm for Maxwell equations
===========================================================================

.. contents::

In this note I benchmark a nodal discontinuous Galerkin solver for
Maxwell equations. This solver will be used as part of the
Vlasov-Maxwell and multi-fluid equations. For notes on Maxwell
equations and divergence cleaning, please see :doc:`JE6
<../je6je6-maxwell-solvers>`.

For Maxwell equations, it can be shown that the *spatial operator* of
the DG scheme conserves the electromagneic (EM) energy exactly when
using *central fluxes*. However, with *upwinding* the EM energy
*decays monotonically*. Note that the time-stepping scheme (in this
case a SSP-RK3 scheme) will add some dissapation, and hence the fully
discrete scheme will decay the energy (not necessarily monotonically),
independent of numerical fluxes used.

