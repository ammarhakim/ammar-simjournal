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
<../je6/je6-maxwell-solvers>`.

For Maxwell equations, it can be shown that the *spatial operator* of
the DG scheme conserves the electromagnetic (EM) energy exactly when
using *central fluxes*. However, with *upwinding* the EM energy
*decays monotonically*. Note that the time-stepping scheme (in this
case a SSP-RK3 scheme) will add some dissipation, and hence the fully
discrete scheme will decay the energy (not necessarily monotonically),
independent of numerical fluxes used.

Gkeyll at present (as of May 19th 2015) has polynomial order 1 and 2
serendipity elements, and arbitrary order Lagrange tensor elements. In
the near future we hope to have higher order serendipity element as
well as other basis sets in which the number of degrees-of-freedom
(DOFs) are minimized.

The nodal DG solver works in 1D, 2D and 3D, but only 1D and 2D studies
are shown below.

Convergence of 1D and 2D scheme
-------------------------------

To test the convergence of the scheme, we can use a general wave
solution on a periodic domain. Let :math:`\mathbf{k} = 2\pi(k,l,m)/L`
be a wave-vector such that :math:`|\mathbf{k}|\ne 0`. Let
:math:`\mathbf{n}` be a non-zero vector such that
:math:`\mathbf{k}\cdot\mathbf{n} = 0`. Then, a general periodic
solution of Maxwell equations on a domain :math:`L\times L\times L` is

.. math::

  \mathbf{E} &= E_0 \cos(\mathbf{k}\cdot\mathbf{x} - |\mathbf{k}|c t +
  \delta)\hat{\mathbf{n}} \\
  \mathbf{B} &= E_0\sqrt{\mu_0\epsilon_0} 
    \cos(\mathbf{k}\cdot\mathbf{x} - |\mathbf{k}|c t + \delta)\hat{\mathbf{k}}\times\hat{\mathbf{n}}

where :math:`\hat{\mathbf{n}}` and :math:`\hat{\mathbf{k}}` are unit
vectors. We can use this solution to generate 1D and 2D initial
conditions that can be used to test the convergence of the scheme.



