:Author: Ammar Hakim
:Date: 22nd May 2016
:Completed: 
:Last Updated:

JE30: Computing moments of a distribution function
==================================================

.. contents::

Introduction
------------

In this note I test the updater that computes moments of a
distribution function. The updater ``DistFuncMomentCalcCDIMFromVDIM``
computes number density, momentum density, total particle energy and
pressure tensor. For coupling to field equations, the number density
and momentum density are needed, but the other moments are useful for
diagnostics. In the near future we will also extend the updater to
compute all independent components of the heat-flux tensor.

Computing moments is tricky in DG schemes as the various convolutions
in each cell need to be done with care. Further, these then need to be
summed across over all velocity space cells (at a given configuration
space cell). In parallel, an "all-reduce" operation is
required. Gkeyll allows decomposition in phase-space (and not just
configuration space), making the parallel code rather
complicated. Hence, careful tests are needed to ensure that the moment
calculator works correctly.

The updater computes the following moments

.. math::

   n &= \int_{-\infty}^{\infty} f(\mathbf{v}) d\mathbf{v} \\
   nu_i &= \int_{-\infty}^{\infty} v_j f(\mathbf{v}) d\mathbf{v} \\
   \mathcal{P}_{ij} &= \int_{-\infty}^{\infty} v_i v_j f(\mathbf{v}) d\mathbf{v} \\
   \mathcal{E} &= \frac{1}{2} \int_{-\infty}^{\infty} v^2
   f(\mathbf{v}) d\mathbf{v}

Note that the moment calculator computes the pressure tensor in the
*lab frame* and hence contains the contribution from the Reynolds
stresses:

.. math::

   \mathcal{P}_{ij} = P_{ij} + mn u_i u_j

where

.. math::

   {P}_{ij} = \int_{-\infty}^{\infty} (v_i-u_i) (v_j-u_j)
   f(\mathbf{v}) d\mathbf{v}

is the pressure tensor in the fluid frame.   


A Multi-modal Gaussian
----------------------

To test the moment calculator, I initialize the distribution function
with a multi-modal Gaussian. This allows the pressure tensor to be
anisotropic, with each of the six components specified
arbitrarily. The multi-modal Gaussian in :math:`N` velocity space
dimensions is given by

.. math::

   \mathcal{G}_N =
   \frac{n}{(2\pi)^{N/2}\triangle^{1/2}}\exp(-\frac{1}{2}\Theta^{-1}_{ij}c_ic_j)

where :math:`\triangle = \mathrm{det}(\Theta_{ij})`,
:math:`\Theta_{ij} = P_{ij}/mn` and :math:`c_i = v_i-u_i`.
