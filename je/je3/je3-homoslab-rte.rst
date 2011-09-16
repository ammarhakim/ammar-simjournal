:Author: Ammar Hakim
:Date: September 15th 2011

JE3: Testing the radiation transport equation solver in a homogeneous slab
==========================================================================

Problem formulation and Fourier decomposition
---------------------------------------------

In this entry I test the radiation transport equation (RTE) solver
implemented in Lucee, designed to solve the RTE in a homogeneous
slab. This solver can be be used as a building block for solving the
RTE in an inhomogeneous slab. The RTE is a linear integro-differential
equation and is written as

.. math::

  \mu\frac{\partial L(\tau,\mu,\phi)}{\partial \tau} + L(\tau,\mu,\phi)
  =
  \frac{\varpi}{4\pi}
  \int_{-1}^1 \int_0^{2\pi}
  p(\cos\Theta) L(\tau,\mu,\phi) d\mu d\phi

where

- :math:`L(\tau,\mu,\phi)` is the radiance in units of Watt
  :math:`\mathrm{m}^{-2}` :math:`\mathrm{sr}^{-1}`
  :math:`\mathrm{nm}^{-1}`,

- :math:`\tau` is the optical depth,

- :math:`\varpi` is the albedo of single scattering,

- :math:`\mu` is the cosine of the polar angle measured with the
  positive :math:`Z`-axis and :math:`\phi` is the azimuthal angle,

- :math:`p(\cos\Theta) = \sum_{l=0}^L\beta_lP_l(\cos\Theta)` is the
  phase function, where :math:`\Theta` is the scattering angle and
  :math:`\beta_0=1`.