:Author: Ammar Hakim
:Date: September 23st 2011

JE5: Hyperbolic balance laws with dispersive source terms
========================================================

While solving the :doc:`two-fluid Riemann problems
<../je4/je4-twofluid-shock>` small oscillation are observed in the
solution. On first sight it might be thought that these oscillations
are numerical and do not reflect a correct solution to the ideal
two-fluid equation system. In this note I show that such solutions can
(and do) occur when certain types of source terms are present and the
homogeneous system is hyperbolic.

Linearization and hyperbolic balance laws with dispersive sources
-----------------------------------------------------------------

Linearizing a hyperbolic balance law can lead to important insights
into the structure of solutions, although linearization can not reveal
the rich non-linear phenomena (shocks, rarefactions and contact
discontinuities) that are described by these equations. We start from
the non-conservative form of the equations

.. math::

  \frac{\partial \mathbf{v}}{\partial t} 
  + \mathbf{A}_p\frac{\partial \mathbf{v}}{\partial x} = \mathbf{s}_p

where :math:`\mathbf{v}` are the primitive variables,
:math:`\mathbf{A}_p` is a Jacobian matrix and :math:`\mathbf{s}_p` are
source terms. Linearize about a uniform equilibrium
:math:`\mathbf{v}_0` to write :math:`\mathbf{v} = \mathbf{v}_0 +
\mathbf{v}_1`, where :math:`\mathbf{v}_1` is a small
perturbation. Using a Taylor series expansion to first-order to write
:math:`\mathbf{s}_p(\mathbf{v}) = \mathbf{s}_p(\mathbf{v}_0) + \left(
{\partial \mathbf{s}_p}/{\partial \mathbf{v}} \right)_{\mathbf{v}_0}
\mathbf{v}_1` and :math:`\mathbf{A}_p(\mathbf{v}) =
\mathbf{A}_p(\mathbf{v}_0) + \left( {\partial \mathbf{A}_p}/{\partial
\mathbf{v}} \right)_{\mathbf{v}_0} \mathbf{v}_1` and letting
:math:`\mathbf{M}_p \equiv {\partial \mathbf{s}_p}/{\partial
\mathbf{v}}`, the linear form of the non-conservative equation becomes

.. math::

  \frac{\partial \mathbf{v}_1}{\partial t} 
  + \mathbf{A}_p(\mathbf{v}_0)\frac{\partial \mathbf{v}_1}{\partial x} 
  = \mathbf{M}_p(\mathbf{v}_0)\mathbf{v}_1.

To understand the mathematical structure of the linear equation a
Fourier representation of the solution is assumed and each mode is
represented as :math:`\mathbf{v}_1 = \mathbf{\hat{v}}_1 e^{i\omega t}
e^{i k x}`, where :math:`\omega` is the frequency and :math:`k` is the
wavenumber. Using this we obtain

.. math::

  \left[
    i\omega\mathbf{I} + ik\mathbf{A}_p(\mathbf{v}_0) - \mathbf{M}_p(\mathbf{v}_0)
    \right] \mathbf{v}_1 = 0,

where :math:`\mathbf{I}` is a unit matrix. For non-trivial solutions
the determinant of the matrix in the square brackets must
vanish. Another way to state this condition is that the frequency and
wavenumbers must be related by the *dispersion relations*
:math:`\omega = \omega(k)` which are the eigenvalues of the matrix
:math:`-k\mathbf{A}_p(\mathbf{v}_0) -
i\mathbf{M}_p(\mathbf{v}_0)`. Thus if
:math:`\lambda^p(k,\mathbf{v}_0)` is the :math:`p^{\textrm{th}}`
eigenvalue of this matrix then the :math:`p^{\textrm{th}}` branch of
the dispersion relation is :math:`\omega = \lambda^p(k,\mathbf{v}_0)`.

If the dispersion relation is purely real then the equation system
will support non-decaying waves. Further, if the dispersion relation
is *linear*, then a wave packet will propagate without distortion. If,
however, if the dispersion relation is non-linear (but still real),
the wave packet will suffer dispersion, i.e. waves with different
wave-numbers will propagate with different group and phase speeds.

For the simple case of vanishing sources (:math:`\mathbf{s}_p=0`) the
dispersion relation reduces to :math:`\omega = \lambda^p k`, where
:math:`\lambda^p` are the eigenvalues of the Jacobian matrix
:math:`\mathbf{A}_p`. As the homogeneous equation is assumed to be
hyperbolic the eigenvalues are all real, indicating that for the
homogeneous case waves will propagate without dispersion or decay.

This simple analysis indicates that the linear solution will depend on
the nature of the eigenvalues of the source Jacobian matrix
:math:`\mathbf{M}_p`. In case the eigenvalues of this matrix are
*purely imaginary*, the dispersion relation will be real but waves
will sufferer dispersion as they propagate in the background uniform
equilibrium. In this case the system of equations is called
*hyperbolic balance laws with dispersive source terms*. (This is my
own terminology and I have not seen such equations discussed in the
literature).

The Dispersive Euler Equations
------------------------------

Several simple hyperbolic balance laws with dispersive source terms can
be constructed. However, a particularly useful system that has
properties similar to the ideal two-fluid equations is the following
*dispersive Euler* system

.. math::

  &\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho\mathbf{u}) = 0 \\
  &\frac{\partial \mathbf{u}}{\partial t} + 
  \mathbf{u}\cdot\nabla\mathbf{u} =
  -\nabla p/\rho + \lambda\mathbf{u}\times\mathbf{b} \\
  &\frac{\partial p}{\partial t} + \mathbf{u}\cdot\nabla p = 
  -\gamma p \nabla\cdot\mathbf{u}

where :math:`\mathbf{b}(\mathbf{x})` is a time-independent vector
field and :math:`\lambda` is a constant.

Consider solving the linearized system in one-dimension. Linearizing
around stationary uniform initial state :math:`\rho = \rho_0`,
:math:`p = p_0` and assuming :math:`\mathbf{b} = (0,0,b_z)`, where
:math:`b_z` is constant leads to the linear system