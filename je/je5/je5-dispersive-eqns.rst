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

Linearization of a hyperbolic balance law
-----------------------------------------

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

