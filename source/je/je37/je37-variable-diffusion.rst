:Author: Ammar Hakim
:Date: September 8th 2025
:Completed: 
:Last Updated:

JE37: Recovery Schemes for Variable Coefficient Diffusion
=========================================================

Consider computing the weak-form of the term

.. math::
   
   g = \frac{\partial}{\partial x}\left( D \frac{\partial f}{\partial x} \right)

Here :math:`D = D(x)` is the position dependent diffusion
coefficient. This term appears, for example, in diffusion equation and
Poisson equation with variable coefficient. We will at first assume
:math:`D` is continuous and smooth. Of course, this *does not mean*
that it's DG representation is continuous. 

Multiply :math:`g` by a weight :math:`w` and integrate by parts *once*
over a single cell :math:`I_j = [x_{j+1/2}, x_{j-1/2}]` to get the
weak-form:

.. math::
   
   \int_{I_j} w g \thinspace dx
   =
   \left. \left(w D \frac{\partial f}{\partial x} \right) \right\rvert_{j-1/2}^{j+1/2}
   -
   \int_{I_j} \frac{dw}{dx} D  \frac{\partial f}{\partial x} \thinspace dx.

From this weak-form we see that we have to compute interface values of
:math:`f`, :math:`D` and :math:`f_x`. To do this we will recover
:math:`D` across the interfaces :math:`x_{j\pm 1/2}`. We then have one
of two choices for :math:`f`:

#. Recover :math:`f` across interfaces :math:`x_{j\pm 1/2}` for use in
   the surface terms, but use the DG representation of :math:`f` in
   the volume term

#. Recover :math:`f` across three cells :math:`[x_{j-1}, x_j,
   x_{j+1}]` and use this in both the interface values and the volume
   term.

Note that in each of the above cases we will use the DG representation
of :math:`D` for the volume term.
