JE1: Solving Poisson Equation on 2D Periodic Domain
===================================================

The problem and solution technique
----------------------------------

With periodic boundary conditions, the Poisson equation in 2D

.. math::
  :label: eq:poisson-eqn

  \frac{\partial^2 \psi}{\partial x^2} + 
  \frac{\partial^2 \psi}{\partial y^2} = s(x,y)

can be solved using discrete Fourier transforms. We assume that the
domain is :math:`\Omega = [-L_x/2 \times L_x/2] \times [-L_y/2 \times
L_y/2]`, discretized using :math:`N_x \times N_y` cells. Integrating
:eq:`eq:poisson-eqn` over :math:`\Omega` and using periodicity shows
that the source :math:`s(x,y)` must satisfy

.. math::

  \int_\Omega s(x,y) dx dy = 0

In the solver implemented in Lucee the source is modified by
subtracting the integrated source from the RHS of :eq:`eq:poisson-eqn`
to ensure that this condition is met.

Once the source term is adjusted, its 2D Fourier transform is
computed. The Fourier transform of the solution is then

.. math::
  :label: eq:poisson-eqn

  \overline{\psi}(k_x, k_y) = \frac{\overline{s}(k_x,k_y)}{k_x^2+k_y^2}

where overbars indicate 2D Fourier transforms and :math:`k_x` and
:math:`k_y` are the wave-numbers. The zero wave-numbers are replaced
by :math:`10^{-8}` to prevent divide-by-zero. The `FFTW
<http://fftw.org/>`_ library is used to compute the transforms.

The algorithm is implemented in the class
`Lucee::PeriodicPoisson2DUpdater` class in the `proto` directory. 

This updater is for use in testing finite-volume/finite-difference 2D
incompressible flow algorithms. In this entry the stand-alone updater
is tested on and verified a few test problems.

Test problems
-------------

The domain is selected to be a square with :math:`L_x=L_y=2`. For each
problem second-order finite-differences are used to verify that the
solution satisfies the Poisson equation.

The first problem is for a Gaussian source term of the form

.. math::

  s(x,y) = e^{-10(x^2+y^2)}

The numerical solution is shown below.