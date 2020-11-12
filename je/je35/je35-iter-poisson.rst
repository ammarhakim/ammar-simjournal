:Author: Ammar Hakim
:Date: November 11th 2020
:Completed: 
:Last Updated:

JE35: Benchmarking an iterative discontinuous Galerkin Poisson solver
=====================================================================

.. contents::

Introduction
------------

In this note I document and test an iterative solver for Poisson
equation:

.. math::

  \nabla^2 f = -s

The solver uses a recovery DG scheme (RDG). This details of the solver
will be presented later, however, the essential idea is to multiply
the equation by basis functions and perform two integration by
parts. Each surface term then is evaluated by recovering a continuous
function across the two cells, and taking the derivatives as
needed. The recovery process is subtle and needs care, specially in
higher dimensions. Also, boundary conditions need to be incorporated
in the recovery to obtain full accuracy. For now, I am only using
periodic BCs and will show results with other BCs later.

Gkeyll already has a direct solver for the RDG-based Poisson
equation. However, it is not parallel (and is difficult to
parallelize) and being a direct solver scales poorly with problem
size. Hence, in this note I explore two iterative solvers that in
their inner loop use the discrete RDG poisson operator to march a
pseudo-time-dependent problem to convergence. As the iterative solver
uses the *local* stencil to update the Laplacian it is trivial to
parallelize.

Super Time-Stepping Scheme
--------------------------

In the first scheme a diffusion equation is marched to steady-state:

.. math::

   \frac{\partial f}{\partial t} = \nabla^2 f + s

An explicit scheme would be very slow for this equation: the time-step
would be go as :math:`\Delta x^2` and so the number of iterations to
converge to steady-state would go as :math:`N^2`, where :math:`N` is
the number of cells in each direction. Hence, instead one can use the
Super-Time Stepping Scheme (STS) to accelerate the convergence. In
particular, I use the first-order (in time) Runge-Kutta-Legendre
(RKL1) scheme as described in [Meyer2014]_. In this scheme the
equations are marched in (pseudo) time using a :math:`s`-stage RK
scheme, where the number of stages are determined by how large a
time-step we wish to take. The RK scheme is given by

.. math::

   Y_0 &= f^n \\
   Y_1 &= Y_0 + \tilde{\mu}_1\Delta \tau L[Y_0] \\
   Y_j &= \mu_j Y_{j-1} + \nu_j Y_{j-2} + \tilde{\mu}_j \Delta \tau L[Y_{j-1}]

for :math:`2\le j \le s` and where

.. math::

   L[f] = \nabla^2_h f + s

where :math:`\nabla^2_h` is the *discrete* Laplacian operator computed
using the RDG scheme. The various coefficients are given by

.. math::

   \mu_{j} &=\frac{2 j-1}{j} ; \quad v_{j}=\frac{1-j}{j} \\
   \tilde{\mu}_{j} &=\frac{2 j-1}{j} \frac{2}{s^{2}+s}

If :math:`\Delta t_e` is the maximum time-step for an explicit scheme
and we wish to take a time-step that is :math:`W` times larger, then
the number of stages is computed from

.. math::

   s = \frac{1}{2} \lceil\sqrt{1+8W} - 1 \rceil + s_0

where :math:`s_0` are extra stages. Note that :math:`\Delta \tau = W
\Delta t_e`. Hence, in this scheme, the free parameters one must
choose are :math:`W` and :math:`s_0`. In Gkeyll one chooses this
scheme using the "RKL1" parameter in the `IterPoisson` updater. On top
of the "inner" STS stages one can also use a sequence acceleration
method that drives the system to steady-state even faster.

Second order Richardson Iteration
---------------------------------

In the second scheme we use a damped wave-equation to march to
steady-state:

.. math::

   \frac{\partial^2 f}{\partial t^2}
   + 2\nu \frac{\partial f}{\partial t}
   = 
   \nabla^2 f + s.

This scheme is a form of second order Richardson iterative method. See
[Golub1961]_ for general analysis of such schemes. For this specific
scheme see [Jardin2010]_ Chapter 3, Section 5.2. (Jardin uses
:math:`\tau = 1/\nu` instead, and calls this scheme "Dynamic
Relaxation"). The time-derivatives are replaced by simple central
differences to get the time-marching scheme:

.. math::

   \frac{f^{n+1} - 2f^n + f^{n-1}}{\Delta t^2}
   + \nu \frac{f^{n+1} - f^{n-1}}{\Delta t}
   = \nabla^2_h f^n + s.

The free parameter in this scheme is :math:`\nu`. To pick this, first
write the solution as :math:`f = f_0 + f_1`, where :math:`f_0` is the
steady-state solution and :math:`f_1` is the error. The error itself
satisfies the equation

.. math::

   \frac{\partial^2 f_1}{\partial t^2}
   + 2\nu \frac{\partial f_1}{\partial t}
   = 
   \nabla^2 f_1.

Now assume a mode :math:`f_1 \sim e^{-i\omega t}e^{i k
\mathbf{x}}`. Then we get the dispersion relation

.. math::

   \omega(\mathbf{k}) = -i\nu \pm \sqrt{ k^2 - \nu^2 }.

Hence, to ensure that all modes damp (and the second term remains
purely oscillatory) we choose

.. math::

   \nu = k_{min}

where :math:`k_{min}` is the smallest wavenumber that can be
represented on the grid. Typically, in 1D we have :math:`k_{min} =
2\pi/L`, where :math:`L` is the domain size. Note that the fastest
*frequency* we must resolve is approximately :math:`k_{max}` and the
time-step we must choose is :math:`k_{max} \Delta t \lt 2`. For simple
spectral approximation to the Laplacian we have :math:`k_{max} =
\pi/\Delta x`, or that the maximum stable time-step must be

.. math::

   \Delta t \lt \frac{2 \Delta x}{\pi}.

For central difference approximation :math:`\Delta t \lt \Delta
x/\sqrt{2}`.

In either case, as the (pseudo) time-step is *linearly* dependent on
the cell spacing, indicates that the scheme will converge *linearly*
with the number of cells in each direction. So, doubling the number of
cells in each direction in 3D will lead to twice as many
iterations. As there are 8 times more cells now, the scheme will hence
take 16 times longer to converge. This is scaling is dramatically
better than a direct solver, which would be :math:`8^3 = 512` times
more expensive due to the cost scaling of the LU decomposition.


Nishikawa's First Order Scheme
------------------------------

In [Nishikawa2007]_ studied a system of first-order relaxation
equations that reduce to the Poisson equation at steady-state:

.. math::

   \frac{\partial f}{\partial t} &= \alpha
   \left(
   \nabla\cdot\mathbf{g} + s
   \right) \\
   \frac{\partial \mathbf{g}}{\partial t} &= -\frac{1}{T_r}
   \left(
   \mathbf{g} - \nabla f
   \right)

where :math:`\alpha` and :math:`T_r` are parameters. In 3D, for
example, this is a system of 4 first-order equations. At steady-state
:math:`\mathbf{g} = \nabla f` and hence the system will converge to
the solution of the Poisson equation.

Now, write :math:`f = f_0 + f_1` and :math:`\mathbf{g} =
\mathbf{g}_0 + \mathbf{g}_1`, where :math:`f_0` and
:math:`\mathbf{g}_0` are steady-state solution. Then the errors
satisfy

.. math::

   \frac{\partial f_1}{\partial t} &= \alpha \nabla\cdot\mathbf{g}_1 \\
   \frac{\partial \mathbf{g}_1}{\partial t} &= -\frac{1}{T_r}
   \left(
   \mathbf{g}_1 - \nabla f_1
   \right).

Consider the 1D case and write this as

.. math::

   \frac{\partial }{\partial t}
   \left[
    \begin{matrix}
      f_1 \\
      g_x
    \end{matrix}
   \right]    
    +
    \left[
    \begin{matrix}
      0 & -\alpha \\
      -1/T_r & 0
    \end{matrix}    
   \right]
   \frac{\partial }{\partial x}
   \left[
    \begin{matrix}
      f_1 \\
      g_x
    \end{matrix}
   \right]
   =
   -\frac{1}{T_r}
   \left[
    \begin{matrix}
      0 \\
      g_x
    \end{matrix}
   \right].

As is easily seen, the eigenvalues of the Jacobian matrix are simply

.. math::

   \lambda_{1,2} = \pm \sqrt{\frac{\alpha}{T_r}}.

What this means is that the errors propagate at a finite speed and,
due to the relaxation term, damp away as they propagate.

Now, take the time-derivative of the first of these equations, use the
second equation and then the first equation to see that

.. math::

   \frac{\partial^2 f_1}{\partial t^2}
   + \frac{1}{T_r} \frac{\partial f_1}{\partial t}
   = 
   \frac{\alpha}{T_r}\nabla^2 f_1.

Hence, Nishikawa's scheme is identical the second order Richardson
iteration if we choose :math:`\alpha = T_r` and :math:`T_r =
1/2\nu`. Other choices are also possible, of course, and could lead to
iterative schemes with different properties.

As Nishikawa's scheme essentially reduces to solving a system of
hyperbolic (plus relaxation source) equations, the time-step for
stability will also be linearly proportional to :math:`\Delta x`, and
hence will have the same cost scaling as the two schemes described
above. In fact, for the choice :math:`\alpha = T_r` we will have
:math:`\lambda_{1,2} = \pm 1` and hence :math:`\Delta t = \Delta x`
(in 1D).

However, one serious disadvantage of this scheme is that it involves
solving *four* first-order equations in 3D, while the scheme in the
previous section has only a single second-order equation. The RDG
implementation for the second-order system in Gkeyll has the *same
cost* as the cost of a single first-order equation, and hence
Nishikawa's scheme will be approximately four times more expensive (in
3D) if the number of iterations are approximately the same.

Residual norm, updater structure
--------------------------------

To check convergence of the solution we use the *residual norm*
computed as

.. math::

   R_2[f,s] = \frac{\lVert \nabla_h f + s \rVert_2 }{\lVert s
   \rVert_2}

where :math:`\lVert \cdot \rVert_2` is the :math:`l_2`-norm of the
discrete solution. See `this note
<https://gkeyll.readthedocs.io/en/latest/dev/modalbasis.html#convolution-of-two-functions>`_
on how to compute :math:`l_2`-norm of the from the Gkeyll
representation of the DG solution.  For all tests below I use the
initial guess of zero, and hence the initial residual norm is
always 1. Typically, I set the condition of :math:`R_2 \lt 10^{-8}` as
the discretization error is typically larger than this. For some
:math:`p=2` tests with high resolution one needs a more stringent
error criteria.

An example of the use of the updater is below:

.. code:: lua

  local iterPoisson = Updater.IterPoisson {
     onGrid = grid,
     basis = basis,
     errEps = 1e-8, -- maximum residual error
     stepper = 'richard2',
     verbose = true,
  }
  iterPoisson:advance(0.0, {fIn}, {fOut})

Note the parameter `stepper` is set to "richard2" to select the second
order Richardson iteration scheme.  When the `verbose` flag is set the
updater will show messages on the console. You can also save the error
history by calling the `writeDiagnostics()` method after the updater
has converged:

.. code:: lua

  iterPoisson:writeDiagnostics()

This will produce a DynVector BP file which can be plotted in the
usual way. For example::

  pgkyl -f f1-r2-iter-periodic_errHist.bp pl --logy

Note that the `IterPoisson` updater is not really restricted to only
DG discretization of the Poisson equation. In fact, any equation
system and discretization can be used. For example, density weighted
diffusion or FEM discretization. The updater simply calls the
appropriate equation object to compute the residual and does not use
any equation or discritization specific information.
  
Convergence tests in 1D
-----------------------



    
References
----------

.. [Meyer2014] C.D. Meyer, D.S. Balsara, T.D. Aslam. "A stabilized
   Runge–Kutta–Legendre method for explicit super-time-stepping of
   parabolic and mixed equations". Journal of Computational Physics,
   **257** (PA), 594–626. http://doi.org/10.1016/j.jcp.2013.08.021,
   (2014).

.. [Golub1961] G.H Golub, R.S. Varga. "Chebyshev semi-iterative
   methods, successive overrelaxation iterative methods and second
   order Richardson iterative methods", Numerische Mathematik **3**,
   147–156. (1961)

.. [Jardin2010] S. Jardin. "Computational Methods in Plasma Physics",
   Chapman & Hall/CRC Computational Science Series (2010).

.. [Nishikawa2007] H. Nishikawa. "A first-order system approach for
   diffusion equation. I: Second-order residual-distribution
   schemes". Journal of Computational Physics, **227** (1),
   315–352. http://doi.org/10.1016/j.jcp.2007.07.029 (2007)
