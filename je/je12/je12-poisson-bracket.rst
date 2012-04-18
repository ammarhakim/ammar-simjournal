:Author: Ammar Hakim
:Date: April 13th 2012
:Completed: 
:Last Updated:  

JE12: Studies with a discontinuous Galerkin Poisson bracket solver
==================================================================

.. contents::

In this entry I study a discontinuous Galerkin (DG) algorithm to
discretize the Poisson bracket operator. In particular, I look at the
spatial and temporal convergence properties of this algorithm and
study the ability of the algorithm to solve complicated problems.

The algorithm updates the equation

.. math::

  \frac{\partial \chi}{\partial t} + \{\chi,\psi\} = 0

the Poisson bracket :math:`\{\chi,\psi\}` is defined as

.. math::

  \{\chi,\psi\} = 
  \frac{\partial \chi}{\partial x}\frac{\partial \psi}{\partial y} 
  -
  \frac{\partial \chi}{\partial y}\frac{\partial \psi}{\partial x}.

This equation describes the advection of :math:`\chi` with the
advection velocity :math:`\mathbf{u} = \nabla\psi\times \mathbf{e}_z`
or :math:`u_x = \partial \psi/ \partial y` and :math:`u_y = -\partial
\psi/ \partial x`. Hence, although :math:`\chi` can be discontinuous,
:math:`\psi` must be be continuous.

Convergence of Poisson bracket algorithm
----------------------------------------

Temporal Convergence
++++++++++++++++++++

In the first set of tests the temporal convergence is tested. For this
a Gaussian pulse :math:`\chi(x,y,0) = \exp(-75(x-x_c)^2)` is
initialized, where :math:`x_c = 0.5` and :math:`x \in [0,1]` with
periodic boundary conditions. The potential is set to
:math:`\psi(x,y)=y` which advects the pulse with constant speed in the
X-direction. Simulations were run on a :math:`32\times 4` grid to
:math:`t=1.0`. The time-steps were adjusted to CFL numbers of 0.2,
0.1, 0.05 and 0.025. To isolate the errors from the temporal
discretization alone the differences in the solution were computed
between successive results and their convergence computed. 

The results with Runge-Kutta second-order schemes is presented in the
following table.

.. list-table:: Poisson bracket convergence for RK-2 time-stepping.
  :header-rows: 1
  :widths: 20,40,20,20

  * - CFL
    - Change in error
    - Order
    - Simulation
  * - :math:`0.1`
    - :math:`4.9346\times 10^{-3}`
    - 
    - :doc:`s105 <../../sims/s105/s105-pb-advection-1d>`
  * - :math:`0.05`
    - :math:`1.2315\times 10^{-3}`
    - 2.0
    - :doc:`s106 <../../sims/s106/s106-pb-advection-1d>`
  * - :math:`0.025`
    - :math:`3.0780\times 10^{-4}`
    - 2.0
    - :doc:`s107 <../../sims/s107/s107-pb-advection-1d>`

The results with Runge-Kutta third-order schemes is presented in the
following table.

.. list-table:: Poisson bracket convergence for RK-3 time-stepping.
  :header-rows: 1
  :widths: 20,40,20,20

  * - CFL
    - Change in error
    - Order
    - Simulation
  * - :math:`0.1`
    - :math:`1.9146\times 10^{-4}`
    - 
    - :doc:`s109 <../../sims/s109/s109-pb-advection-1d>`
  * - :math:`0.05`
    - :math:`2.4022\times 10^{-5}`
    - 2.99
    - :doc:`s110 <../../sims/s110/s110-pb-advection-1d>`
  * - :math:`0.025`
    - :math:`3.0023\times 10^{-6}`
    - 3.00
    - :doc:`s111 <../../sims/s111/s111-pb-advection-1d>`

Spatial Convergence
+++++++++++++++++++

To test the spatial convergence of the algorithms, a Gaussian pulse is
initialized and propogated diagonally across a unit square with
periodic boundary conditions. The pulse returns to its starting
position after unit time has elapsed. Note that diagonal propagation
is a harder problem than propagation parallel to grid lines: it not
only tests the isotropy of the scheme but also the ability of the
scheme to capture features propagating across grid lines.

The Gaussian pulse is

.. math::

  \chi(x,y,0) = \exp(-75 r^2)

where :math:`r = \sqrt{(x-x_c)^2+(y-y_c)^2}` and :math:`(x_c,y_c)` are
the corrdinates of the center of the pulse. The potential is selected
as

.. math::

  \psi(x,y) =y - x

giving an advection speed of :math:`\sqrt{2}` top right corner of the
domain. For all problems, an RK2 time-stepping scheme with fixed
time-step (for all spatial resolutions) was used.

In the first set of tests, the convergence of the second-order scheme
is tested. This scheme uses the second-order 4-node Lobatto
elements. Grids of :math:`32\times 32`, :math:`64\times 64` and
:math:`128\times 128` were used and convergence computed by comparing
to the initial conditions. Results are shown in the following table.

.. list-table:: Poisson bracket convergence for second-order spatial scheme
  :header-rows: 1
  :widths: 20,40,20,20

  * - Cell size
    - Average Error
    - Order
    - Simulation
  * - :math:`1/32`
    - :math:`1.4036 \times 10^{-3}`
    - 
    - :doc:`s112 <../../sims/s112/s112-pb-advection-2d>`
  * - :math:`1/64`
    - :math:`2.0966\times 10^{-4}`
    - 2.74
    - :doc:`s113 <../../sims/s113/s113-pb-advection-2d>`
  * - :math:`1/128`
    - :math:`4.6609\times 10^{-5}`
    - 2.17
    - :doc:`s114 <../../sims/s114/s114-pb-advection-2d>`

The solution computed on the :math:`32\times 32` grid is shown below.

.. figure:: s112-projected-solution.png
  :width: 100%
  :align: center

  Solution computed on a :math:`32\times 32` with the 2D Poisson
  bracket updater (left) with a slice in the X-direction (red, right)
  compared to exact solution (black) at :math:`t=0`. See :doc:`s112
  <../../sims/s112/s112-pb-advection-2d>` for input file.

In the second set of tests, the convergence of the third-order scheme
is tested. This scheme uses the third-order 8-node Serendipity
elements. Grids of :math:`12\times 12`, :math:`16\times 16`, and
:math:`32\times 32` were used and convergence computed by comparing to
the initial conditions. Results are shown in the following table.

