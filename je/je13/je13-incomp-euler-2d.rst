:Author: Ammar Hakim
:Date: May 15th 2012
:Completed: 
:Last Updated:  

JE13: 2D Incompressible Euler Solver
====================================

.. contents::

The 2D incompressible Euler equations can be written in
vorticity-streamfunction form

.. math::

  \frac{\partial \chi}{\partial t} + \nabla\cdot(\mathbf{u}\chi) = 0

where :math:`\chi` is the fluid vorticity and :math:`\mathbf{u} =
\nabla\psi\times\mathbf{e}_z` is the fluid velocity. Here :math:`\psi`
is the streamfunction determined from a Poisson equation

.. math::

  \nabla^2 \psi = -\chi. 

As the flow is incompressible (:math:`\nabla\cdot\mathbf{u}=0`) we can
rewrite the fluid equation in the form

.. math::

  \frac{\partial \chi}{\partial t} + \{\chi,\psi\} = 0

where :math:`\{\chi,\psi\}` is the Poisson bracket operator defined by

.. math::

  \{\chi,\psi\} = 
  \frac{\partial \chi}{\partial x}\frac{\partial \psi}{\partial y} -
  \frac{\partial \chi}{\partial y}\frac{\partial \psi}{\partial  x}.

In this entry I use the FE Poisson solver tested in :doc:`JE11
<../je11/je11-fem-poisson.rst>` and combine it with the Poisson
bracket algorithm tested in :doc:`JE12 <../je12/je12-poisson-bracket>`
to solve this set of equations.

Problem 1: A double shear flow
------------------------------

In this problem the simulation is initialized with two shear
layers. The initially planar shear layers are perturbed slightly due
to which they roll around each other, forming finer and finer
vortex-like features. The initial conditions for this problem are

.. math::
  \chi(x,y,0) = 
  \left\{
    \begin{array}{1 1}
      \delta\cos(x) - \frac{1}{\rho}\mathrm{sech}^2((y-\pi/2)/\rho) \quad y\le\pi \\
      \delta\cos(x) + \frac{1}{\rho}\mathrm{sech}^2((3\pi/2-y)/\rho) \quad y\gt\pi
    \end{array}
  \right.

For the results show below :math:`\rho = \pi/15` and :math:`delta =
0.05`.
