:Author: Ammar Hakim
:Date: July 15th 2014
:Completed: 
:Last Updated:

JE24: Tests for embedded boundary Euler solver
==============================================

.. contents::

.. warning::

  There is a subtle bug in the present implementation, which I have
  not been able to track down. Basically, when running in parallel, if
  the MPI domain boundary and the embedded boundaries are coincident
  or separated by one cell, the results differ very slightly from
  serial results. Note that for all other cases in parallel the bug
  does not manifest itself.

In this note I test the finite-volume embedded boundary updater on the
Euler equations. At present, Gkeyll uses a stair-stepped mesh to
represent the surface of an object. Although crude, this is sufficient
to tackle many problems, including magnetosphere modeling. We are
interested in the latter problem, and most global magnetosphere codes
use a stair-stepped mesh to represent the planet/moon surface. The
representation of the boundary may be improved in the future by using
either a conformal cut-cell mesh, or using a multi-block body-fitted
mesh. Note that for viscous flows, using stair-stepped or even
cut-cell BCs is not a good idea as near a wall, to capture the
boundary layer accurately, one needs meshes aligned with the surface.

For now, stair-stepped meshes are sufficient.

Note on geometry representation and embedded boundary conditions
----------------------------------------------------------------

The key step in doing embedded boundary simulations is to setup the
geometry. For this, one needs to create an Gkeyll field with a single
component. This field should store a positive number when the
corresponding point is inside the domain and a negative number when it
is outside the domain. Rather complex objects can be created by
combining a set of elementary shapes (circles, boxes, ellipses, ...)
using shift/rotate and logical operators. For example, to represent a
circle one can use the Lua function

.. code-block:: lua

  function inCircle(x,y, xc,yc,rad)
    return (x-xc)^2 + (y-yc)^2 < rad^2 and 1.0 or -1.0
  end

This function evaluates to 1.0 when inside the sphere, and to -1.0
otherwise. Note that for flow *around* a sphere one needs to flip the
signs in the above equation.

Given two in/out functions :math:`d_A(x,y)` and :math:`d_B(x,y)` can
can create a new functions representing the union, intersection and
subtraction as follows

.. math::

  d_U(x,y) &= \max(d_A(x,y), d_B(x,y)) \qquad &\mathrm{union} \\
  d_I(x,y) &= \min(d_A(x,y), d_B(x,y)) \qquad &\mathrm{intersection} \\
  d_S(x,y) &= \min(d_A(x,y), -d_B(x,y)) \qquad &\mathrm{substraction}.

In three dimensions several other operations are possible: rotation of
a 2D shape about some axis, extrusion of a 2D shape along a curve,
etc. Arbitrary shapes can be represented in this way, however, the
final function representing a complex shape may not be very readable.

Once the object geometry has been created, the boundary condition
update needs to be created (in 2D) using the `StairSteppedBc2D`
updater. This updater takes in the in/out field and a list of BCs to
apply. 

Note that due to the way in which BCs are applied this updater also
needs the direction in which to apply the BC (using a setDir
method). In essence, I am assuming that this updater is used as part
of a dimensionally-split algorithm in which the BC is applied in a
particular direction and then the update in that direction is
performed. For unsplit algorithms, embedded BCs require more
complicated data-structures.

.. note::

  In most of these simulations, I need to use HLLE fluxes (enabled
  using numericalFlux = "lax" and useIntermediateWave = true in
  HyperEquation.Euler block). The Roe fluxes causes severe problems,
  including carbuncle problem and launch of spurious shocks from
  stair-stepped cells. So, the solutions are perhaps rather diffusive.
