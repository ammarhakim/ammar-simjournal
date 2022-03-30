:Author: Ammar Hakim
:Date: March 29th 2022
:Completed: 
:Last Updated:

JE36: Benchmarking general geometry Maxwell solver with annular (cylindrical) waveguides
========================================================================================

.. contents::

Introduction
------------

In this note I present a series of benchmark problems to ensure that
the Maxwell solver in GkeyllZero (G0) works on non-Cartesian
geometries. G0 allows an arbitrary mapping from computational to
physical space. Here I assume that the physical domain is a coaxial
cylinder with inner and outer radii :math:`r_0` and :math:`r_1` and
height :math:`L_z`. In this geometry I study two sets of problems:
transverse magnetic (TM) modes and transverse electric (TE) modes. For
a derivation of these solutions `see this document <../../_static/files/maxwell-cyl.pdf>`_.

Transverse Magnetic Modes
-------------------------

Transverse Magentic modes are modes in which the magnetic field is
perpendicular to the direction of propagation. In this case, the
electric field is in the :math:`Z`-direction and is given by

.. math::

   E_r(r,\phi,z,t) = E_z(r)e^{im\phi}e^{ik_nz}e^{-i\omega t}

where :math:`m` is the azimuthal mode number, :math:`k_n = 2\pi/L_z`
describes the axial mode, :math:`\omega` is the frequency. The radial
dependence of the field is given by

.. math::

   E_z(r) = a J_m\big[r\sqrt{\omega^2-k_n^2}\big] + b Y_m\big[r\sqrt{\omega^2-k_n^2}\big]

where :math:`J_m` and :math:`Y_m` are Bessel functions of the first
and second kind, and :math:`a` and :math:`b` are constants. The
boundary conditions :math:`E_z(r_0) = E_z(r_1) = 0` determine the
possible frequencies as roots of the equation

.. math::

  J_m\big[r_0\sqrt{\omega^2-k_n^2}\big]
  Y_m\big[r_1\sqrt{\omega^2-k_n^2}\big]
  -
  J_m\big[r_1\sqrt{\omega^2-k_n^2}\big]
  Y_m\big[r_0\sqrt{\omega^2-k_n^2}\big]
  = 0.   

Once we compute the desired root (for propagating modes we must have
:math:`\omega>k_n`) we can choose set :math:`a` to whatever we like
and then compute :math:`b` from

.. math::

   b = -a \frac{J_m\big[r_0\sqrt{\omega^2-k_n^2}\big]}{Y_m\big[r_0\sqrt{\omega^2-k_n^2}\big]}.

To compute the frequencies we can use a few lines of Maxmia. In
general, it is not possible to "guess" the interval in which the roots
lie and so it is best to first make a plot to determine the
interval. Then Maxima's root finder can compute the root to machine
precision. The function we need to find the root of can be written in
Maxima as

.. code-block::

   F(m,kn,r0,r1,w)  := bessel_j(m,r0*sqrt(w^2-kn^2))*bessel_y(m,r1*sqrt(w^2-kn^2))
       - bessel_j(m,r1*sqrt(w^2-kn^2))*bessel_y(m,r0*sqrt(w^2-kn^2))$   

Once this is define we can plot the specific case (say :math:`m=2`,
:math:`k_n=0`, :math:`r_0 = 2` and :math:`r_1 = 5`) we are interested
as

.. code-block::

   draw2d( grid=true, explicit( F(2,0,2,5,w), w, 0, 5) )$

This will show the following figure:

.. figure:: maxima-root-1.svg
  :width: 100%
  :align: center

  Plot of the nonlinear function whose roots (zero crossings) are the
  allowed frequencies. Maxima root-finder requires we find the
  interval in which the root is desired. We also need to ensure that
  the function changes sign only once in the interval.

Using this figure we can choose the interval :math:`[1,2]` and find
the root as

.. code-block::

   w1 : find_root( F(2,0,2,5,w), w, 1, 2 )$

This will yield :math:`1.19318673737701`. We can also find
higher-frequency roots by passing other intervals to the above
command. Once we have the frequency we can determine :math:`a` and
:math:`b` as described above, thus completing the solution.


2D Modes
++++++++

First consider the case in which :math:`k_n = 0`. This a 2D standing
mode inside an annular disk (i.e. there is no variation in the
$Z$-direction). We will choose :math:`r_0 = 2` and :math:`r_1 = 5` and
:math:`m=4`. For this the first two roots are :math:`\omega =
1.557919724821651` and :math:`\omega = 2.430327042902498`. The
following plot shows the solution at :math:`t=0` and :math:`t=T_0`,
where :math:`T_0 = 2\pi/omega` is the mode period. 
