:Author: Ammar Hakim
:Date: May 15th 2013
:Completed: 
:Last Updated:

JE18: Five-moment two-fluid reconnection on open domains
========================================================

In this note I study the process of magentic reconnection with the
five-moment model. This model consists of a set of fluid equations for
each of the species in the plasma (electron and ions in this case),
coupled to the Maxwell equations via Lorentz force and currents. With
a scalar pressure, the fluid equations can be written in
non-conservative form as

.. math::

  \partial_t{n} + n \partial_j{u_j} + u_j \partial_j{n} &= 0 \\
  \partial_t{u_i}
  + u_k \partial_k{u_i}
  + \frac{1}{mn}\partial_i{p}
   &=
  \frac{q}{m}\left(E_i + \epsilon_{kmi}u_kB_m\right) \\
  \partial_t{p} + u_k\partial_k{p}
  &= -\gamma p \partial_k u_k

where :math:`n` is the number density, :math:`u_i` is the fluid
velocity, :math:`p` is the scalar pressure and :math:`E_i` and
:math:`B_i` are the electric and magnetic fields. The electromagnetic
fields are computed by solving the full Maxwell equations.

The simulations are performed on an open domain, with the plasma
initialized with a Harris current sheet with initial magentic field
given by

.. math::

  B_x(y) = B_0 \tanh{(y/L)}

where :math:`L` is the current sheet half-thickness. With this initial
magnetic field the conservation of total pressure (fluid plus
magnetic)

.. math::

  p_e(y) + p_i(y) + \frac{B(y)^2}{2\mu_0} = \mathrm{const.}

can be used to determine the number density and pressure profiles
(initial fluid temperatures are assumed constant). Using a background
number density :math:`n_b`, this gives

.. math::

  n(y) = n_0\mathrm{sech}^2{(y/L)} + n_b

The other parameters are taken from [Daughton2006]_ as

.. math::

  \frac{\rho_i}{L} = 1,\quad
  \frac{m_i}{m_e} = 25,\quad
  \frac{T_i}{T_e} = 5,\quad
  \frac{\omega_{pe}}{\Omega_{ce}} = 3,\quad
  \frac{n_b}{n_0} = 0.3,

where :math:`\rho_i=v_{thi}/\Omega_{ci}` is the ion gyroradius,
:math:`v_{thi}=\sqrt{2T_i/m_i}` is the ion thermal speed,
:math:`\Omega_{cs}=e B_0/m_s` is the species gyrofrequency and
:math:`\omega_{pe} = \sqrt{e^2n_0/\epsilon_0 m_e}` is the electron
plasma frequency.

Picking a normalization as :math:`n_0=\epsilon_0=\mu_0=m_i=1`, gives
:math:`B_0=1/15`, :math:`v_{the}=1/3\sqrt{6}\approx 0.361`, :math:`T_e
= B_0^2/12` and :math:`\rho_i=\sqrt{10/12}`,

The BCs are open (zero normal derivative of all quantities). The
domain is :math:`25d_i \times 25d_i`, where :math:`d_i=c/\omega_{pi}`
is the ion inertial length.

Different grid sizes were used: :math:`256\times 256`,
:math:`512\times 512` and :math:`768\times 768`. This gives about 10
(20, 30) cells per :math:`d_i` and 2 (4, 6) per :math:`d_e`,
respectively. For a complete description of the simulation see the Lua
program [:doc:`s238 <../../sims/s238/s238-gemguide-5m>`].

XXX

.. _fig:

  .. image:: s238-ne-diff.png
     :width: 100%
     :align: center

  .. image:: s238-bx-diff.png
     :width: 100%
     :align: center

  Number density (top) and magnetic field (bottom) along vertical
  slice at :math:`x=12.5d_i`. Inset in top plot shows number density
  in the middle of the slice, showing a small dip (probably numerical)
  also seen in Fig. 7 of the PIC paper. At the upstream edge of the
  diffusion region the magnetic field is :math:`B_x/B_0=0.81`. See
  [:doc:`s238 <../../sims/s238/s238-gemguide-5m>`].

.. _fig:

  .. image:: s238-ne.png
     :width: 100%
     :align: center

  .. image:: s238-uiz.png
     :width: 100%
     :align: center

  .. image:: s238-uix.png
     :width: 100%
     :align: center

  .. image:: s238-uey.png
     :width: 100%
     :align: center

  .. image:: s238-uex.png
     :width: 100%
     :align: center

  Number density, inflow ion velocity, outflow ion velocity,
  out-of-plane electron velocity and outflow electron velocity. Strong
  outflows are seen in the electron fluid with flow changing
  directions across the seperatrix. This leads to Kelvin-Helmholtz
  instabilities.
  
References
----------

.. [Daughton2006] William Daughton, Jack Scudder and Homa Karimabadi,
   "Fully kinetic simulations of undriven magnetic reconnection with
   open boundary conditions", *Physics of Plasmas*, **13**, 072101,
   2006.
