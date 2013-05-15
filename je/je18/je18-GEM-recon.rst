:Author: Ammar Hakim
:Date: May 15th 2013
:Completed: 
:Last Updated:

JE18: Five-moment two-fluid reconnection with zero guide field
==============================================================

.. note::

  These are *very* preliminary notes (as of May 15th 2013) on using
  five-moment two-fluid model to do reconnection with zero guide
  field. I have done many more simulations that the one showed here,
  and will put those up soon.

In this initial comparison I have used the same parameters as in the
original GEM problem, which are slightly different than in the PIC
simulations [Daughton2006]_, which use :math:`n_b/n_0=0.3`, while I
use :math:`n_b/n_0=0.2`. The BCs in my fluid simulations are open, and
the domain is :math:`25d_i \times 25d_i` on a grid of :math:`768\times
768` cells. This gives about 30 cells per :math:`d_i` and 6 per
:math:`d_e`.

For a complete description of the simulation see the Lua program
[:doc:`s238 <../../sims/s238/s238-gemguide-5m>`].

I do not completely understand the PIC open BCs. It seems like there
is a constant background density in the ghost cells in the PIC
simulation that serves as a source of particles. This can potentially
be done in a fluid simulation also, but I think there might be
reflections from the domain boundaries of plasma waves (though not of
EM waves).

These plots are very preliminary (as of 5/15/2013) as my simulation
was killed by the batch script as I "exceeded requested resources". I
have now restarted the simulation again.

The overall structures in the PIC and five-moment model look
similar. There is strong outflow of the electrons from the diffusion
region, the flow direction changing sign across the seperatrix.  In
the five-moment calculations this flow shear leads to Kelvin-Helmholtz
instabilities, which are visible early in time. These instabilities do
not appear in PIC, as far I can tell.

The other key difference between the PIC and five-moment results is in
the electron out-of-plane velocity. The five-moment results show a
much thinner and less extended flow structure than the PIC
results. However, there is a thin extended current sheet, although
much less extended than in the PIC simulations.

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
