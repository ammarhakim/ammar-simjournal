


Simulation Index
================

1D electrostatic shock problem.

.. list-table::
  :header-rows: 1
  :widths: 10,90

  * - Number
    - Description
  * - :doc:`s1 <s1/s1-es-shock>` 
    - :math:`Te/Ti=9.0`, Mach number 1.5, polyOrder 2 serendipity elements
  * - :doc:`s2 <s2/s2-5m-es-shock>` 
    - 5-moment version of s1, for comparison.
  * - :doc:`s3 <s3/s3-es-shock>` 
    - :math:`Te/Ti=1.0`, Mach number 0.5, polyOrder 2 serendipity elements
  * - :doc:`s4 <s4/s4-es-shock>` 
    - :math:`Te/Ti=1.0`, Mach number 3.0, m_i=40, polyOrder 2 serendipity elements
  * - :doc:`s5 <s5/s5-es-shock>` 
    - Same as s4, except :math:`Te/Ti=9.0`

Redoing the above in the proper, clean way.

.. list-table::
  :header-rows: 1
  :widths: 10,90

  * - Number
    - Description
  * - :doc:`m1 <m1/m1-es-shock>` 
    - :math:`Te/Ti=9.0`, Mach number 1.5, polyOrder 2 serendipity elements.
  * - :doc:`m2 <m2/m2-es-shock>` 
    - Same as m1, Mach number 2.0
  * - :doc:`m3 <m3/m3-es-shock>` 
    - Same as m1 Mach number 3.0
  * - :doc:`m4 <m4/m4-es-shock>` 
    - Same as m1, Mach number 5.0
  * - :doc:`m5 <m5/m5-es-shock>` 
    - Same as m1, Mach number 1.0
  * - :doc:`m6 <m6/m6-es-shock>` 
    - Same as m1, Mach number 0.75
  * - :doc:`m7 <m7/m7-es-shock>` 
    - Same as m1, Mach number 0.5
  * - :doc:`m8 <m8/m8-es-shock>` 
    - Same as m1, Mach number 0.25
  * - :doc:`m9 <m9/m9-es-shock>` 
    - Same as m1, Mach number 4.0
  * - :doc:`m10 <m10/m10-es-shock>` 
    - Same as m1, except run for longer, and 2X velocity mesh (to see if ions go unstable)
  * - :doc:`m11 <m11/m11-es-shock>` 
    - :math:`Te/Ti=100.0`, Mach number 1.5, polyOrder 2 serendipity elements.
  * - :doc:`m12 <m12/m12-es-shock>` 
    - Same as m1, except Mach number of 4, asymmetric beams, :math:`T_R/T_L = 4`.

Fluid sims corresponding to m series simulations [Source solver takes
about 1/2 the total run time! Also, Maxwell solver takes longer than
each fluid solver. This is very strange. Perhaps one can imagine a
hard-coded Maxwell solver (even Balsara based) and manually inverting
the source term. One way to speed up the fluid solver is to use ES
equations and not Maxwell equations at all.]

.. list-table::
  :header-rows: 1
  :widths: 10,90

  * - Number
    - Description
  * - :doc:`f1 <f1/f1-5m-es-shock>` 
    - Corresponds to m1
  * - :doc:`f2 <f2/f2-5m-es-shock>` 
    - Corresponds to m2
  * - :doc:`f3 <f3/f3-5m-es-shock>` 
    - Corresponds to m3
  * - :doc:`f4 <f4/f4-5m-es-shock>` 
    - Corresponds to m4
  * - :doc:`f5 <f5/f5-5m-es-shock>` 
    - Corresponds to m5
  * - :doc:`f6 <f6/f6-5m-es-shock>` 
    - Corresponds to m6
  * - :doc:`f7 <f7/f7-5m-es-shock>` 
    - Corresponds to m7
  * - :doc:`f8 <f8/f8-5m-es-shock>` 
    - Corresponds to m8
    