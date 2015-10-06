:Author: Ammar Hakim
:Date: 27th September 2015
:Completed: 
:Last Updated:

JE27: Neutral and multi-fluid simulations of Rayleigh-Taylor instability
========================================================================

.. contents::

In this note I study the Rayleigh-Taylor (RT) instability with neutral
and multi-fluid models. The goal is to understand the effects of Hall
and other physics not captured in ideal MHD models in 3D RT
instability. In 3D, in a :math:`beta=1` plasma one can imagine that
there will be formation of current sheaths and fast reconnection,
leading to new and novel insights into turbulent reconnection and
other physics [1]_.

I am using the dimensionally split, positivity preserving
wave-propagation algorithm. At some point we may redo these
simulations using a DG scheme, but at present the WAVE algorithm
appears to be very robust and efficient, potentially allowing good
results to be obtained in reasonable amount of resources.

The algorithm is described in :doc:`JE22 <../je22/je22-euler-2d>`. A
RT simulation benchmark, comparing to Liska and Wendroff [Liska2003]_
was performed. This initial benchmark showed that the Gkeyll algorithms
compare very well with other schemes.

Footnotes
---------

.. [1] Apparently, the case of :math:`\beta=1` plasmas (known as
   "Parker instability") has been studied extensively in the
   astrophysical literature. See for example, J. Kim,
   S.S. Hong, D. Ryu, and T.W. Jones, "Three-dimensional evolution of
   the Parker instability under uniform gravity", *Astrophys. J*.,
   **506**, L139 (1998). However, note that these are ideal MHD
   simulations, and do not capture the correct physics needed for a
   proper understanding of the current sheath dynamics and fast
   magnetic reconnection.

References
----------

.. [Liska2003] Liska, R., & Wendroff, B. "Comparison of Several
   Difference Schemes on 1D and 2D Test Problems for the Euler
   Equations", *SIAM Journal on Scientific Computing*, **25** (3),
   995–1017. doi:10.1137/S1064827502402120


