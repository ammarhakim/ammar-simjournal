:Author: Ammar Hakim
:Date: September 8th 2020
:Completed: 
:Last Updated:

JE34: Multi-moment multifluid linear dispersion solvers
=======================================================

.. contents::

Introduction
------------

In this note I benchmark and document the multi-moment multifluid
dispersion solver tool in Gkeyll. This solver allows arbitrary number
of species, each of which can be either an isothermal fluid, a
five-moment fluid or a ten-moment fluid. The fields can be computed
from Maxwell equations or Poisson equations, with the option of some
species "ignoring" the background fields. Certain forms of closures,
including non-local Hammett-Perkins Landau fluid closures can be used.

For the list of equations and a brief overview of the algorithm used,
please see `this technical note
<../../_static/files/gkyl-mom-lin.pdf>`_. Essentially, the key idea of
this algorithm is to convert the problem of finding the dispersion
relation to an *eigenvalue problem* and then use a standard linear
algebra package (`Eigen
<http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ in this case)
to compute the eigensystem. This allows great flexibility as there is
no need to directly find complex nonlinear polynomial roots or even
formulate the dispersion relation explicitly. I note that this
technique was described in my 2008 paper on the ten-moment model
[Hakim2008]_. Since then, it has been signifianctly developed and used
by Hausheng Xie and others. See, for example [Xie]_ and references
therein. As far as I know, however, the inclusion of the ten-moment
model and Landau closures is unique to Gkeyll.

Eventually, the goal is to perform a careful comparison of the linear
physics contained in the ten-moment system (especially with various
Landau and other simplified closures) with full kinetic
equations. Other applications are to computing RF wave-propagation and
potential extension to retain quadratic nonlinearities to allow study
of wave-wave coupling physics.

References
----------

.. [Hakim2008] A.H. Hakim. "Extended MHD Modelling with the Ten-Moment
    Equations". *Journal of Fusion Energy*, **27** (1-2),
    36–43. http://doi.org/10.1007/s10894-007-9116-z

.. [Xie] H. Xie, Y. Xiao. "PDRK: A General Kinetic Dispersion Relation
    Solver for Magnetized Plasma". *Plasma Science and Technology*,
    **18** (2), 97–107. http://doi.org/10.1088/1009-0630/18/2/01
