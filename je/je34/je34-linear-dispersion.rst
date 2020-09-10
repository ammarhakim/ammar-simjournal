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
equations. Other applications are to computing initial conditions that
excite specific modes, computing RF wave-propagation and potential
extension to retain quadratic nonlinearities to allow study of
wave-wave scattering/coupling physics.

Note on running the dispersion solver tool
------------------------------------------

To run the tool create an input file (Lua script, see examples below
in various figure captions) and run it as::

  gkyl multimomlinear -i inp.lua

This will create a output file with the eigenvalues stored in a Gkeyll
"DynVector" object. For each element in the dynvector, the first three
components are the components of the wave-vector and the rest the
corresponding :math:`\omega_n(\mathbf{k})` with real and imaginary
parts stored separately (next to each other). You can plot the real
part of the frequencies as function of wave-vector (say :math:`k_x`)
as::

  pgkyl -f inp_frequencies.bp val2coord -x0 -y 3::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_r$"

And the imaginary parts as::

  pgkyl -f inp_frequencies.bp val2coord -x0 -y 4::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_i$"

Often, it is useful to plot the eigenvalues in the complex plane (real
part on X-axis and imaginary part on the Y-axis). For this do::

  pgkyl -f inp_frequencies.bp val2coord -x3::2 -y 4::2 pl -s -f0 --no-legend --xlabel "$\omega_r$" --ylabel "$\omega_i$"  

Note that the frequencies are not outputed in any particular
order. Hence it is not possible to easily extract a single "branch" of
the dispersion relation from the output. Please see pgkyl help to
understand what the ``val2coord`` and ``pl`` (short for ``plot``) do
and how to use them.

The boolean flag ``calcEigenVectors`` can be set to ``true`` to
optionally compute the eigenvectors. These are stored in a
``CartField`` object.

Waves in a uniform plasma
-------------------------

In the first test I look at the problem of waves in a uniform
plasma. The first test assumes the plasma is composed of cold
electrons and ions (:math:`m_i/m_e = 25`) with the background magnetic
field :math:`\mathbf{B}_0 = (1.0, 0.0, 0.75)` transverse to the
wave-vector pointing in the X-direction . The following plot shows the
various plasma waves that appear in the system.

.. figure:: iso-cold-waves.png
  :width: 100%
  :align: center

  Waves in a cold magnetized plasma, with wave-vector transverse to
  the background magnetic field. Seen are the right (R) and left (L)
  polarized modes that asymptote to light waves at large
  :math:`k`. Also seen the second branch of the R mode which contains
  the Whister mode at low-:math:`k` and also the Alfven mode. See
  input file :doc:`iso-1-wave <iso-1-wave>`.

In the next problem I use the five-moment model with a finite pressure
:math:`p = 0.1` for both electrons and ions. The dispersion relation
is compared to the corresponding cold case in the following plot.

.. figure:: iso-5m-cmp-cold.png
  :width: 100%
  :align: center

  Waves in a cold (black) and five-moment finite-pressure (red)
  magnetized plasma, with wave-vector transverse to the background
  magnetic field. The finite pressure effect changes the sound-wave
  mediated branches, allowing them to propagating at high :math:`k`.
  See input file :doc:`5m-1-wave <5m-1-wave>`.

For the next test I use the ten-moment model with a diagonal pressure
tensor :math:`P_{xx} = P_{yy} = P_{zz} = 0.1` and all other components
set to zero. This gives the same scalar pressure as the previous
five-moment test. The ten-moment model has a a lot of modes. To
illustrate the differences between five-moments the following plot
shows the high-frequency branches of the dispersion relation. Note the
existence of *two* cyclotron harmonics marked with :math:`\omega_{ce}`
and :math:`2 \omega_{ce}` in the plot. These are missing in the
five-moment model. Also, the R- and L-mode structure is different,
with the lower cyclotron harmonic transitioning to the L-mode and the
upper cyclotron harmomic to the R-mode at larger :math:`k` values.

.. figure:: 10m-5m-cmp-elc.png
  :width: 100%
  :align: center

  Comparison of high-frequency branches of the dispersion relation for
  ten-moment (black) and five-moment (red) models. The ten-moment
  model contains the first two electron cyclotron harmonics which
  transition to the L- and R-mode at higher :math:`k`. See input file
  :doc:`10m-1-wave <10m-1-wave>`.

The following plot shows the low-frequency branches of the ten-moment
model compared to the five-moment model. Again, the two ion cyclotron
harmonics are seen as well as the various ion acoustic mediated
branches.

.. figure:: 10m-5m-cmp-ion.png
  :width: 100%
  :align: center

  Comparison of low-frequency branches of the dispersion relation for
  ten-moment (black) and five-moment (red) models. The first two ion
  cyclotron harmonics are seen. See input file :doc:`10m-1-wave
  <10m-1-wave>`.

The tests in these sections were done without any closures, and in
particular, without Landau closures. In the kinetic system many of the
modes seen above are damped.
  
References
----------

.. [Hakim2008] A.H. Hakim. "Extended MHD Modelling with the Ten-Moment
    Equations". *Journal of Fusion Energy*, **27** (1-2),
    36–43. http://doi.org/10.1007/s10894-007-9116-z

.. [Xie] H. Xie, Y. Xiao. "PDRK: A General Kinetic Dispersion Relation
    Solver for Magnetized Plasma". *Plasma Science and Technology*,
    **18** (2), 97–107. http://doi.org/10.1088/1009-0630/18/2/01
