The Fokker-Planck collision operator in Gkeyll
----------------------------------------------

`PDF of note <./_static/files/fpo-notes.pdf>`_.

Also see the note:
`Dimensionally independent indexing <./_static/files/dim-independent-idx.pdf>`_

There are several texts in which the Fokker-Planck operator (FPO) is
described and studied in detail. See, for example, the original paper
of `Rosenbluth, MacDonald and Judd from 1957
<./_static/files/RMD-1957-FPO.pdf>`_ (RMJ) or Chapter 3 of
"Collisional Transport in Magnetized Plasmas" by Helander and
Dieter. For the list of equation `see this page
<./_static/files/Formulary-FPO.pdf>`_ from the NRL Plasma Formulary.

Historically, the first derivation of the FPO in the case of Coulomb
potential (inverse square law) was by `Lev Landau in 1936
<./_static/files/Landau-FPO-1936.pdf>`_. However, the 1957 paper by
RMJ does not mention Landau's work at all.  As is usually the case,
the original papers by these Masters of the field remain highly
readable and still provide the best derivations.

My goal in this note is to list the FPO equations as implemented in
Gkeyll, with especial emphasis on properties needed to discretize them
using finite-volume and discontinuous Galerkin schemes. Note that we
use the *Rosenbluth form* of the equations and not the Landau
form. The FPO is a very complicated equation: it is a nonlinear
integro-differential equation in 3D velocity space and has a rich
structure, specially when combined with the particle motion in
self-consistent electromagnetic fields. Designing a general production
solver for the case of multi-species plasmas poses a formidable
challenge. (This is an understatement). In fact, I am not aware of any
widely usable production implementation of the operator for general,
multi-species problems. (Though there are some excellent
special-purpose solvers that are used in fusion physics).
