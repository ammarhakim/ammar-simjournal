<%!
  import subprocess

  def makeRst(lua):
    subprocess.Popen(["python", "./code/mkluarst.py", "-l", lua+".lua"])
    
%>
# Make document link
<%def name="mdl(sim, luaNm)"><% makeRst(luaNm) %>:doc:`${sim} <${luaNm}>` </%def>

Simulation Index
================

Following is a list of simulation numbers with one-line descriptions.

.. list-table::
  :header-rows: 1
  :widths: 10,90

  * - Number
    - Description
  * - ${mdl("s1", "s1/s1-periodic-poisson")}
    - Poisson solve on a 2D periodic grid, with Gaussian source
  * - s2
    - Poisson solve on a 2D periodic grid, with anisotropic Gaussian source
  * - s3
    - Poisson solve on a 2D periodic grid, with source a sum of two Gaussians
  * - s4
    - Same as s3, but with different grid spacing in X and Y.
  * - s5
    - 1D Euler Sod-shock with sonic point in rarefaction. Wave-propagation algorithm.
  * - s6
    - Exact solution to s5.
  * - s7
    - Same as s5, except using the MUSCL-Hancock scheme.
  * - s8
    - 1D Euler shock with low density/pressure region. Wave-propagation scheme.
  * - s9
    - Exact solution to s8.
  * - s10
    - Same as s8, except using the MUSCL-Hancock scheme.
  * - s11
    - 1D Noh problem. Wave-propagation algorithm.
  * - s12
    - Exact solution to s11.
  * - s13
    - Same as s11, except using the MUSCL-Hancock scheme.
  * - s14
    - 1D Euler shock with a stationary contact discontinuity. Wave-propagation scheme.
  * - s15
    - Exact solution to s14
  * - s16
    - Same as s14, except using the MUSCL-Hancock scheme.
  * - s17
    - 1D Euler shock with two strong shocks. Wave-propagation scheme.
  * - s18
    - Exact solution to s17
  * - s19
    - Same as s17, except using the MUSCL-Hancock scheme.
  * - s20
    - 1D Euler with stationary contact discontinuity. Wave-propagation scheme.
  * - s21
    - Exact solution to s20
  * - s22
    - Same as s20, except using the MUSCL-Hancock scheme.
  * - s23
    - 1D Euler with slowly moving contact discontinuity. Wave-propagation scheme.
  * - s24
    - Exact solution to s23
  * - s25
    - Same as s23, except using the MUSCL-Hancock scheme.
  * - s26
    - 1D Euler with sharp spike in density. Wave-propagation scheme.
  * - s27
    - Exact solution to s26
  * - s28
    - Same as s26, except using the MUSCL-Hancock scheme.
  * - s29
    - 1D Euler Woodward-Collela blast wave problem. Wave-propagation scheme.
  * - s30
    - Same as s29 run with higher-resolution to serve as an "exact" solution.
  * - s31
    - Same as s29, except using the MUSCL-Hancock scheme.
