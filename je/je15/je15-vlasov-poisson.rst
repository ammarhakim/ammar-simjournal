:Author: Ammar Hakim
:Date: July 16th
:Completed:  
:Last Updated:  

JE14: Studies with a DG electrostatic Vlasov solver
===================================================

.. contents::

In this document I study various standard "text-book" problems with a
discontinuous Galerkin algorithm for the solution of the electrostatic
Vlasov system. This system is written in Hamiltonian form as 

.. math::

  \frac{\partial f}{\partial t} = \{H,f\}

where :math:`f(x,v,t)` is the distribution function, :math:`H(x,v)` is
a Hamiltonian function and where :math:`\{H,f\}` is the Poisson
bracket operator defined by

.. math::

  \{H,f\} = 
  \frac{\partial H}{\partial x}\frac{\partial f}{\partial v} -
  \frac{\partial H}{\partial v}\frac{\partial f}{\partial x}.

The Hamiltonian takes the form

.. math::

  H = \frac{1}{2}v^2 + \frac{q}{m}\phi

where :math:`\phi(x,t)` is a scalar potential and :math:`q` and
:math:`m` are the particle charge and mass respectively. 

Two classes of problems will be considered: the first, in which the
potential is determined by solving a Poisson equation

.. math::

  \frac{\partial^2 \phi}{\partial x^2} = -\frac{\rho_c}{\epsilon_0}

where, :math:`\rho_c` is the total charge density. In this class of
problems the ion distribution function is not evolved: the ions simply
provide a neutralizing background to the electron motion. Hence, the
charge density is determined from

.. math::

 \rho_c = |e| \left(Zn_{io}(x) - \int_{-\infty}^{\infty} f(x,v,t) dv\right)

where :math:`n_{io}(x)` is the fixed ion number-density profile,
:math:`e` is the charge on an electron and :math:`Z` the ion charge
number.

In the second class of problems the electrons are assumed to be a
massless, isothermal fluid. In this case, the electron motion is
determined from the condition of quasi-neutrality that, for small
deviations, can be written as

.. math::

  n_{i}(x) = n_{eo}\left(1 + \frac{|e|\phi}{T_e}\right)

where :math:`n_{eo}` is the constant electron initial density and
:math:`T_e` is the fixed electron temperature. This allows the
determination of the potential once the ion number density is known.

Notice that in the first class of problems the distribution function
for the *electrons* is evolved, while in the second class of problems
the distribution function for the *ions* is evolved. Hence, they
represent the high-frequency (plasma frequency) and low-frequency (ion
sound speed) branch of the electrostatic plasma dispersion relation
respectively.

Tests for the Vlasov-Poisson system
-----------------------------------

Linear Landau damping
+++++++++++++++++++++

In this test the ability of the algorithm to capture the phenomena of
Landau damping is shown. For this the electrons are initialized with a
perturbed Maxwellian given by

.. math::

    f(x,v,0) = \frac{1}{\sqrt{2\pi v_t}} \exp(-v^2/2v_t^2)
    (1+\alpha\cos(kx))

where :math:`v_t` is the thermal velocity, :math:`k` is the
wave-number and :math:`\alpha` controls the perturbation. The
computational domain is :math:`(x,v) \in [-\pi/k,\pi/k] \times
[-6,6]`. This means the domain contains a single wave. Periodic
boundary conditions imposed in the spatial direction and open
boundary conditions in the velocity direction. The ion density is
simply set to :math:`n_i(x) = 1`. For all these tests we set
:math:`m=1`, :math:`\epsilon_0=1`, :math:`Z=1`, :math:`k=1/2` and
:math:`\alpha=0.01`. With these settings, the plasma frequency is
:math:`\omega_{pe}=1` and Debye length is :math:`\lambda_D =
\sqrt{T_e}`.

To test the dependence of the damping rate on normalized wavenumber
(:math:`k\lambda_D`) simulations were performed with different
:math:`T_e`. The field energy was computed and recorded as a function
of time. An example plot of field energy is shown in the following
figure for the case :math:`T_e=1.0`.

.. figure:: s151-field-energy.png
  :width: 100%
  :align: center

  Field energy (blue) as a function of time for linear Landau damping
  problem with :math:`k = 0.5` and :math:`T_e = 1.0`. The red dots
  represent the maxima in the field energy which are used to compute a
  linear least-square fit. The slope of the black line gives the
  damping rate. See [:doc:`s151
  <../../sims/s151/s151-landau-damping-vp>`] for the input file.

To understand the dependence of the damping rate on Debye length, the
electron temperature was varied. The results are shown in the
following figure.

.. figure:: damping-rates-elc-osc.png
  :width: 100%
  :align: center

  Damping rate from Landau damping for electron plasma oscillations as
  a function of normalized Debye length. The damping rate increases
  rapidly with electron temperature (Debye length), eventually
  saturating.

Nonlinear Landau damping
------------------------

For this problem the we set :math:`\alpha = 0.5`, rapidly driving the
system nonlinear. Other parameters are the same as for the linear
Landau damping problem with :math:`k=0.5` and :math:`T_e=1.0`. The
field energy history and distribution function at various times are
shown in the following figures. Full details of the evolution of the
distribution function can be seen in `this movie
<../../_static/s162-distf.mov>`_.

.. figure:: s162-field-energy.png
  :width: 100%
  :align: center

  Field energy as a function of time for nonlinear Landau damping
  problem with :math:`k = 0.5`, :math:`T_e = 1.0` and
  :math:`\alpha=0.5`. The initial perturbation decays at a rate of
  :math:`\gamma = -0.2916`, after which the damping is halted from
  particle trapping. The growth rate of this phase is :math:`\gamma =
  0.0879`. See [:doc:`s162 <../../sims/s162/s162-landau-damping-vp>`]
  for the input file.

.. figure:: s162-landau-damping-vp_distf.png
  :width: 100%
  :align: center

  Distribution function at different times for the nonlinear Landau
  damping problem. The initial perturbation undergoes shearing in
  phase space, leading to Landau damping from the phase mixing (see
  previous plot for damping rate). Starting at around :math:`t=20` the
  damping is halted due to particle trapping, finally leading to
  saturation. Phase-space holes are clearly visible.

Tests for the Vlasov-Quasineutral system
----------------------------------------

Etc. etc.
