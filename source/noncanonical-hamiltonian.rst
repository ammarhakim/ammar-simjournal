Noncanonical Hamiltonian mechanics
==================================

.. contents::

Lagrangian dynamics
-------------------

Let :math:`\mathbf{x}(t)` describe the evolution of a mechanical
system.  Then is well known that the equations of motions determining
the evolution can be derived by finding the minimum of the *action
integral*

.. math::

  S[\mathbf{x}] = \int_{t_0}^{t_1} L(\mathbf{x},\dot{\mathbf{x}},t) dt

where :math:`L(\mathbf{x},\dot{\mathbf{x}},t)` is the Lagrangian of
the system. Note that the action can be computed for any given path
:math:`\mathbf{x}(t)`. However, if :math:`\mathbf{x}(t)` is the path
that minimizes the action integral then perturbing it to
:math:`\mathbf{x}(t)+\delta\mathbf{x}(t)` will not change the action
to first order. This condition, combined with the boundary conditions
:math:`\delta\mathbf{x}(t_0)=\delta\mathbf{x}(t_1)=0`, leads to the
*Euler-Lagrange* equations

.. math::
 :label: eq:euler-lagrange

 \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{x}^i}\right) 
    = \frac{\partial L}{\partial x^i}

which, given the initial conditions :math:`\mathbf{x}(0)` and
:math:`\mathbf{\dot{x}}(0)`, determine the time-evolution of the
system.

For example, the Lagrangian for the motion of a particle in a
time-dependent electromagnetic field is given by

.. math::
 :label: eq:ptcl-lagrangian

  L(\mathbf{x},\dot{\mathbf{x}},t) = \frac{m}{2}|\dot{\mathbf{x}}|^2
   + \frac{e}{c}\dot{\mathbf{x}}\cdot\mathbf{A}(\mathbf{x},t)
   - e \Phi(\mathbf{x},t)

where :math:`m` and :math:`e` are the particle mass and charge
respectively and the scalar and vector potentials, :math:`\Phi` and
:math:`\mathbf{A}`, determine the electromagnetic fields via

.. math::

  \mathbf{E} &= -\nabla\Phi - \frac{1}{c}\frac{\partial
    \mathbf{A}}{\partial t} \\
  \mathbf{B} &= \nabla\times\mathbf{A}.

As the Lagrangian is a scalar it can be transformed to generalized
coordinates :math:`\mathbf{q}(t)` via the transformation
:math:`\mathbf{x}(\mathbf{q},t)` and

.. math::

  \dot{\mathbf{x}} = \frac{\partial \mathbf{x}}{\partial t} 
  +
  \dot{q}^i \frac{\partial \mathbf{x}}{\partial q^i}.

Here and below summation over repeated indices is assumed. This
transformation changes the Lagrangian to
:math:`L(\mathbf{q},\dot{\mathbf{q}},t)`, however, does not change the
form of the Euler-Lagrange equations, Eq. :eq:`eq:euler-lagrange`,
themselves.

Hamiltonian dynamics
--------------------

Instead of using a Lagrangian formulation a *Hamiltonian* formulation
can be used. For this the *canonical momentum* is introduced via

.. math::

  p_i \equiv \frac{\partial L}{\partial \dot{q}^i}.

This can be used to eliminate the :math:`\dot{q}^i` in favor of the
canonical momentum to introduce the *Hamiltonian* function via a
*Legendre transformation* as follows

.. math::
  :label: eq:lt-hamiltonian

  H(\mathbf{q},\mathbf{p},t) = \mathbf{p}\cdot\dot{\mathbf{q}}
    - L(\mathbf{q},\dot{\mathbf{q}}(\mathbf{q},\mathbf{p},t),t)

For example, for the single-particle motion :eq:`eq:ptcl-lagrangian`
the canonical momentum is

.. math::

  p_i = m \dot{x}_i + \frac{e}{c}A_i

and the Hamiltonian becomes

.. math::
 :label: eq:ptcl-hamiltonian

  H(\mathbf{x},\mathbf{p},t) =
  \frac{1}{2m}
  \left|
    \mathbf{p} - \frac{e}{c}\mathbf{A}
  \right|^2
  + e \Phi(\mathbf{x},t).

To derive the equations of motion in terms of the Hamiltonian consider
a variation of :eq:`eq:lt-hamiltonian`

.. math::

   \delta H &= 
     \frac{\partial H}{\partial q^i}\delta q^i
     + \frac{\partial H}{\partial p_i}\delta p_i
     + \frac{\partial H}{\partial t}\delta t \\
   &=
     -\frac{\partial L}{\partial q^i}\delta q^i
     + \dot{q}^i \delta p_i
     - \frac{\partial L}{\partial t}\delta t

From which we get *Hamilton's equations* for the canonical coordinates
:math:`q^i` and :math:`p_i`

.. math::
  :label: eq:hamilton-equations

  \dot{q}^i = \frac{\partial H}{\partial p_i} 
  \quad\mathrm{and}\quad
  \dot{p}_i = -\frac{\partial H}{\partial q^i}.

Phase-space Lagrangian
----------------------

An advantage of the Lagrangian formulation is that one can make
arbitrary coordinate transformations without changing the form of the
Euler-Lagrange equations. Hence, it would be useful to look for a
Lagrangian that yields, as the Euler-Lagrange equations, Hamilton's
equations. This would allow determining Hamilton's equations under an
transform of *both* :math:`\mathbf{q}` and :math:`\mathbf{p}`. The
*phase-space Lagrangian* provides just that. It reads

.. math::
  :label: eq:phase-space-lag
  
  \mathcal{L}(\mathbf{q},\mathbf{p},\dot{\mathbf{q}},\dot{\mathbf{p}},t)
  = \mathbf{p}\cdot\dot{\mathbf{q}} - H(\mathbf{q},\mathbf{p},t).

It is easy to verify that the Euler-Lagrange equations for this
Lagrangian yield Hamilton's equations, Eq. :eq:`eq:hamilton-equations`.

For single-particle motion in an electromagnetic field the phase-space
Lagrangian is

.. math::

  \mathcal{L} = \mathbf{p}\cdot\dot{\mathbf{x}}
  - \frac{1}{2m}
  \left|
    \mathbf{p} - \frac{e}{c}\mathbf{A}
  \right|^2
  - e \Phi(\mathbf{x},t).

As an example of a transformation, consider using the particle
velocity instead of the canonical momentum

.. math::

  \mathbf{v} \equiv \frac{1}{m}\left(
    \mathbf{p} - \frac{e}{c}\mathbf{A}
  \right).

This gives the transformed phase-space Lagrangian

.. math::

  \mathcal{L} = \left(
    m\mathbf{v} + \frac{e}{c}\mathbf{A}(\mathbf{x},t)
  \right)\cdot\dot{\mathbf{x}}
  - \frac{m}{2}|\mathbf{v}|^2 - e\Phi(\mathbf{x},t).

As :math:`\partial \mathcal{L}/\partial \dot{\mathbf{v}} = 0` we get
the kinematic relation :math:`\dot{\mathbf{x}}=\mathbf{v}`. The other
Euler-Lagrange equation leads to

.. math::

  \frac{d}{dt}\left(
    m\mathbf{v} + \frac{e}{c}\mathbf{A}
  \right)
  =
  \nabla\cdot
  \left(
    \frac{e}{c}\mathbf{A} - e\Phi
  \right).
  
Using :math:`d \mathbf{A}/dt =  \partial \mathbf{A}/\partial t +
\mathbf{v}\cdot\nabla\mathbf{A}` and some vector identities leads to
the well known equation of motion

.. math::

  m\dot{\mathbf{v}} = e\mathbf{E} + \frac{e}{c}\mathbf{v}\times\mathbf{B}.

Transformation of the phase-space Lagrangian
--------------------------------------------

The phase-space Lagrangian Eq. :eq:`eq:phase-space-lag` leads to
Hamilton's equations, Eq. :eq:`eq:hamilton-equations`, which have a
very specific form. Consider a general set of coordinates
:math:`z^\alpha`, :math:`\alpha=1,\ldots,2N`, in terms of which we can
write :math:`q^i(\mathbf{z},t)` and :math:`p_i(\mathbf{z},t)`. Note
that for a N degree of freedom system we need 2N general coordinates
:math:`z^\alpha`. The question we now ask is: what are the equations
of motion in these new coordinates?  Note that the transformation to
the coordinates :math:`\mathbf{z}` is completely arbitrary: we can
not, in general, pick out half of the coordinates as "positions" and
other half as "generalized momentum", and the phase-space may be
completely mixed in the new coordinate system.

In these new coordinates we have

.. math::

  \dot{q}^i = \frac{\partial q^i}{\partial t} 
    + \dot{z}^\alpha \frac{\partial q^i}{\partial z^\alpha}.

Using this in the phase-space Lagrangian we get

.. math::
  :label: eq:phase-space-lag-z

  \mathcal{L}(\mathbf{z},t)
  = \Lambda_\alpha \dot{z}^\alpha - \mathcal{H}

where

.. math::

  \Lambda_\alpha \equiv p_i\frac{\partial q^i}{\partial z^\alpha}
  
and

.. math::

  \mathcal{H} \equiv H - p_i \frac{\partial q^i}{\partial t}.

The first term Eq. :eq:`eq:phase-space-lag-z` is called the
*symplectic* part and the second term is called the Hamiltonian part.

The Euler-Lagrange equations corresponding to this transformed
Lagrangian is

.. math::

 \frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{z}^\alpha}\right) 
    = \frac{\partial \mathcal{L}}{\partial z^\alpha}.

We have

.. math::

  \frac{\partial \mathcal{L}}{\partial \dot{z}^\alpha}
  &=
  \Lambda_\alpha \\
  \frac{\partial \mathcal{L}}{\partial z^\alpha}
  &=
  \frac{\partial \Lambda_\beta}{\partial z^\alpha}\dot{z}^\beta
  - 
  \frac{\partial \mathcal{H}}{\partial z^\alpha}.

From this the Euler-Lagrange equations give

.. math::

  \frac{d\Lambda_\alpha}{dt}
  \equiv 
  \frac{\partial \Lambda_\alpha}{\partial t} + \dot{z}^\beta
  \frac{\partial \Lambda_\alpha}{\partial z^\beta}
  =
  \frac{\partial \Lambda_\beta}{\partial z^\alpha}\dot{z}^\beta
  - 
  \frac{\partial \mathcal{H}}{\partial z^\alpha}.

This gives the equations of motion

.. math::

  \omega_{\alpha\beta}\thinspace\dot{z}^\beta
  =
  \frac{\partial \mathcal{H}}{\partial z^\alpha}
  +
  \frac{\partial \Lambda_\alpha}{\partial t}.

Here the *Lagrange matrix* is defined as

.. math::

  \omega_{\alpha\beta} \equiv
  \frac{\partial \Lambda_\beta}{\partial z^\alpha}
  -
  \frac{\partial \Lambda_\alpha}{\partial z^\beta}.

Assuming that :math:`\det(\boldsymbol{\omega})` is non-singular the
explicit form of the equations of motion can be written as

.. math::
  :label: eq:zevolve

  \dot{z}^\beta
  =
  \Pi^{\beta\alpha}\left(\frac{\partial \mathcal{H}}{\partial z^\alpha}
  +
  \frac{\partial \Lambda_\alpha}{\partial t}
  \right)

where the *Poisson structure* :math:`\mathbf{\Pi} =
\boldsymbol{\omega}^{-1}`. This are the equations of motions we have
been seeking.

Canonical transforms, Symplectic matrices
-----------------------------------------

Let the vector :math:`\mathbf{M} = (q^1,\ldots,q^N, p_1,\ldots,p_N)`
represent canonical coordinates. Then the equations of motion
Eq. :eq:`eq:hamilton-equations` can be written in the compact form

.. math::

 \dot{M}^\alpha = \sigma^{\alpha \beta}\frac{\partial H}{\partial M^\beta}

where the *fundamental symplectic matrix* :math:`\boldsymbol{\sigma}`
is defined as

.. math::

  \boldsymbol{\sigma}
  =
  \left(
   \begin{matrix}
    \mathbf{0} & \mathbf{I} \\
    -\mathbf{I} & \mathbf{0}
   \end{matrix}
  \right)

where :math:`\mathbf{I}` is a :math:`2N\times 2N` unit matrix. For a
time-independent transformation :math:`M^\alpha(\mathbf{z})`
Eq. :eq:`eq:zevolve` shows that

.. math::
 
  \dot{z}^\alpha
  =
  \Pi^{\alpha\beta}\frac{\partial H}{\partial z^\beta}.

Substituting :math:`M^\alpha(\mathbf{z})` in the canonical equations
of motion and comparing with the above equation it is clear that if

.. math::
  :label: eq:symplectic-def

  \mathbf{D}\thinspace\mathbf{\sigma}\thinspace\mathbf{D}^T 
  = \mathbf{\sigma},

where :math:`D^\alpha_\beta = \partial z^\alpha/\partial M^\beta` is
the Jacobian matrix of the transformation, then the *form* of
Hamilton's equations remains unchanged. The class of all such
transformation that preserve the form of the Hamilton's equations is
called *canonical transforms*. All matrices :math:`\mathbf{D}` that
satisfy :eq:`eq:symplectic-def` are called *symplectic matrices*. For
arbitrary time-independent transforms, however, the relation
Eq. :eq:`eq:symplectic-def` will not hold, i.e. the Jacobian matrix of
the transformation will not be symplectic.

Poisson brackets
----------------

With the matrices :math:`\Pi^{\alpha\beta}` we can define the *Poisson
bracket* of two functions :math:`f(\mathbf{z},t)` and
:math:`g(\mathbf{z},t)` as

.. math::

  \{f,g\} \equiv 
  \frac{\partial f}{\partial z^\alpha}
  \Pi^{\alpha\beta}
  \frac{\partial g}{\partial z^\beta}.

For canonical transforms this reduces to

.. math::

  \{f,g\} &=
  \frac{\partial f}{\partial z^\alpha}
  \sigma^{\alpha\beta}
  \frac{\partial g}{\partial z^\beta} \\
  &=
  \frac{\partial f}{\partial q^i}\frac{\partial g}{\partial p_i}
  -
  \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q^i}.

Starting from the equations of motion in canonical coordinates
:math:`\mathbf{M}` (defined in the previous section) for a
time-independent transformation :math:`M^\alpha(\mathbf{z})` we can
show that

.. math::
 
  \dot{z}^\alpha
  =
  \Pi^{\alpha\beta}\frac{\partial H}{\partial z^\beta}
  = 
  \frac{\partial z^\alpha}{\partial M^\delta}
  \sigma^{\delta\gamma}
  \frac{\partial z^\beta}{\partial M^\gamma}
  \frac{\partial H}{\partial z^\beta}
  = 
  \{z^\alpha,z^\beta\} \frac{\partial H}{\partial z^\beta}.

which indicates that

.. math::

  \Pi^{\alpha\beta} = \{z^\alpha,z^\beta\}

where the Poisson bracket is defined with respect to the canonical
coordinates. Note that once the expression for the Poisson brackets
are known the equation of motion can be written in the compact form

.. math::

  \dot{z}^\alpha = \{z^\alpha,H\}.

Liouville theorem
-----------------

Let :math:`\mathcal{J}=\det(\mathbf{D}^{-1})` be the Jacobian of the
transformation, where :math:`D^\alpha_\beta = \partial
z^\alpha/\partial M^\beta`. For a time-dependent transformation the
Jacobian satisfies

.. math::

 \frac{\partial \mathcal{J}}{\partial t}
 +
 \frac{\partial} {\partial z^\alpha}
   \left(\dot{z}^\alpha\mathcal{J}\right)
 =
 0.

This indicates that the equations of motion satisfy the Liouville
theorem, that is, the Hamiltonian flow conserves phase-space volume
:math:`d\mathbf{M} = \mathcal{J}d\mathbf{z}`.

For time-independent transforms this gives

.. math::

 0
 =
  \frac{\partial} {\partial z^\alpha}\left(
   \mathcal{J}\Pi^{\alpha\beta}
   \frac{\partial H}{\partial z^\beta}
  \right)
 =
  \frac{\partial} {\partial z^\alpha}\left(
   \mathcal{J}\Pi^{\alpha\beta}
  \right)
   \frac{\partial H}{\partial z^\beta}

where the antisymmetry of the :math:`\Pi^{\alpha\beta}` was used, and
which yields the *Liouville identities*

.. math::

 \frac{\partial} {\partial z^\alpha}\left(
   \mathcal{J}\Pi^{\alpha\beta}
  \right)
  =
  0.

As can be verified, this allows writing the noncanonical Poisson
bracket as a phase-space divergence

.. math::

 \{f,g\}
 =
 \frac{1}{\mathcal{J}}
 \frac{\partial} {\partial z^\alpha}\left(
   f \mathcal{J}\Pi^{\alpha\beta}
   \frac{\partial g}{\partial z^\beta}
  \right).

