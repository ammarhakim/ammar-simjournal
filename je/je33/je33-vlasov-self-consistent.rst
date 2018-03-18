:Author: Ammar Hakim
:Date: March 5th
:Completed: 
:Last Updated:

JE33: Self consistent evolution (continuation of JE 32)
=======================================================

.. contents::

This note is a continuation of the JE32. Here, I account for the
feedback from the currents that are generated as the particles move
about in the specified EM field. The assumption here is that there is
some external process that generates the applied field. As the
particles move in this applied field they generate currents that feed
into Maxwell equations and create fields that linearly combine with
the applied field and hence modify them. The net EM fields are hence

.. math::

  \mathbf{E} &= \mathbf{E}^m + \mathbf{E}^a \\
  \mathbf{B} &= \mathbf{B}^m + \mathbf{B}^a

where :math:`\mathbf{E}^m, \mathbf{B}^m` are fields computed from
Maxwell equations and :math:`\mathbf{E}^a, \mathbf{B}^a` are ones that
are externally appied fields. Note that the external fields, as
generated from currents and charges far away (not included in the
simulation domain) *do not contribute* to the curl terms in Maxwell
equations. Hence, they are not used when updating Maxwell equations.

Uniform time-dependent electric field: resonant case
----------------------------------------------------

In the resonante case (see J
