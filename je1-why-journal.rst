JE1: What is a simulation journal?
====================================================

.. epigraph::

  An article about computational science in a scientific publication
  is not the scholarship itself, it is merely advertising of the
  scholarship. The actual scholarship is the complete software
  development environment and the complete set of instructions which
  generated the figures.

  -- J.B. Buckheit, and D. L. Donoho

*Reproducible* research is an essential component of scientific
discovery. Unfortunately, a lot computational physics research is
*not* reproducible. Very often results and figures in published papers
are hard or even impossible to recreate. This makes it difficult to
compare algorithms, understand and recreate published results and
slows down the process of scientific discovery.

Maintaining a simulation journal, analogous to a laboratory notebook,
may promote reproducibility. A simulation journal can contain more
details than a published paper, provide annotated input files and code
needed to produce tables, figures and results of the calculation.

What this means in practice is that for all computations we need to
store the following information.

- For *each* simulation, the complete input needed to run the
  calculation.

- For *each* simulation, the scripts/programs and data required to
  create the results (tables, figures, etc).

- A journal entry describing initial conditions, boundary conditions
  and any special notes on running the simulations. There need not be
  one entry per simulation, but each simulation must be described in
  *some* entry.

- The results presented and described in as much detail as needed.

The journal entry should written in a manner that it can be eventually
incorporated into a publication and should be explicitly written with
reproducibility in mind. The entry need not be long and drawn out. It
should be precise and just detailed enough for a dedicated reader to
reproduce the results presented.
