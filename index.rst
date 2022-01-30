.. Simulation Journal: See daily-notes entry for August 5th 2011

SimJournal: Ammar Hakim's Simulation Journal
============================================

.. note::

  Each note below also has links to the Lua script used to run the
  simulation. Usually, the **links are in figure caption or in
  tables**. They have (unhelpful at first) names like [:doc:`s5
  <../../sims/s5/s5-euler-shock-wave>`]. If you want the exact initial
  conditions, boundary conditions and other simulation details, please
  click those links and look at the Lua script. The initial conditions
  are in (obviously named) functions like init(). The script also
  contains other details like exact setup (resolution, algorithms,
  limiters, time-steps, etc). My goal is that others can reproduce
  these results. Hence, there is no "hidden hand-of-god", as one
  commonly finds in some papers, etc.

.. warning::

   Gkeyll has undergone (and continues to undergo) major changes. The
   old-style Gkeyll 1.0 (G1) input files here will not work with
   Gkeyll 2.0 (G2). **G2 and G1 are completely different codes**. Even
   some G1 input files will need some (minor) mods to work. It is
   nearly impossible for me to keep these input files up-to-date as
   there are, literally, hundreds of them linked to my simulation
   journal.

Below are a set of useful (to me) technical notes and a list of
journal entries, documenting various problems that have I used to
benchmark features in Gkeyll. I can not guarantee that everything here
is correct or accurate. I am very careful in testing and teasing out
physics, but to err is human and some humans err more than others. If
you find any typos or errors please let me know!

For Gkeyll documentation please `see
<https://gkeyll.readthedocs.io/en/latest/>`_.

Technical Notes
---------------

.. toctree::
  :numbered:
  :maxdepth: 1

  moment-eqns
  maxwell-eigensystem
  euler-eigensystem
  tenmom-eigensystem
  twofluid-sources
  fpo-notes
  fluid-kin-drifts
  geometry-metric-symplectic
  hancock-muscl
  multifluid-equilibrium

Programming Notes
-----------------

.. toctree::
   :maxdepth: 1

   pn/pe0/pe0-minimalism.rst

