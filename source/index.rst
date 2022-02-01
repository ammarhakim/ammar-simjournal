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

   pn/pn0/pn0-minimalism.rst  

Simulation Journal
------------------

.. toctree::
  :maxdepth: 1

  je/je0/je0-repro-research.rst
  je/je1/je1-periodic-poisson.rst
  je/je2/je2-euler-shock.rst
  je/je3/je3-homoslab-rte.rst
  je/je4/je4-twofluid-shock.rst
  je/je5/je5-dispersive-eqns.rst
  je/je6/je6-maxwell-solvers.rst
  je/je7/je7-dual-yee.rst
  je/je8/je8-plasmabeach.rst
  je/je9/je9-cyclotron-tunneling.rst
  je/je10/je10-icw.rst
  je/je11/je11-fem-poisson.rst
  je/je12/je12-poisson-bracket.rst
  je/je13/je13-incomp-euler-2d.rst
  je/je14/je14-vlasov-fixed-pot.rst
  je/je15/je15-vlasov-poisson.rst
  je/je16/je16-ldg.rst
  je/je17/je17-hasegawa-wakatani.rst
  je/je18/je18-GEM-recon.rst
  je/je19/je19-rdg.rst
  je/je20/je20-vlasov-bounded.rst
  je/je21/je21-lin-gke-em.rst
  je/je22/je22-euler-2d.rst
  je/je23/je23-euler-3d.rst
  je/je24/je24-euler-embedded-bc.rst
  je/je25/je25-three-wave-bra-bba.rst
  je/je26/je26-dg-maxwell.rst
  je/je28/je28-boltz-bgk.rst
  je/je29/je29-es-shock.rst
  je/je30/je30-coupled-hw.rst
  je/je31/je31-antilim-adv.rst
  je/je32/je32-vlasov-test-ptcl.rst
  je/je33/je33-buneman.rst
  je/je34/je34-linear-dispersion.rst
  je/je35/je35-iter-poisson.rst   
