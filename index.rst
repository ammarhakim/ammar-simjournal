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

Below are a list of journal entries, documenting various problems that
have been attempted with Gkeyll. The eventual goal of Gkeyll is to
solve the gyrokinetic equations in the edge region of tokamaks,
including the scrap-off-layer. However, Gkeyll provides a powerful
framework to study various physics problems as well as test different
algorithms in a modular way. Some of these are documented below, with
links to the Lua scripts to run those problems.

Technical Notes
---------------

.. toctree::
  :maxdepth: 1

  maxwell-eigensystem
  euler-eigensystem
  tenmom-eigensystem
  twofluid-sources
  hancock-muscl
  noncanonical-hamiltonian
  multifluid-equilibrium

Gkeyll Simulation Journal
-------------------------

.. toctree::
  :maxdepth: 1

  sims/simindex.rst
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

