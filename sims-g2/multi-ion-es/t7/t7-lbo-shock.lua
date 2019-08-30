-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell
local Logger = require "Lib.Logger"

log = Logger { logFile = "info" }

-- physical constants
elcMass = 1.0
ionMass = 1836.2
elcCharge = -1.0
ionCharge = 1.0
epsilon0 = 1.0
mu0 = 1.0

gasGamma = 1.66666666667 -- gas adiabatic constant
LX = 50.0 -- length of domain
mfp = 0.1 -- mean-free path

Te = 1.0e-2 -- temperature
pM = 2.0 -- plasma Mach number
M = pM*math.sqrt(2) -- Mach number for flow

-- left state
nl = 1.0 -- left density
pl = nl*Te -- left pressure
ul = M*math.sqrt(gasGamma*pl/(ionMass*nl))

-- right state (creates a stationary shock in ions)
factor = (gasGamma+1.0)*M*M*0.5/((gasGamma-1)*M*M*0.5 + 2.)
nr = nl*factor
pr = (factor*(gasGamma+1.0) - (gasGamma-1.))/((gasGamma+1.0) - factor*(gasGamma-1.0))*pl
ur = ul/factor

-- thermal velocity of ions
vThermalIon_l = math.sqrt(pl/(nl*ionMass))
vThermalIon_r = math.sqrt(pr/(nr*ionMass))

-- thermal velocity electrons
vThermalElc_l = math.sqrt(pl/(nl*elcMass))
vThermalElc_r = math.sqrt(pr/(nr*elcMass))

nuIon = vThermalIon_l/mfp -- ion collision frequency
nuElc = nuIon*math.sqrt(ionMass/elcMass) -- ion collision frequency

tEnd = 0.1*LX/vThermalIon_l -- end time for sim

-- print some information
lambdaD = math.sqrt(epsilon0/elcCharge^2/(nl/Te+nl/Te))

log(string.format("Electron collision frequency: %g\n", nuElc))
log(string.format("Ion collision frequency: %g\n", nuIon))
log(string.format("Ion mfp/LX: %g\n", mfp/LX))
log(string.format("Ion mfp/lambdaD: %g\n", mfp/lambdaD))
log(string.format("ul: %g\n", ul))
log(string.format("ur: %g\n", ur))
log(string.format("vThermalIon_l: %g\n", vThermalIon_l))
log(string.format("vThermalIon_r: %g\n", vThermalIon_r))
log(string.format("vThermalElc_l: %g\n", vThermalElc_l))
log(string.format("vThermalElc_r: %g\n", vThermalElc_r))
log("\n\n")

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = tEnd, -- end time
   nFrame = 50, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {LX}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {16}, -- cuts in each configuration direction
   useShared = true, -- if to use shared memory

   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-0.75, -0.75, -0.75}, 
      upper = {0.75, 0.75, 0.75},
      cells = {12, 12, 12},

      -- initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local n = nl
	    if xn[1]>0.5*LX then
	       n = nr
	    end
	    return n
	 end,
	 driftSpeed = function (t, xn)
	    local ux = ul
	    if xn[1]>0.5*LX then
	       ux = ur
	    end
	    return { ux, 0., 0. }
	 end,
         temperature = function (t, xn)
	    local T = vThermalElc_l^2*elcMass
	    if xn[1]>0.5*LX then
	       T = vThermalElc_r^2*elcMass
	    end
	    return T
	 end,
         exactLagFixM012 = true,
	 isReservoir = true,
      },

      evolveCollisionless = true,
      evolveCollisions = true,
      -- collisions
      lbo = Plasma.LBOCollisions {
	 collideWith = {"elc","ion"},
	 frequencies = {nuElc, nuElc},
      },

      bcx = { Plasma.Species.bcReservoir, Plasma.Species.bcReservoir },

      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },

   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-1.0/math.sqrt(ionMass/elcMass), -1.0/math.sqrt(ionMass/elcMass), -1.0/math.sqrt(ionMass/elcMass)},
      upper = {1.0/math.sqrt(ionMass/elcMass), 1.0/math.sqrt(ionMass/elcMass), 1.0/math.sqrt(ionMass/elcMass)},
      cells = {12, 12, 12},

      -- initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local n = nl
	    if xn[1]>0.5*LX then
	       n = nr
	    end
	    return n
	 end,
	 driftSpeed = function (t, xn)
	    local ux = ul
	    if xn[1]>0.5*LX then
	       ux = ur
	    end
	    return { ux, 0., 0. }
	 end,
         temperature = function (t, xn)
	    local T = vThermalIon_l^2*ionMass
	    if xn[1]>0.5*LX then
	       T = vThermalIon_r^2*ionMass
	    end
	    return T
	 end,
         exactLagFixM012 = true,
	 isReservoir = true,
      },
      
      evolveCollisionless = true,
      evolveCollisions = true,
      -- collisions
      lbo = Plasma.LBOCollisions {
	 collideWith = {"ion", "elc"},
	 frequencies = {nuIon, elcMass/ionMass*nuElc},
      },

      bcx = { Plasma.Species.bcReservoir, Plasma.Species.bcReservoir },

      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },   

   -- field solver
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      bcx = {Plasma.Field.bcCopy, Plasma.Field.bcCopy},
      evolve = true, -- evolve field?
   },   
}
-- run application
plasmaApp:run()
