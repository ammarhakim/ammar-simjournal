-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local Logger = require "Lib.Logger"

log = Logger { logFile = "info" }

elcMass = 1.0
ionMass = 1836.2

gasGamma = 3.0 -- gas adiabatic constant
LX = 50.0 -- length of domain
mfp = 0.1 -- mean-free path

Te = 1.0e-2 -- temperature
M = 2.0 -- Mach number for flow

-- left state
nl = 1.0 -- left density
pl = nl*Te -- left pressure
ul = M*math.sqrt(gasGamma*pl/(ionMass*nl))

-- right state (creates a stationary shock in ions)
nr = nl*(gasGamma+1)*M^2/((gasGamma-1)*M^2+2)
pr = pl*(2*gasGamma*M^2/(gasGamma+1) - (gasGamma-1)/(gasGamma+1))
ur = ul*(nl/nr)

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
log(string.format("Electron collision frequency: %g\n", nuElc))
log(string.format("Ion collision frequency: %g\n", nuIon))
log(string.format("Ion mfp/LX: %g\n", mfp/LX))
log("\n\n")

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = tEnd, -- end time
   nFrame = 50, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {LX}, -- configuration space upper right
   cells = {1000}, -- configuration space cells
   timeStepper = "fvDimSplit",

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   elc = Plasma.Species {
      charge = -1.0, mass = elcMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,

      -- initial conditions
      init = function (t, xn)
	 local n, u, T = nl, ul, vThermalElc_l^2*elcMass
	 if xn[1]>0.5*LX then
	    n, u, T = nr, ur, vThermalElc_r^2*elcMass
	 end
	 return n*elcMass, n*elcMass*u, 0.0, 0.0, 0.5*n*elcMass*u^2 + n*T/(gasGamma-1)
      end,

      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcCopy, Plasma.Species.bcCopy }, -- boundary conditions in X
   },

   ion = Plasma.Species {
      charge = 1.0, mass = ionMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,

      -- initial conditions
      init = function (t, xn)
	 local n, u, T = nl, ul, vThermalIon_l^2*ionMass
	 if xn[1]>0.5*LX then
	    n, u, T = nr, ur, vThermalIon_r^2*ionMass
	 end
	 return n*ionMass, n*ionMass*u, 0.0, 0.0, 0.5*n*ionMass*u^2 + n*T/(gasGamma-1)
      end,

      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcCopy, Plasma.Species.bcCopy }, -- boundary conditions in X
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      bcx = {Plasma.Field.bcCopy, Plasma.Field.bcCopy},
      evolve = true, -- evolve field?
      bcx = { Plasma.Field.bcCopy, Plasma.Field.bcCopy }, -- boundary conditions in X
   },

   emSource = Plasma.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "time-centered-direct",
   },   
}
-- run application
plasmaApp:run()
