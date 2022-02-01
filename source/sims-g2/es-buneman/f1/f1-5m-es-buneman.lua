-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

-- electron parameters
vDriftElc = 0.159
vtElc = 0.02
-- ion parameters
vDriftIon = 0.0
vtIon = 0.001
-- mass ratio
massRatio = 25.0

knumber = 1.0 -- wave-number
perturbation = 1.0e-6 -- distribution function perturbation

-- gas gamma
gasGamma = 3.0

momentApp = Moments.App {
   logToFile = true,

   tEnd = 35.0, -- end time
   nFrame = 10, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {128}, -- configuration space cells
   timeStepper = "fvDimSplit",
   suggestedDt = 0.005,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions
   
   elc = Moments.Species {
      charge = -1.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]
	 local rho = 1+perturbation*math.cos(2*math.pi*knumber*x)
	 local rhou = rho*vDriftElc
	 local pr = rho*vtElc^2
	 local Er = 0.5*rhou^2/rho + pr/(gasGamma-1)
	 return rho, rhou, 0.0, 0.0, Er
      end,
      evolve = true, -- evolve species?
   },

   ion = Moments.Species {
      charge = 1.0, mass = massRatio,

      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
	 local rho = massRatio
	 local rhou = rho*vDriftIon
	 local pr = rho*vtIon^2
	 local Er = 0.5*rhou^2/rho + pr/(gasGamma-1)
	 return rho, rhou, 0.0, 0.0, Er
      end,
      evolve = true, -- evolve species?
   },

   -- ghost electrons (carries opposite current to electrons)
   elcGhost = Moments.Species {
      charge = -1.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      init = function (t, xn)
         local x = xn[1]
	 local rho = 1.0
	 local rhou = -rho*vDriftElc
	 local pr = rho*vtElc^2
	 local Er = 0.5*rhou^2/rho + pr/(gasGamma-1)
	 return rho, rhou, 0.0, 0.0, Er
      end,
      evolve = false, -- evolve species?
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x = xn[1]

	 local Ex = -perturbation*math.sin(2*math.pi*knumber*xn[1])/(2*math.pi*knumber)
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion", "elcGhost"},
      timeStepper = "analytic",
   },   

}
-- run application
momentApp:run()
