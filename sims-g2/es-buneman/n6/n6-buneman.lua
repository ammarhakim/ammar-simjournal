-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- electron parameters
vDriftElc = 0.159
vtElc = 0.02
-- ion parameters
vDriftIon = 0.0
vtIon = 0.001
-- mass ratio
massRatio = 1836.2

knumber = 1.0 -- wave-number
perturbation = 1.0e-6 -- distribution function perturbation

local function maxwellian1v(v, vDrift, vt)
   return 1/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vDrift)^2/(2*vt^2))
end

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 2000.0, -- end time
   nFrame = 200, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {128}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 0.9,

   -- decomposition for configuration space
   decompCuts = {4}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Vlasov.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vDriftElc},
      upper = {6.0*vDriftElc},
      cells = {256},
      decompCuts = {1},

      -- initial conditions
      init = Vlasov.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 1+perturbation*math.cos(2*math.pi*knumber*x)
	 end,
	 driftSpeed = {vDriftElc},
         temperature = function (t, zn)
	    return vtElc^2
	 end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- electrons
   ion = Vlasov.Species {
      charge = 1.0, mass = massRatio,
      -- velocity space grid
      lower = {-128.0*vtIon},
      upper = {128.0*vtIon},
      cells = {256},
      decompCuts = {1},

      -- initial conditions
      init = Vlasov.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 1
	 end,
	 driftSpeed = {vDriftIon},
         temperature = function (t, zn)
	    return vtIon^2*massRatio
	 end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },   

   -- field solver
   field = Vlasov.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      useGhostCurrent = true,
      init = function (t, xn)
	 local Ex = -perturbation*math.sin(2*math.pi*knumber*xn[1])/(2*math.pi*knumber)
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
