-- Gkyl ------------------------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- electron parameters
vDriftElc = 0.159
vtElc = 0.02
-- ion parameters
vDriftIon = 0.0
vtIon = 0.001
-- mass ratio
massRatio = 100.0

knumber = 1.0 -- wave-number
perturbation = 1.0e-6 -- distribution function perturbation

local function maxwellian1v(v, vDrift, vt)
   return 1/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vDrift)^2/(2*vt^2))
end

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 100.0, -- end time
   nFrame = 20, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {16}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Vlasov.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vDriftElc},
      upper = {6.0*vDriftElc},
      cells = {64},
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

   -- ghost electrons
   elcGhost = Vlasov.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vDriftElc},
      upper = {6.0*vDriftElc},
      cells = {64},
      decompCuts = {1},

      -- initial conditions
      init = Vlasov.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 1
	 end,
	 driftSpeed = {-vDriftElc},
         temperature = function (t, zn)
	    return vtElc^2
	 end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = false, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },   

   -- electrons
   ion = Vlasov.Species {
      charge = 1.0, mass = massRatio,
      -- velocity space grid
      lower = {-32.0*vtIon},
      upper = {32.0*vtIon},
      cells = {64},
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
      init = function (t, xn)
	 local Ex = -perturbation*math.sin(2*math.pi*knumber*xn[1])/(2*math.pi*knumber)
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
