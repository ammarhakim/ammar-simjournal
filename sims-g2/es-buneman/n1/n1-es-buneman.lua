-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

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

local function maxwellian1v(v, vDrift, vt)
   return 1/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vDrift)^2/(2*vt^2))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 150.0, -- end time
   nFrame = 150, -- number of output frames
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
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vDriftElc},
      upper = {6.0*vDriftElc},
      cells = {256},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local fv = maxwellian1v(v, vDriftElc, vtElc)
	 return fv*(1+perturbation*math.cos(2*math.pi*knumber*x))
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- ghost electrons: balances the electrons to cancel initial currents
   elcGhost = Plasma.FuncSpecies {
      charge = -1.0, mass = 1.0,
      -- momentum provided by this species
      momentumDensity = function (t, xc)
   	 local n0 = 1.0
   	 return -n0*vDriftElc
      end,
      
      evolve = false, -- evolve species?
   },   

   -- electrons
   ion = Plasma.Species {
      charge = 1.0, mass = massRatio,
      -- velocity space grid
      lower = {-128.0*vtIon},
      upper = {128.0*vtIon},
      cells = {256},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian1v(v, vDriftIon, vtIon)
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },   

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local Ex = -perturbation*math.sin(2*math.pi*knumber*xn[1])/(2*math.pi*knumber)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
