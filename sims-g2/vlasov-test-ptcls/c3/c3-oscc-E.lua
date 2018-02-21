-- Gkyl ------------------------------------------------------------------------
local Vlasov = require "App.VlasovOnCartGrid"

-- Maxwellian in 1x2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 20.0, -- end time
   nFrame = 20, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {2*math.pi}, -- configuration space upper right
   cells = {2}, -- configuration space cells
   cflFrac = 0.5,
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- ions
   ions = Vlasov.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-8.0, -8.0},
      upper = {8.0, 8.0},
      cells = {20, 20},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(1.0, vx, vy, 0.0, 0.0, 1.0)
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M1i" }
   },

   -- field solver
   field = Vlasov.FuncField {
      emFunc = function (t, xn)
	 local omega = 1.0
	 local B0 = 1.0
	 local Ex = 0.5*math.cos(omega*t)
	 return Ex, 0.0, 0.0, 0.0, 0.0, B0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
