-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

knumber = 0.5 -- wave-number
perturbation = 1.0e-6 -- distribution function perturbation
fbeam = 1.0 -- beam bag
fmid = 1e-6 -- middle bag
fout = 1e-6 -- outside bags

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = 200.0, -- end time
   nFrame = 400, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   decompCuts = {1},
   useShared = true, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-16.0},
      upper = {16.0},
      cells = {96},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local fv = fout
	 local alpha = perturbation
	 local k = knumber

	 if math.abs(v) < 2 then
	    fv = fmid
	 elseif math.abs(v) < 4 then
	    fv = fbeam
	 else
	    fv = fmid
	 end
	 return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local alpha = perturbation
	 local k = knumber
	 return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
