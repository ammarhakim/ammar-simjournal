-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 1.4 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 4.0, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0, 0.0}, -- lower left corner
   upper = {2.0, 2.0}, -- upper right corner
   cells = {100, 100}, -- number of cells
   cflFrac = 0.9/2.0,
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,
      limiter = "no-limiter",

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]

	 local rho = 1 + 0.2*math.sin(math.pi*(x+y))
	 local u = 1.0
	 local v = -0.5
	 local pr = 1.0
	 local Er = 0.5*rho*(u^2+v^2) + pr/(gasGamma-1)	 
	 
	 return rho, rho*u, rho*v, 0.0, Er
      end,
      
      evolve = true, -- evolve species?
   },   
}
-- run application
eulerApp:run()
