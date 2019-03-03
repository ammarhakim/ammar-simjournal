-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 1.4 -- gas adiabatic constant
   
-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 20.0, -- end time
   nFrame = 200, -- number of output frame
   lower = {-0.5, -0.5}, -- lower left corner
   upper = {0.5, 0.5}, -- upper right corner
   cells = {200, 200}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {1, 2}, -- periodic directions

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 
	 local rho, vx = 1.0, 0.5
	 if math.abs(y)<0.25 then
	    rho, vx = 2.0, -0.5
	 end
	 vx = vx + 0.01*2*(0.5*math.random()-1)
	 local vy = 0.01*2*(0.5*math.random()-1)
	 local pr = 2.5
	 return rho, rho*vx, rho*vy, 0.0, pr/(gasGamma-1) + 0.5*rho*(vx^2+vy^2)
      end,
      
      evolve = true, -- evolve species?
   },   
}
-- run application
eulerApp:run()
