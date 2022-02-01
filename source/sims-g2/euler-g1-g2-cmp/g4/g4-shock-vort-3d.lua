-- Gkyl ------------------------------------------------------------------------
local Hyper = require "App.HyperEqnOnCartGrid"
local Euler = require "Eq.Euler"

-- gas adiabatic index
gasGamma = 1.4
-- wall BCs (these are custom wall BCs for Euler equations)
bcWall = Hyper.bcCustom(Euler.bcWall)

-- create app
eulerApp = Hyper.App {
   logToFile = true, -- false if no log file is desired

   -- basic parameters
   tEnd = 0.8, -- end time
   nFrame = 8, -- number of output frame
   lower = {0.0, 0.0, 0.0}, -- lower left corner
   upper = {1.5, 1.0, 0.5}, -- upper right corner
   cells = {300, 200, 100}, -- number of cells
   cfl = 0.9, -- CFL number
   limiter = "monotonized-centered", -- limiter
   equation = Euler { gasGamma = gasGamma }, -- equation to solve

   -- decomposition stuff
   decompCuts = {1, 1, 1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- initial condition
   init = function (t, xn)
      local rho, pr = 1.0, 1.0
      local rho1, pr1 = 1.0, 10.0
      local rho2, pr2 = 0.1, 1.0
      local x, y, z = xn[1], xn[2], xn[3]
      if (math.sqrt(x^2+y^2) < 0.2) then
	 rho, pr = rho1, pr1
      elseif (math.sqrt((x-0.4)^2+z^2) < 0.2) then
	 rho, pr = rho2, pr2
      end
      return rho, 0.0, 0.0, 0.0, pr/(gasGamma-1)
   end,
   
   -- boundary conditions
   periodicDirs = {}, -- periodic directions
   bcx = { bcWall, Hyper.bcCopy }, -- boundary conditions in X
   bcy = { bcWall, Hyper.bcCopy }, -- boundary conditions in Y
   bcz = { bcWall, Hyper.bcCopy }, -- boundary conditions in Z

   -- diagnostics
   diagnostics = { }
}
-- run application
eulerApp:run()
