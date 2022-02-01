-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

gasGamma = 3.0 -- gas adiabatic constant

-- left state
M = 2.0 -- mach number for flow
rhol = 1.0 -- left density
pl = 1.0 -- left pressure
ul = M*math.sqrt(gasGamma*pl/rhol)

-- right state
rhor = rhol*(gasGamma+1)*M^2/((gasGamma-1)*M^2+2)
pr = pl*(2*gasGamma*M^2/(gasGamma+1) - (gasGamma-1)/(gasGamma+1))
ur = ul*(rhol/rhor)

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of output frame
   lower = {0.0}, -- lower left corner
   upper = {1.0}, -- upper right corner
   cells = {128}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {1}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
	 local xs = 0.5
	 local rhol, ul, pl = rhol, ul, pl
	 local rhor, ur, pr = rhor, ur, pr

	 local rho, u, press = rhor, ur, pr
	 if xn[1]<xs then
	    rho, u, press = rhol, ul, pl
	 end
	 return rho, rho*u, 0.0, 0.0, press/(gasGamma-1) + 0.5*rho*u*u
      end,
      evolve = true, -- evolve species?

      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy }, -- boundary conditions in X
   },   
}
-- run application
eulerApp:run()
