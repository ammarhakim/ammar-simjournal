-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Constants = require "Lib.Constants"

-- Physical parameters
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)

L = 1.0
kwave = 2
freq = 2*math.pi/L*math.sqrt(kwave^2)*lightSpeed
tperiod = L/lightSpeed

-- create app
maxwellApp = Moments.App {
   logToFile = true,

   tEnd = 8*tperiod, -- end time
   nFrame = 16, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {L}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   cflFrac = 0.25,
   timeStepper = "fvDimSplit",

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      limiter = "no-limiter",

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], 0.0
	 local c = lightSpeed
	 local phi = 2*math.pi/L*(kwave*x)
	 local knorm = math.sqrt(kwave^2)
	 local kxn = kwave/knorm
	 local E0 = 1.0/math.sqrt(2.0)
	 local Ex = 0.0
         local Ey = E0*math.cos(phi)
	 local Ez = 0.0
	 local Bx = 0.0
	 local By = 0.0
	 local Bz = E0/c*math.cos(phi)*kxn
	 return Ex, Ey, Ez, Bx, By, Bz
      end,
            
      evolve = true, -- evolve species?
   },   
}
-- run application
maxwellApp:run()
