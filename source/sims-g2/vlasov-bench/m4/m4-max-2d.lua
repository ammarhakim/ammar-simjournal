-- Gkyl ------------------------------------------------------------------------
--
-- 

local Vlasov = require "App.VlasovOnCartGrid"
local Constants = require "Lib.Constants"

-- Physical parameters
epsilon0 = Constants.EPSILON0
mu0 = Constants.MU0
lightSpeed = Constants.SPEED_OF_LIGHT

L = 1.0
kwave = 2
lwave = 2
freq = 2*math.pi/L*math.sqrt(kwave^2+lwave^2)*lightSpeed
tperiod = 2*math.pi/freq

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = tperiod, -- end time
   nFrame = 2, -- number of output frames
   lower = {0.0, 0.0}, -- configuration space lower left
   upper = {L, L}, -- configuration space upper right
   cells = {16, 16}, -- configuration space cells
   basis = "maximal-order", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {1, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   -- field solver
   field = Vlasov.EmField {
      epsilon0 = epsilon0, mu0 = mu0,
      elcErrorSpeedFactor = 0,
      mgnErrorSpeedFactor = 0,
      
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 local cos = math.cos
	 local pi = math.pi
	 local c = lightSpeed
	 local phi = 2*pi/L*(kwave*x+lwave*y)
	 local knorm = math.sqrt(kwave^2+lwave^2)
	 local kxn, kyn = kwave/knorm, lwave/knorm
	 local E0 = 1.0
	 local Ex, Ey = 0.0, 0.0
	 local Ez = E0*cos(phi)
	 local Bx = 0.0
	 local By = -E0/c*cos(phi)*kxn
	 local Bz = E0/c*cos(phi)*kyn
	 return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- evolve field?
   },   
   
}
-- run application
vlasovApp:run()
