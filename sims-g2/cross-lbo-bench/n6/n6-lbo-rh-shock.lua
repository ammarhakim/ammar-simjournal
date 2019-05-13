-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

gasGamma = 2.0 -- gas adiabatic constant
LX = 10.0 -- length of domain
mfp = 0.1 -- mean-free path

-- left state
M = 2.0 -- mach number for flow
nl = 1.0 -- left density
pl = 1.0 -- left pressure
ul = M*math.sqrt(gasGamma*pl/nl)

-- right state
nr = nl*(gasGamma+1)*M^2/((gasGamma-1)*M^2+2)
pr = pl*(2*gasGamma*M^2/(gasGamma+1) - (gasGamma-1)/(gasGamma+1))
ur = ul*(nl/nr)

-- thermal velocity to give same energy as in fluid internal energy
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- use left state as reference
nu = vThermal/mfp -- collision frequency

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 0.5*LX/vThermal_l, -- end time
   nFrame = 10, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {LX}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {2}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory


   neut1 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal, -8.0*vThermal},
      upper = {10.0*vThermal, 8.0*vThermal},
      cells = {8, 8},

      -- initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local n = nl
	    if xn[1]>0.5*LX then
	       n = nr
	    end
	    return n
	 end,
	 driftSpeed = function (t, xn)
	    local u = ul
	    if xn[1]>0.5*LX then
	       u = ur
	    end
	    return { u, 0.0 }
	 end,
         temperature = function (t, xn)
	    local T = vThermal_l^2
	    if xn[1]>0.5*LX then
	       T = vThermal_r^2
	    end
	    return T
	 end,
         exactLagFixM012 = true,
      },

      evolveCollisionless = true,
      evolveCollisions = true,
      -- collisions
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut1"},
	 frequencies = {nu},
      },

      bcx = { Plasma.Species.bcOpen, Plasma.Species.bcOpen },

      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i", "u", "vtSq" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },

   neut2 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal, -8.0*vThermal},
      upper = {10.0*vThermal, 8.0*vThermal},
      cells = {8, 8},

      -- initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    local n = nl
	    if xn[1]>0.5*LX then
	       n = nr
	    end
	    return n
	 end,
	 driftSpeed = function (t, xn)
	    local u = ul
	    if xn[1]>0.5*LX then
	       u = ur
	    end
	    return { u, 0.0 }
	 end,
         temperature = function (t, xn)
	    local T = vThermal_l^2
	    if xn[1]>0.5*LX then
	       T = vThermal_r^2
	    end
	    return T
	 end,
         exactLagFixM012 = true,
      },

      -- collisions
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut2"},
	 frequencies = {nu},
      },      

      bcx = { Plasma.Species.bcOpen, Plasma.Species.bcOpen },

      -- diagnostics
      diagnosticMoments = { "M0", "M1i", "M2", "M3i", "u", "vtSq" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },   
}
-- run application
plasmaApp:run()
