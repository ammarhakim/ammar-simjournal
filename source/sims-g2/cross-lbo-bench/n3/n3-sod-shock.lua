-- Gkyl --------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- left/right state for shock
nl, ul, pl = 1.0, 0.0, 1.0
nr, ur, pr = 0.125, 0.0, 0.1

Lx = 1.0 -- domain size
mfp = Lx/100 -- mean-free path

-- thermal velocity to give same energy as in fluid internal energy
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- use left state as reference
nu = vThermal/mfp -- collision frequency

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   local ux, uy = vdrift[1], vdrift[2]
   local vx, vy = v[1], v[2]
   return n0/(2*math.pi*vt^2)*math.exp(-((vx-ux)^2+(vy-uy)^2)/(2*vt^2))
end

sim = Plasma.App {
   logToFile = true,

   tEnd = 0.1, -- end time
   nFrame = 1, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {Lx}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"
   cflFrac = 0.9,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = true, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   neut1 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal, -6.0*vThermal},
      upper = {6.0*vThermal, 6.0*vThermal},
      cells = {8, 8},

      -- initial conditions
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5 then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, {u, 0.0}, vt, {vx, vy})
      end,

      evolveCollisionless = true,
      evolveCollisions = true,
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
   
   neut2 = Plasma.Species {
      charge = 0.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0*vThermal, -6.0*vThermal},
      upper = {6.0*vThermal, 6.0*vThermal},
      cells = {8, 8},

      -- initial conditions
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5 then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, {u, 0.0}, vt, {vx, vy})
      end,

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
}
-- run application
sim:run()
