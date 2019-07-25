-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local Mpi = require "Comm.Mpi"

gasGamma = 1.4 -- gas adiabatic constant
mach = 0.1 -- Mach number
pr0 = 1.0 -- Reference pressure
rho0 = 1.0 -- Reference density
cs0 = math.sqrt(gasGamma*pr0/rho0)
V0 = cs0*mach -- reference speed
pert = 0.0 -- perturbation

-- set random seed based on processor number
math.randomseed(100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())

-- create app
eulerApp = Moments.App {
   logToFile = true,

   tEnd = 20*(1/V0), -- end time
   nFrame = 200, -- number of output frame
   lower = {0.0, 0.0, 0.0}, -- lower left corner
   upper = {2.0*math.pi, 2.0*math.pi, 2.0*math.pi}, -- upper right corner
   cells = {96, 96, 96}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   timeStepper = "fvDimSplit",
   
   -- decomposition stuff
   decompCuts = {6, 6, 6}, -- cuts in each direction
   useShared = false, -- if to use shared memory

   periodicDirs = {1, 2, 3}, -- periodic directions

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         local sin = math.sin
         local cos = math.cos
         local vx = V0*sin(x)*cos(y)*cos(z) + V0*pert*2*(0.5*math.random()-1)
         local vy = -V0*cos(x)*sin(y)*cos(z) + V0*pert*2*(0.5*math.random()-1)
         local vz = V0*pert*2*(0.5*math.random()-1)
	 local pr = pr0 + rho0*V0^2/16*(cos(2*x)+cos(2*y))*(cos(2*z)+2)
	 return rho0, rho0*vx, rho0*vy, rho0*vz, pr/(gasGamma-1) + 0.5*rho0*(vx^2+vy^2+vz^2)
      end,
      
      evolve = true, -- evolve species?
   },   
}
-- run application
eulerApp:run()

