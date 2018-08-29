-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"
local Mpi = require "Comm.Mpi"

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0 = 40*eV 
Ti0 = 40*eV 
B_axis = 0.5 -- [T]
R0 = 0.85  -- [m]
a0  = 0.5 -- [m]
R = R0 + a0
B0 = B_axis*(R0/R) -- [T]
n0 = 7e18 -- [1/m^3]
P_SOL = 8.1e5 -- [W] 
S0 = 5.7691e23
xSource = R - 0.05 -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature
-- derived parameters
vti     = math.sqrt(Ti0/mi)
vte  	= math.sqrt(Te0/me)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s   = c_s/omega_ci

-- box size
Lx = 50*rho_s
Ly = 100*rho_s
Lz = 4 -- [m]

-- source profiles
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local sourceFloor = 0.1
   if math.abs(z) < Lz/4 then
      return 0.90625*S0*math.max(math.exp(-(x-xSource)^2/(2*lambdaSource)^2), sourceFloor)
   else
      return 1e-10
   end
end
sourceTemperature = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < xSource + 3*lambdaSource then
      return 80*eV
   else
      return 30*eV
   end
end

-- initialize a random seed for initial conditions
-- will be used for both ions and electrons
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 1e-4, -- end time
   nFrame = 100, -- number of output frames
   lower = {R - Lx/2, -Ly/2, -Lz/2}, -- configuration space lower left
   upper = {R + Lx/2, Ly/2, Lz/2}, -- configuration space upper right
   cells = {12, 24, 8}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 0.2,
   restartFrameEvery = 0.01,

   -- decomposition for configuration space
   decompCuts = {4, 8, 2}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {2}, -- periodic in y only

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = qe,
      mass = me,
      lower = {-4*vte, 0},
      upper = {4*vte, 12*me*vte^2/(2*B0)},
      cells = {8, 4},
      decompCuts = {1, 1},
      -- initial conditions
      init = {"maxwellian", 
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local Ls = Lz/4
                 local effectiveSource = sourceDensity(t,{x,y,0})
                 local c_ss = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb = 1e-3*math.random(-1,1)
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if (x < xSource + 3*lambdaSource) then 
                    return 50*eV
                 else 
                    return 20*eV
                 end
              end,
      },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperature},
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkM1"}, 
      randomseed = randomseed,
      bcx = {Plasma.GkSpecies.bcZeroFlux, Plasma.GkSpecies.bcZeroFlux},
      bcz = {Plasma.GkSpecies.bcSheath, Plasma.GkSpecies.bcSheath},
   },

   -- gyrokinetic ions
   ion = Plasma.GkSpecies {
      charge = qi,
      mass = mi,
      -- velocity space grid
      lower = {-4*vti, 0},
      upper = {4*vti, 12*mi*vti^2/(2*B0)},
      cells = {16, 8},
      decompCuts = {1, 1},
      -- initial conditions
      init = {"maxwellian", 
              density = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local Ls = Lz/4
                 local effectiveSource = sourceDensity(t,{x,y,0})
                 local c_ss = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb = 1e-3*math.random(-1,1)
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if x < xSource + 3*lambdaSource then 
                    return 50*eV
                 else 
                    return 20*eV
                 end
              end,
              driftSpeed = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local Te
                 if x < xSource + 3*lambdaSource then 
                    Te = 50*eV
                 else 
                    Te = 20*eV
                 end
                 if math.abs(z) <= Lz/4 then
		    return z/(Lz/4)*math.sqrt(Te/mi)
                 else
		    return z/math.abs(z)*math.sqrt(Te/mi)
                 end
              end,
      },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperature},
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkM0", "GkM1"}, 
      randomseed = randomseed,
      bcx = {Plasma.GkSpecies.bcZeroFlux, Plasma.GkSpecies.bcZeroFlux},
      bcz = {Plasma.GkSpecies.bcSheath, Plasma.GkSpecies.bcSheath},
   },

   -- field solver
   field = Plasma.GkField {
      -- dirichlet in x
      phiBcLeft = { T ="D", V = 0.0},
      phiBcRight = { T ="D", V = 0.0},
      -- periodic in y --
      -- no bc in z
      phiBcBack = { T ="N", V = 0.0},
      phiBcFront = { T ="N", V = 0.0},
      evolve = true, -- evolve fields?
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R/x
      end,

      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
