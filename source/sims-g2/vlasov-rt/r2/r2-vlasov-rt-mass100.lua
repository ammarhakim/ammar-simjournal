-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local prng = require "sci.prng"
local Mpi = require "Comm.Mpi"

local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
local rng = prng.mrg32k3a(rank+1234)

-- Constants
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)

elcCharge = -1.0
ionCharge = 1.0
ionMass = 100.0
elcMass = 1.0

-- problem parameters
Atwood = 0.17
n1 = 1.0 -- heavy-side number density
betaIon = 0.071
gHat = 0.09 -- gHat = g/Omega_i v_A
theta = 0.0 -- inclination angle of field to Z-direction
Te_Ti = 0.1 -- Te/Ti
B0 = 0.1 -- Left wall magnetic field (this may not be correct)
maxPert = 0.01 -- Perturbation amplitude

-- secondary quantities
n2 = n1*(1-Atwood)/(1+Atwood) -- light-side number density
T_ion = (betaIon/n1)*B0^2/(2*mu0) -- ion temperature (uniform)
T_elc = Te_Ti*T_ion -- electron temperature (uniform)

vtElc = math.sqrt(2.0*T_elc/elcMass)
vtIon = math.sqrt(2.0*T_ion/ionMass)

wpe = math.sqrt(n1*elcCharge^2/(epsilon0*elcMass))
wpi = math.sqrt(n1*ionCharge^2/(epsilon0*ionMass))
de = lightSpeed/wpe
di = lightSpeed/wpi
Valf = B0/math.sqrt(mu0*n1*ionMass)
OmegaCe0 = elcCharge*B0/elcMass
OmegaCi0 = ionCharge*B0/ionMass

grav = gHat*OmegaCi0*Valf -- gravitational acceleration

-- domain size
Lx = 3.0*di
Ly = 3.75*di

-- resolution and time-stepping
NX = 64
NY = 64
endTime = 60/OmegaCi0

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, temp, mass)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2.0*math.pi*temp/mass)*math.exp(-mass*v2/(2.0*temp))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = endTime, -- end time
   nFrame = 100, -- number of output frames
   lower = {0.0, 0.0}, -- configuration space lower left
   upper = {Lx, Ly}, -- configuration space upper right
   cells = {NX, NY}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {16, 16}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {2}, -- periodic directions

   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,
   restartFrameEvery = 0.01,

   -- electrons 
   elc = Plasma.VlasovSpecies {
      nDistFuncFrame = 10,
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-8.0*vtElc, -8.0*vtElc},
      upper = {8.0*vtElc, 8.0*vtElc},
      cells = {16, 16},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
       local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
       local xloc = 0.5*Lx
       local numDens, Bz
       local B1 = B0
       local maxPert = maxPert -- add this pertubation
       local ph = 0.0
       local xregion = Lx/50.0

       if (x<xloc) then
          -- heavy side
          numDens = n1
       else
          -- light side
          numDens = n2
       end
       local dn = 0.0 --n2*h*5.*math.exp(-(x-xloc)^2/(2.*xregion^2))
       local ne = (numDens+dn)
       local uxe = maxPert*Valf*(rng:sample()-0.5)
       local fv= maxwellian2D(ne, vx, vy, uxe, 0.0, T_elc, elcMass)
       return fv
      end,
      -- gravity is in the x direction and has magnitude of grav
      constGravity = { dir = 1, accel = grav },
      -- boundary conditions
      bcx = { Plasma.VlasovSpecies.bcReflect, Plasma.VlasovSpecies.bcReflect },
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
   },
   -- protons
   ion = Plasma.VlasovSpecies {
      nDistFuncFrame = 10,
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-6.0*vtIon, -6.0*vtIon},
      upper = {6.0*vtIon, 6.0*vtIon},
      cells = {16, 16},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
       local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
       local xloc = 0.5*Lx
       local numDens, Bz
       local B1 = B0
       local maxPert = maxPert -- add this pertubation
       local ph = 0.0
       local xregion = Lx/50.0

       if (x<xloc) then
          -- heavy side
          numDens = n1
       else
          -- light side
          numDens = n2
       end
       local dn = 0.0 --n2*h*5.*math.exp(-(x-xloc)^2/(2.*xregion^2))
       local ni = (numDens+dn)
       local uxi = maxPert*Valf*(rng:sample()-0.5)
       local fv= maxwellian2D(ni, vx, vy, uxi, 0.0, T_ion, ionMass)
       return fv
      end,
      -- gravity is in the x direction and has magnitude of grav
      constGravity = { dir = 1, accel = grav },
      -- boundary conditions
      bcx = { Plasma.VlasovSpecies.bcReflect, Plasma.VlasovSpecies.bcReflect },
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
   },

   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local xloc = 0.5*Lx
         local T = (T_elc + T_ion)
         local B1 = B0
         local maxPert = maxPert -- add this pertubation
   
         if (x<xloc) then
            -- heavy side
            numDens = n1
            Bz = math.sqrt(B1^2 + 2*mu0*(n1*T - numDens*T + ionMass*grav*n1*x))
        else
           -- light side
           numDens = n2
           Bz = math.sqrt(B1^2 + 2*mu0*(n1*T - numDens*T + ionMass*grav*n1*xloc + ionMass*grav*n2*(x-xloc)))
         end
	 return 0.0, 0.0, 0.0, 0.0, 0.0, Bz
      end,
      bcx = { Plasma.MaxwellField.bcReflect, Plasma.MaxwellField.bcReflect },
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
