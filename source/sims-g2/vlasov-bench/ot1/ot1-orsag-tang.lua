-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

Te_Ti = 1.0 -- ratio of electron to ion temperature
n0 = 1.0 -- initial number density

elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge
ionMass = 25.0 -- ion mass
ionCharge = 1.0 -- ion charge

vAe = 0.1
B0 = vAe*math.sqrt(mu0*n0*elcMass)
beta = 0.16
vtElc = vAe*math.sqrt(beta)

elcTemp = vtElc^2/2.0
ionTemp = elcTemp/Te_Ti

-- ion velocities
vAi = vAe/math.sqrt(ionMass)
vtIon = vtElc/math.sqrt(ionMass) --Ti/Te = 1.0

-- cyclotron frequencies
omegaCe = ionCharge*B0/elcMass
omegaCi = ionCharge*B0/ionMass

-- inertial length
de = vAe/omegaCe
di = vAi/omegaCi

-- OT initial conditions
u0x = 0.2*vAi
u0y = 0.2*vAi
B0x = 0.2*B0
B0y = 0.2*B0
ni0 = n0 -- Guess for initial value of ni0

-- domain size and simulation time
Lx = 20.48*di
Ly = 20.48*di

-- Maxwellian in 1x2v
local function maxwellian3D(n, vx, vy, vz, ux, uy, uz, temp, mass)
   local v2 = (vx - ux)^2 + (vy - uy)^2 + (vz - uz)^2
   return n/math.sqrt((2*math.pi*temp/mass)^3)*math.exp(-mass*v2/(2.0*temp))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 50.0/omegaCi, -- end time
   nFrame = 500, -- number of output frames
   lower = {0.0, 0.0}, -- configuration space lower left
   upper = {Lx, Ly}, -- configuration space upper right
   cells = {64, 64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {16, 16}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory
   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions
   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,
   restartFrameEvery = 0.01,

   -- electrons
   elc = Plasma.VlasovSpecies {
      nDistFuncFrame = 100,
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-6.0*vtElc, -6.0*vtElc, -6.0*vtElc},
      upper = {6.0*vtElc, 6.0*vtElc, 6.0*vtElc},
      cells = {12, 12, 12},
      decompCuts = {1, 1, 1},
      -- initial conditions
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0
         local ni = n0 - (epsilon0/qi)*(_2pi/(Lx^2.*Ly^2.*n0*qi*mu0))*(B0*Lx*Ly^2.*u0y*ni0*qi*mu0*cos(_2pi*x/Lx) + 8*Pi*B0y^2.*Ly^2.*cos(8*Pi*x/Lx) + Lx*(Ly*(B0*Lx*u0x*ni0*qi*mu0 + 8*Pi*B0x*B0y*cos(_4pi*x/Lx))*cos(_2pi*y/Ly)+_2pi*B0x^2*Lx*cos(_4pi*y/Ly)))
   
         local vdrift_x = -u0x*sin(_2pi*y/Ly)*ni
         local vdrift_y = u0y*sin(_2pi*x/Lx)*ni
         local vdrift_z = -Jz / qi

	 local fv = maxwellian3D(n0, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, elcTemp, elcMass)

	 return fv
      end,
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },
   -- protons
   ion = Plasma.VlasovSpecies {
      nDistFuncFrame = 100,
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-6.0*vtIon, -6.0*vtIon, -6.0*vtIon},
      upper = {6.0*vtIon, 6.0*vtIon, 6.0*vtIon},
      cells = {12, 12, 12},
      decompCuts = {1, 1, 1},
      -- initial conditions
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local ni = n0 - (epsilon0/qi)*(_2pi/(Lx^2.*Ly^2.*n0*qi*mu0))*(B0*Lx*Ly^2.*u0y*ni0*qi*mu0*cos(_2pi*x/Lx) + 8*Pi*B0y^2.*Ly^2.*cos(8*Pi*x/Lx) + Lx*(Ly*(B0*Lx*u0x*ni0*qi*mu0 + 8*Pi*B0x*B0y*cos(_4pi*x/Lx))*cos(_2pi*y/Ly)+_2pi*B0x^2*Lx*cos(_4pi*y/Ly)))
  
         local vdrift_x = -u0x*sin(_2pi*y/Ly)*ni
         local vdrift_y = u0y*sin(_2pi*x/Lx)*ni
         local vdrift_z = 0.0

	 local fv = maxwellian3D(ni, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, ionTemp, ionMass)

	 return fv
      end,
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },

   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local ne = n0
         local ni = n0 - (epsilon0/qi)*(_2pi/(Lx^2.*Ly^2.*n0*qi*mu0))*(B0*Lx*Ly^2.*u0y*ni0*qi*mu0*cos(_2pi*x/Lx) + 8*Pi*B0y^2.*Ly^2.*cos(8*Pi*x/Lx) + Lx*(Ly*(B0*Lx*u0x*ni0*qi*mu0 + 8*Pi*B0x*B0y*cos(_4pi*x/Lx))*cos(_2pi*y/Ly)+_2pi*B0x^2*Lx*cos(_4pi*y/Ly)))
         local Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0

         local Bx = -B0x*sin(_2pi*y/Ly)
         local By = B0y*sin(_4pi*x/Lx)
         local Bz = B0

         -- Assumes qi = abs(qe)
         local u_xe = -u0x*sin(_2pi*y/Ly)*ni/ne
         local u_ye = u0y*sin(_2pi*x/Lx)*ni/ne
         local u_ze = -Jz / (qi*ne)
   
         -- E = - v_e x B ~  (J - u) x B
         local Ex = - (u_ye*Bz - u_ze*By)
         local Ey = - (u_ze*Bx - u_xe*Bz)
         local Ez = - (u_xe*By - u_ye*Bx)

	 return Ex, Ey, Ez, Bx, By, Bz
      end,

      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
