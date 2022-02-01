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

print(endTime)
