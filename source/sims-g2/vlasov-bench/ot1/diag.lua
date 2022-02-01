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

print(20/omegaCi)
