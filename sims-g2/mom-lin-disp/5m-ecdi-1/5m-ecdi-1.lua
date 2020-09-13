-- Input file for multimomlinear Tool

local Species = require "Tool.LinearSpecies"

local elcMass = 1
local elcCharge = -1
local ionMass = 400.0
local ionCharge = 1

local E0 = 0.2 -- Electric field in y-direction
local B0  = 1.0 -- Magentic fiels in z-direction
local v0 = E0*B0/B0^2 -- ExB drift speed
local elcTemp = 0.01 -- electron temperature
local ionTemp = 0.01 -- ion temperature


print(string.format("ExB velocity: %g", v0))
print(string.format("Electron cyclotron frequency: %g", math.abs(elcCharge)*B0/elcMass ))
print(string.format("Electron thermal speed %g", math.sqrt(elcTemp/elcMass)))
print(string.format("Ion thermal speed %g", math.sqrt(ionTemp/ionMass)))

-- Electrons
elc = Species.Euler {
   mass = elcMass, -- mass
   charge = elcCharge, -- charge
   density = 1.0, -- number density
   velocity = {v0, 0.0, 0.0}, -- velocity vector
   pressure = elcTemp,
}

-- Ions
ion = Species.Euler {
   mass = ionMass, -- mass
   charge = ionCharge, -- charge
   density = 1.0, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   pressure = ionTemp,

   ignoreBackgroundField = true, -- ions are demagnetized   
}

-- EM field
field = Species.Poisson {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, E0, 0.0}, -- background electric field
   magneticField = {0.0, 0.0, B0}, -- background magnetic field
}

-- list of species to include in dispersion relation
speciesList = { elc, ion }

-- List of wave-vectors for which to compute dispersion relation
kvectors = {}

local kcurr, kmax, NK = 0.0, 40.0, 3001
dk = (kmax-kcurr)/(NK-1)
for i = 1, NK do
   kvectors[i] = {kcurr, 0.0, 0.0} -- each k-vector is 3D
   kcurr = kcurr + dk
end
