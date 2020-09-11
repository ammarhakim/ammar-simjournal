-- Input file for multimomlinear Tool

local Species = require "Tool.LinearSpecies"

local elcMass = 1
local elcCharge = -1
local ionMass = 1836.2
local ionCharge = 1
local elcTemp = 1e-6
local ud = 0.1

-- Electrons
elc1 = Species.TenMoment {
   mass = elcMass, -- mass
   charge = elcCharge, -- charge
   density = 1.0, -- number density
   velocity = {0.0, ud, 0.0}, -- velocity vector
   pressureTensor = {elcTemp, 0.0, 0.0, elcTemp, 0.0, elcTemp},
}

elc2 = Species.TenMoment {
   mass = elcMass, -- mass
   charge = elcCharge, -- charge
   density = 1.0, -- number density
   velocity = {0.0, -ud, 0.0}, -- velocity vector
   pressureTensor = {elcTemp, 0.0, 0.0, elcTemp, 0.0, elcTemp},
}

-- Ions
ion = Species.Isothermal {
   mass = ionMass, -- mass
   charge = ionCharge, -- charge
   density = 2.0, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   temperature = 0.0, -- temperature
}

-- EM field
field = Species.Maxwell {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, 0.0, 0.0}, -- background electric field
   magneticField = {0.0, 0.0, 0.0}, -- background magnetic field
}

-- list of species to include in dispersion relation
speciesList = { elc1, elc2, ion }

-- List of wave-vectors for which to compute dispersion relation
kvectors = {}

local kcurr, kmax, NK = 0.0, 10.0, 401
dk = (kmax-kcurr)/(NK-1)
for i = 1, NK do
   kvectors[i] = {kcurr, 0.0, 0.0} -- each k-vector is 3D
   kcurr = kcurr + dk
end
