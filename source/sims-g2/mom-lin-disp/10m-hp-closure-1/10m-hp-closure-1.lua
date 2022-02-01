-- Input file for multimomlinear Tool

local Species = require "Tool.LinearSpecies"

-- Electrons
elc = Species.TenMoment {
   mass = 1.0, -- mass
   charge = -1.0, -- charge
   density = 0.75, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   -- pressure tensor: Pxx, Pxy, Pxz, Pyy, Pyz, Pzz
   pressureTensor = {0.1, 0.0, 0.0, 0.1, 0.0, 0.1},
   useClosure = true,
}

ion = Species.Isothermal {
   mass = 1836.2, -- mass
   charge = 1.0, -- charge
   density = 0.75, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   temperature = 0.1,
}

-- EM field
field = Species.Poisson {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, 0.0, 0.0}, -- background electric field
   magneticField = {0.0, 0.0, 0.0}, -- background magnetic field
}

-- list of species to include in dispersion relation
speciesList = { elc, ion }

-- List of wave-vectors for which to compute dispersion relation
kvectors = { }

local kcurr, kmax, NK = 0.0, 1.0, 401
dk = (kmax-kcurr)/(NK-1)
for i = 1, NK do
   kvectors[i] = {kcurr, 0.0, 0.0} -- each k-vector is 3D
   kcurr = kcurr + dk
end



