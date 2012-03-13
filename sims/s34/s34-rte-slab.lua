-- Input file for RTE solution in homogeneous slab

-- define globals for use in simulation
globals = {
   tau0 = 1.0,
   albedo = 0.9,
   mu0 = 0.5,
}

-- Read phase-function coefficients from a file
-- @param fname Name of file to read.
-- @return Table of read coefficients.
function readPfCoeffsFromFile(fname)
   io.input(fname)
   local coeffs = {}
   while true do
      local n1 = io.read("*number")
      if not n1 then break end
      coeffs[#coeffs+1] = n1
   end
   return coeffs
end

-- top-level simulation object
simulation = Solver.RteHomogeneousSlab {
   -- number of phase function coefficients (degree of anisotropy)
   L = 82,
   -- number of quadrature points in each hemisphere
   N = 128,
   -- cosine of incident angle
   mu0 = globals.mu0,
   -- beam flux: downward irradiance is mu0*pi*flux
   flux = 1.0,
   -- optical depth of slab
   tau0 = globals.tau0,
   -- albedo of single scattering
   albedo = globals.albedo,
   -- number of azimuthal modes
   numModes = 64,
   -- phase function: Haze-L read from file
   phaseFunction = RtePhaseFunction.PlCoeffs { 
      coeffs = readPfCoeffsFromFile("hazel")
   },
   -- set of dummy nodes at which to compute radiances
   dummyNodes = {
      0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
   },
   -- optical depths at which outputs are to be computed
   tauRadOut = {
      0.0, globals.tau0/20, globals.tau0/10, globals.tau0/5,
      globals.tau0/2, 3*globals.tau0/4, globals.tau0
   },
   -- irradiance moments to compute
   irradOut = {0, 1},
   -- optical depths at which irradiance are to be computed
   tauIrradOut = {
      0.0, globals.tau0/20, globals.tau0/10, globals.tau0/5,
      globals.tau0/2, 3*globals.tau0/4, globals.tau0
   }
}

-- run simulation (time is irrelevant in this case)
res = simulation:advance(1.0)
-- write data
simulation:write("sol", 1)