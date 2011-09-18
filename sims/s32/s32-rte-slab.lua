-- Input file for RTE solution in homogeneous slab

-- define globals for use elsewhere
globals = {tau0 = 1.0}

-- top-level simulation object
simulation = Solver.RteHomogeneousSlab {
   -- number of phase function coefficients (degree of anisotropy)
   L = 8,
   -- number of quadrature points in each hemisphere
   N = 64,
   -- cosine of incident angle
   mu0 = 0.5,
   -- beam flux: downward irradiance is mu0*pi*flux
   flux = 1.0,
   -- optical depth of slab
   tau0 = globals.tau0,
   -- albedo of single scattering
   albedo = 0.95,
   -- number of azimuthal modes
   numModes = 9,
   -- phase function: Mie scattering from spherical particles.
   -- Size parameter \alpha = 2, index of refraction m=1.33
   phaseFunction = RtePhaseFunction.PlCoeffs { 
      coeffs = {
	 1.0, 2.00916, 1.56339, 0.67407, 
	 0.22215, 0.04725, 0.00671, 0.00068, 0.00005
      }
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

-- initialize simulation
simulation:initialize()
-- run simulation (time is irrelevant in this case)
res = simulation:advance(1.0)
-- write data
simulation:write("sol", 1)