-- Gkyl ------------------------------------------------------------------------
local Vlasov = require "App.VlasovOnCartGrid"

epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

Te_Ti = 1.0 -- ratio of electron to ion temperature
nL = 1.0 -- initial number density
nR = 0.125 -- initial number density on high side

elcTempL = 1.0e-2 -- electron temperature on left
elcTempR = 0.8*elcTempL -- electron temperature on right
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

ionTempL = elcTempL/Te_Ti -- ion temperature on left
ionTempR = elcTempR/Te_Ti -- ion temperature on right
ionMass = 25.0 -- ion mass
ionCharge = 1.0 -- ion charge

-- thermal speeds on left
vtElcL = math.sqrt(elcTempL/elcMass)
vtIonL = math.sqrt(ionTempL/ionMass)
-- plasma frequency and Debye length on left
wpeL = math.sqrt(elcCharge^2*nL/(epsilon0*elcMass))
wpiL = math.sqrt(ionCharge^2*nL/(epsilon0*ionMass))
lambdaDL = vtElcL/wpeL

-- domain size and simulation time
LX = 100*lambdaDL

vlasovApp = Vlasov.App {
   logToFile = true,

   tEnd = 10.0/wpiL, -- end time
   nFrame = 20, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {LX}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- electrons
   elc = Vlasov.Species {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-8.0*vtElcL},
      upper = {8.0*vtElcL},
      cells = {32},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local sloc = 0.5*LX
	 local fv = nL/math.sqrt(2*math.pi*elcTempL/elcMass)*math.exp(-elcMass*v^2/(2.0*elcTempL))
         if x>sloc then
            fv = nR/math.sqrt(2*math.pi*elcTempR/elcMass)*math.exp(-elcMass*v^2/(2.0*elcTempR))
         end
	 return fv
      end,
      -- boundary conditions
      bcx = { Vlasov.Species.bcOpen, Vlasov.Species.bcOpen },
      
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" }
   },
   -- protons
   ion = Vlasov.Species {
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-16.0*vtIonL},
      upper = {16.0*vtIonL},
      cells = {64},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local sloc = 0.5*LX
	 local fv = nL/math.sqrt(2*math.pi*ionTempL/ionMass)*math.exp(-ionMass*v^2/(2.0*ionTempL))
         if x>sloc then
            fv = nR/math.sqrt(2*math.pi*ionTempR/ionMass)*math.exp(-ionMass*v^2/(2.0*ionTempR))
         end
	 return fv
      end,
      -- boundary conditions
      bcx = { Vlasov.Species.bcOpen, Vlasov.Species.bcOpen },
      
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" }
   },

   -- field solver
   field = Vlasov.EmField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      bcx = { Vlasov.EmField.bcOpen, Vlasov.EmField.bcOpen },
      
      evolve = true, -- evolve field?
   },
}
-- run application
vlasovApp:run()
