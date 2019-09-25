-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"

charge = 100.0
gasGamma = 5.0/3.0
elcCharge = -charge
elcMass = 1/1836.2
ionCharge = charge
ionMass = 1.0

momentApp = Moments.App {
   logToFile = true,

   tEnd = 10.0,
   nFrame = 200,
   lower = {-2.0},
   upper = {2.0},
   cells = {8000},
   timeStepper = "fvDimSplit",

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      -- initial conditions
      init = function (t, xn)
	 local x, sloc = xn[1], 0.0
	 if (x<sloc) then
	    return 1.0*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
	 else
	    return 0.125*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
	 end
      end,
      evolve = true, -- evolve species?
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      -- initial conditions
      init = function (t, xn)
	 local x, sloc = xn[1], 0.0
	 if (x<sloc) then
	    return 1.0, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
	 else
	    return 0.125, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
	 end
      end,
      evolve = true, -- evolve species?
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x, sloc = xn[1], 0.0
	 if (x<sloc) then
	    return 0.0, 0.0, 0.0, 0.75e-2, 0.0, 1.0e-2
	 else
	    return 0.0, 0.0, 0.0, 0.75e-2, 0.0, -1.0e-2
	 end
      end,
      evolve = true, -- evolve field?
      bcx = { Moments.Field.bcCopy, Moments.Field.bcCopy },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
   },   

}
-- run application
momentApp:run()
