-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"
local ffi = require "ffi"

gasGamma = 5.0/3.0

-- create app
mhdApp = Moments.App {
   logToFile = true,

   tEnd = 0.06, -- end time
   nFrame = 10, -- number of output frame
   lower = { 0.0, 1.0 }, -- lower left corner
   upper = { 0.0, 1.0 }, -- upper right corner
   cells = { 200, 200 }, -- number of cells

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.MHD {
	 gasGamma = gasGamma,
	 rp_type = "roe", -- one of "roe", "hlld", "lax"
	 divergenceConstraint = "eight_waves", -- one of "none", "eight_waves", "glm"
      },

      -- initial conditions
      init = function (t, xn)
	 local x, y = xn[1], xn[2]
	 
	 local gas_gamma = gasGamma

	 local rho = 25/(36*math.pi)
	 local p = 5/(12*math.pi)
	 local vx = math.sin(2*math.pi*y)
	 local vy = -math.sin(2*math.pi*x)
	 local vz = 0
	 local B0 = 1/math.sqrt(4*math.pi)
	 local Bx = B0*math.sin(2*math.pi*y)
	 local By = B0*math.sin(4*math.pi*x)
	 local Bz = 0

	 -- gkyl_mhd_cons_vars is a C function, expecting 0-indexed C arrays
	 local pv = ffi.new("double[8]", rho, vx, vy, vz, p, Bx, By, Bz)
	 local q = ffi.new("double[9]");

	 Moments.gkyl_mhd_cons_vars(gasGamma, pv, q)

	 return q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7]
      end,
   },
}
-- run application
mhdApp:run()
