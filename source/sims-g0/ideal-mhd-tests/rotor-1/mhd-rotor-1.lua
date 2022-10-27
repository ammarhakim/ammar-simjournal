-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"
local ffi = require "ffi"

gasGamma = 1.4

-- create app
mhdApp = Moments.App {
   logToFile = true,

   tEnd = 0.15, -- end time
   nFrame = 10, -- number of output frame
   lower = { 0.0, 0.0 }, -- lower left corner
   upper = { 1.0, 1.0 }, -- upper right corner
   cells = { 400, 400 }, -- number of cells

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
	 local r = math.sqrt((x-0.5)^2 + (y-0.5)^2)
	 
	 local gas_gamma = gasGamma
	 local rho, vx, vy, vz

	 local v0 = 2.0
	 local r0, r1 = 0.1, 0.115
	 local fr = (r1-r)/(r1-r0)

	 local vz = 0.0
	 local p, Bx, By, Bz = 1.0, 5/math.sqrt(4*math.pi), 0.0, 0.0
	 if r<=r0 then
	    rho = 10.0
	    vx = -v0*(y-0.5)/r0
	    vy = v0*(x-0.5)/r0
	 elseif r>r0 and r<=r1 then
	    rho = 1.0 + 9*fr
	    vx = -fr*v0*(y-0.5)/r
	    vy = fr*v0*(x-0.5)/r
	 else
	    rho = 1.0
	    vx = 0.0
	    vy = 0.0
	 end

	 -- gkyl_mhd_cons_vars is a C function, expecting 0-indexed C arrays
	 local pv = ffi.new("double[8]", rho, vx, vy, vz, p, Bx, By, Bz)
	 local q = ffi.new("double[9]");

	 Moments.gkyl_mhd_cons_vars(gasGamma, pv, q)

	 return q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7]
      end,

      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
      bcy = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },
}
-- run application
mhdApp:run()
