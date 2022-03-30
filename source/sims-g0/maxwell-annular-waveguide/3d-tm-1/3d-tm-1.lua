-- Gkyl ------------------------------------------------------------------------
local Moments = require "Moments"
local ffi = require "ffi"

local r0, r1 = 2.0, 5.0 -- inner and outer radii
local m = 4 -- aziumuthal mode number
local n = 1 -- axial mode number
local Lz = 5.0 -- lenght of cylinder
local kn = 2*math.pi*n/Lz -- axial mode
local w = 2.73598725136604 -- frequency of mode

local tperiod = 2*math.pi/w

-- for Bessel functions
ffi.cdef [[
  double jn(int, double);
  double yn(int, double);
]]

-- create app
app = Moments.App {
   tEnd = 2*tperiod, -- end time
   nFrame = 1, -- number of output frame
   lower = {r0, 0.0, 0.0}, -- lower left corner
   upper = {r1, 2*math.pi, Lz}, -- upper right corner
   cells = {48, 48*2, 48*2}, -- number of cells
   cflFrac = 0.9,

   mapc2p = function(t, xn)
      local r, th, z = xn[1], xn[2], xn[3]
      return r*math.cos(th), r*math.sin(th), z
   end,

   periodicDirs = { 2, 3 },

   -- electrons
   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,      

      init = function (t, xn)
	 local r, phi, z = xn[1], xn[2], xn[3]
  
	 local a = 1.0
	 local wkn = math.sqrt(w^2-kn^2)
	 local b = -a*ffi.C.jn(m,r0*wkn)/ffi.C.yn(m,r0*wkn)
	 local Ez_r = a*ffi.C.jn(m,r*wkn) + b*ffi.C.yn(m,r*wkn)
	 local Ez = Ez_r*math.cos(m*phi)*math.cos(kn*z)

	 return 0.0, 0.0, Ez, 0.0, 0.0, 0.0
      end,

      bcx = { Moments.Field.bcPEC, Moments.Field.bcPEC },
   },
}
-- run application
app:run()
