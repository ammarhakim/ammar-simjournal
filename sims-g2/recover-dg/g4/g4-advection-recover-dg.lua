-- load app code
local App = dofile("../code/advection.lua")

advection = App {
   polyOrder = 2, -- polynomial order
   cflFrac = 1.0/128, -- cflFrac (defaults to 1.0)
   extents = {-1, 1}, -- domain size
   nCell = 12, -- number of cells
   tEnd = 2.0, -- time

   -- initial conditions
   init = function(t, xn)
      return 1.0+math.exp(-xn[1]^2/(2*0.2^2))
   end,
}
advection()
