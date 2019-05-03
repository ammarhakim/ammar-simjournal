-- load app code
local App = dofile("../code/advection.lua")

advection = App {
   polyOrder = 1, -- polynomial order
   cflFrac = 4.0/128, -- cflFrac (defaults to 1.0)
   extents = {-1, 1}, -- domain size
   nCell = 64, -- number of cells
   tEnd = 2.0, -- time

   -- initial conditions
   init = function(t, xn)
      return 1+math.exp(-xn[1]^2/(2*0.2^2))
   end,
}
advection()
