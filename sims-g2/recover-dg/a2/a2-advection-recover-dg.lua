-- load app code
local App = dofile("../code/advection.lua")

advection = App {
   polyOrder = 1, -- polynomial order
   cflFrac = 1.0/8, -- cflFrac (defaults to 1.0)
   extents = {-math.pi, math.pi}, -- domain size
   nCell = 8, -- number of cells
   tEnd = 2*math.pi, -- time

   -- initial conditions
   init = function(t, xn)
      return math.sin(xn[1])
   end,
}
advection()
