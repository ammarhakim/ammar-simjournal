-- load app code
local App = dofile("../code/advection.lua")

advection = App {
   polyOrder = 2, -- polynomial order
   cflFrac = 8.0/128/4.0, -- cflFrac (defaults to 1.0)
   extents = {-math.pi, math.pi}, -- domain size
   nCell = 32, -- number of cells
   tEnd = 2*math.pi, -- time

   -- initial conditions
   init = function(t, xn)
      return math.sin(xn[1])
   end,
}
advection()
