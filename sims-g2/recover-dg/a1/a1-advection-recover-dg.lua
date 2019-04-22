-- load app code
local App = dofile("../code/advection.lua")

advection = App {
   polyOrder = 2, -- polynomial order
   cflFrac = 0.5, -- cflFrac (defaults to 1.0)
   extents = {-math.pi, math.pi}, -- domain size
   nCell = 8, -- number of cells
   tEnd = 2*math.pi*80, -- time

   -- initial conditions
   init = function(t, xn)
      return math.sin(xn[1])
   end,
}
advection()
