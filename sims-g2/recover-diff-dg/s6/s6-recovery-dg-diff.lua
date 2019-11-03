-- load app code
local App = dofile("../code/diffusion.lua")

diffusion = App {
   polyOrder = 1,
   cflFrac = 0.9,
   lower = {0.0, 0.0},
   upper = {2*math.pi, 2*math.pi},
   cells = {12, 12},

   -- diffusion coefficient
   D = { Dxx = 1.0, Dyy = 1.0, Dxy = 0.9 },
   
   tEnd = 0.5,
   numFrames = 1,

   -- initial conditions
   init = function(t, xn)
      local x, y = xn[1]-math.pi, xn[2]-math.pi
      return math.exp(-2*(x^2+y^2))
   end,
}
diffusion()
