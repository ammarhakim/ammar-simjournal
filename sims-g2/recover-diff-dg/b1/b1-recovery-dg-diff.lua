-- load app code
local App = dofile("../code/diffusion.lua")

diffusion = App {
   polyOrder = 1,
   cflFrac = 0.9,
   lower = {0.0, 0.0},
   upper = {2*math.pi, 2*math.pi},
   cells = {8, 8},
   tEnd = 1.0,
   numFrames = 1,
   useFivePointStencil = true,

   -- diffusion coefficient
   D = { Dxx = 1.0, Dyy = 1.0, Dxy = 1/6.0 },
   
   -- initial conditions
   init = function(t, xn)
      local x, y = xn[1], xn[2]
      local kx, ky = 1.0, 1.0
      return math.cos(x*kx + y*ky)
   end,
}
diffusion()
