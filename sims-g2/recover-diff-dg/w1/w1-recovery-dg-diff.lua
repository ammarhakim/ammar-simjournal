-- load app code
local App = dofile("../code/diffusion.lua")

-- diffusion coefficients
Dxx = 1.0
Dyy = 1.0
Dxy = 0.5

diffusion = App {
   polyOrder = 1,
   cflFrac = 0.9,
   lower = {0.0, 0.0},
   upper = {2*math.pi, 2*math.pi},
   cells = {8, 8},
   tEnd = 40.0,
   numFrames = 40,

   -- diffusion coefficient
   D = { Dxx = Dxx, Dyy = Dyy, Dxy = Dxy },
   
   -- initial conditions
   init = function(t, xn)
      local x, y = xn[1], xn[2]
      return 0.0
   end,

   -- source term
   source = function (t, xn)
      local x, y = xn[1], xn[2]
      local kx, ky = 1.0, 1.0
      local cos, sin = math.cos, math.sin
      return 2*Dxy*cos(x)*sin(y)+Dyy*sin(x)*cos(y)+Dxx*sin(x)*cos(y)
   end,
}
diffusion()
