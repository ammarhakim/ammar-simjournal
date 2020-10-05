-- load app code
local App = dofile("../code/sts-diffusion.lua")

diffusion = App {
   polyOrder = 2,
   lower = {-1.0, -1.0},
   upper = {1.0, 1.0},
   cells = {128, 128},
   errEps = 1e-8,
   factor = 2000,
   extraStages = 15,
   cflFrac = 0.8,
   
   -- initial conditions
   init = function (t, xn)
      local x, y = xn[1], xn[2]
      return 0.0
   end,

   source = function (t, xn)
      local x, y = xn[1], xn[2]
      return math.exp(-10*(x^2+y^2))
   end,
}
diffusion()
