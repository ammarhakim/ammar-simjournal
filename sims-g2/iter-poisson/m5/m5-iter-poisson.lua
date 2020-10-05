-- load app code
local App = dofile("../code/sts-diffusion.lua")

diffusion = App {
   polyOrder = 2,
   lower = {-1.0, -1.0},
   upper = {1.0, 1.0},
   cells = {64, 64},
   errEps = 1e-8,
   cflFrac = 0.85,
   factor = 600,
   extraStages = 6,
   stepper = 'RKL1',
   
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
