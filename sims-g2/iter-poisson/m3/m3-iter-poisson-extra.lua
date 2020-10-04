-- load app code
local App = dofile("../code/sts-diffusion.lua")

diffusion = App {
   polyOrder = 1,
   lower = {-1.0, -1.0},
   upper = {1.0, 1.0},
   cells = {128, 128},
   errEps = 1e-8,
   factor = 1000,
   extraStages = 9,
   stepper = 'RKL1',
   extrapolateInterval = 2,
  
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
