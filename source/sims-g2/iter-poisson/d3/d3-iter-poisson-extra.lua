-- load app code
local App = dofile("../code/sts-diffusion.lua")

diffusion = App {
   polyOrder = 1,
   lower = {0.0, 0.0},
   upper = {2*math.pi, 2*math.pi},
   cells = {32, 32},
   errEps = 1e-8,
   factor = 160,
   extraStages = 4,
   stepper = 'RKL1',
   extrapolateInterval = 2,
   
   -- initial conditions
   init = function (t, xn)
      local x, y = xn[1], xn[2]
      return 0.0
   end,

   source = function (t, xn)
      local x, y = xn[1], xn[2]
      local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local t1, t2 = 0.0, 0.0
      local f = 0.0
      for m = 0,2 do
	 for n = 0,2 do
	    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
	    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
	    f = f + -(m*m+n*n)*(t1+t2)
	 end
      end
      return -f/50.0
   end,

   exact = function (t, xn)
      local x, y = xn[1], xn[2]
      local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local t1, t2 = 0.0, 0.0
      local f = 0.0

      for m = 0,2 do
	 for n = 0,2 do
	    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
	    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
	    f = f + (t1+t2)
	 end
      end
      return f/50.0
   end,
}
diffusion()
