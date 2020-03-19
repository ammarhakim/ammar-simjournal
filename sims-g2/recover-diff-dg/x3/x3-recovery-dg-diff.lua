-- load app code
local App = dofile("../code/diffusion.lua")

-- magnetic field angle
angle = 15*math.pi/180
bx = math.cos(angle)
by = math.sin(angle)
alpha = 10^9 -- Dpar/Dperp
Dpar = 1.0 -- parallel heat conduction

Dperp = Dpar/alpha -- perpendicular heat conduction

-- diffusion coefficients
Dxx = Dpar*bx^2 + Dperp*by^2
Dyy = Dperp*bx^2 + Dpar*by^2
Dxy = (Dpar-Dperp)*bx*by

diffusion = App {
   polyOrder = 2,
   cflFrac = 0.7,
   lower = {0.0, 0.0},
   upper = {2*math.pi, 2*math.pi},
   cells = {8, 8},
   tEnd = 400.0,
   numFrames = 100,

   -- diffusion coefficient
   D = { Dxx = Dxx, Dyy = Dyy, Dxy = Dxy },
   
   -- initial conditions
   init = function(t, xn)
      local x, y = xn[1], xn[2]
      local cos, sin = math.cos, math.sin
      return -cos(x)*sin(y)
   end,

   -- exact solution
   exact = function(t, xn)
      local x, y = xn[1], xn[2]
      local cos, sin = math.cos, math.sin
      return -cos(x)*sin(y)
   end,   

   -- source term
   source = function (t, xn)
      local x, y = xn[1], xn[2]
      local cos, sin = math.cos, math.sin
      return -Dyy*cos(x)*sin(y)-Dxx*cos(x)*sin(y)-2*Dxy*sin(x)*cos(y)
   end,
}
diffusion()
