-- load app code
local App = dofile("../code/anisopoisson.lua")

local x0, y0 = 0.0, 0.0 -- location of O
local zet = 1e9 -- Dpar/Dperp
local Dpar = 1.0 -- parallel heat conduction
local Dperp = Dpar/zet -- perpendicular heat conduction

local bx = function(x, y)
   return -(y-y0)/math.sqrt((x-x0)^2+(y-y0)^2)
end
local by = function(x, y)
   return (x-x0)/math.sqrt((x-x0)^2+(y-y0)^2)
end

aPoisson = App {
   polyOrder = 1,
   cflFrac = 0.9,
   lower = {-0.5, -0.5},
   upper = {0.5, 0.5},
   cells = {8, 8},

   bcLower = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} },
   bcUpper = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} },

   -- diffusion coefficients
   Dxx = function (t, z)
      local x, y = z[1], z[2]
      return Dpar*bx(x,y)^2 + Dperp*by(x,y)^2
   end,
   Dyy = function (t, z)
      local x, y = z[1], z[2]
      return Dperp*bx(x,y)^2 + Dpar*by(x,y)^2
   end,
   Dxy = function (t, z)
      local x, y = z[1], z[2]
      return (Dpar-Dperp)*bx(x,y)*by(x,y)
   end,

   -- source (RHS)
   src = function (t, z)
      local x, y = z[1], z[2]
      return (8*Dperp*x^2*y^4)/(y^4+2*x^2*y^2+x^4)-(Dperp*y^4)/(y^4+2*x^2*y^2+x^4)+(8*Dperp*x^4*y^2)/(y^4+2*x^2*y^2+x^4)-(2*Dperp*x^2*y^2)/(y^4+2*x^2*y^2+x^4)-(Dperp*x^4)/(y^4+2*x^2*y^2+x^4)+(2*Dperp*y^2)/(2*y^2+2*x^2)+(2*Dperp*x^2)/(2*y^2+2*x^2)-(2*Dpar*y^4)/(y^2+x^2)-(24*Dperp*x^2*y^2)/(y^2+x^2)+(12*Dpar*x^2*y^2)/(y^2+x^2)+(Dperp*y^2)/(y^2+x^2)-(2*Dpar*x^4)/(y^2+x^2)+(Dperp*x^2)/(y^2+x^2)
   end,

   -- optional exact solution (App will compute projection)
   sol = function (t, z)
      local x, y = z[1], z[2]
      return (x-1/2)*(x+1/2)*(y-1/2)*(y+1/2)
   end
}
aPoisson()
