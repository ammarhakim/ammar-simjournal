-- load app code
local App = dofile("../code/anisopoisson.lua")

local x0, y0 = 1.0, 0.0 -- location of O-point
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
   cells = {16, 16},

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
      local xc, yc = 0.2, 0.0
      return math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
   end,
}
aPoisson()
