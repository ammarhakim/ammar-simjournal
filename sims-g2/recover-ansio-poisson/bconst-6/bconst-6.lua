-- load app code
local App = dofile("../code/anisopoisson.lua")

local x0, y0 = 0.0, 0.0 -- location of O
local zet = 1e9 -- Dpar/Dperp
local Dpar = 1.0 -- parallel heat conduction
local Dperp = Dpar/zet -- perpendicular heat conduction
local theta = 15.0*math.pi/180

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
      local bx, by = math.cos(theta), math.sin(theta)      
      return Dpar*bx^2 + Dperp*by^2
   end,
   Dyy = function (t, z)
      local x, y = z[1], z[2]
      local bx, by = math.cos(theta), math.sin(theta)      
      return Dperp*bx^2 + Dpar*by^2
   end,
   Dxy = function (t, z)
      local x, y = z[1], z[2]
      local bx, by = math.cos(theta), math.sin(theta)
      return (Dpar-Dperp)*bx*by
   end,

   -- source (RHS)
   src = function (t, z)
      local x, y = z[1], z[2]
      local xc, yc = -0.25, -0.25
      return math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
   end,
}
aPoisson()
