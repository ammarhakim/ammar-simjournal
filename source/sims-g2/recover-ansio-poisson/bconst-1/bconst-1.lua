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
   cells = {4, 4},

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
      local bx, by = math.cos(theta), math.sin(theta)
      return (-6.0*Dperp*by^2*x*y^3)-6.0*Dpar*bx^2*x*y^3+0.5*Dperp*by^2*y^3+0.5*Dpar*bx^2*y^3+18.0*Dperp*bx*by*x^2*y^2-18.0*Dpar*bx*by*x^2*y^2-1.5*Dperp*by^2*x*y^2-3.0*Dperp*bx*by*x*y^2+3.0*Dpar*bx*by*x*y^2-1.5*Dpar*bx^2*x*y^2+0.125*Dperp*by^2*y^2-1.5*Dperp*bx*by*y^2+1.5*Dpar*bx*by*y^2+0.125*Dpar*bx^2*y^2-6.0*Dpar*by^2*x^3*y-6.0*Dperp*bx^2*x^3*y+1.5*Dpar*by^2*x^2*y+3.0*Dperp*bx*by*x^2*y-3.0*Dpar*bx*by*x^2*y+1.5*Dperp*bx^2*x^2*y+1.5*Dperp*by^2*x*y+1.5*Dpar*by^2*x*y-0.5*Dperp*bx*by*x*y+0.5*Dpar*bx*by*x*y+1.5*Dperp*bx^2*x*y+1.5*Dpar*bx^2*x*y-0.125*Dperp*by^2*y-0.375*Dpar*by^2*y-0.25*Dperp*bx*by*y+0.25*Dpar*bx*by*y-0.375*Dperp*bx^2*y-0.125*Dpar*bx^2*y-0.5*Dpar*by^2*x^3-0.5*Dperp*bx^2*x^3+0.125*Dpar*by^2*x^2-1.5*Dperp*bx*by*x^2+1.5*Dpar*bx*by*x^2+0.125*Dperp*bx^2*x^2+0.375*Dperp*by^2*x+0.125*Dpar*by^2*x+0.25*Dperp*bx*by*x-0.25*Dpar*bx*by*x+0.125*Dperp*bx^2*x+0.375*Dpar*bx^2*x-0.03125*Dperp*by^2-0.03125*Dpar*by^2+0.125*Dperp*bx*by-0.125*Dpar*bx*by-0.03125*Dperp*bx^2-0.03125*Dpar*bx^2 
   end,

   -- optional exact solution (App will compute projection)
   sol = function (t, z)
      local x, y = z[1], z[2]
      return (x-1/2)*(x+1/2)*(x-1/4)*(y-1/2)*(y+1/2)*(y+1/4)
   end  
}
aPoisson()
