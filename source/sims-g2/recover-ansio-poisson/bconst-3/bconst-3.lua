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
   cells = {8, 8},

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
      return (-20.0*Dperp*by^2*x^3*y^5)-20.0*Dpar*bx^2*x^3*y^5+3.0*Dperp*by^2*x^2*y^5+3.0*Dpar*bx^2*x^2*y^5+1.875*Dperp*by^2*x*y^5+1.875*Dpar*bx^2*x*y^5-0.15625*Dperp*by^2*y^5-0.15625*Dpar*bx^2*y^5+50.0*Dperp*bx*by*x^4*y^4-50.0*Dpar*bx*by*x^4*y^4-5.0*Dperp*by^2*x^3*y^4-10.0*Dperp*bx*by*x^3*y^4+10.0*Dpar*bx*by*x^3*y^4-5.0*Dpar*bx^2*x^3*y^4+0.75*Dperp*by^2*x^2*y^4-9.375*Dperp*bx*by*x^2*y^4+9.375*Dpar*bx*by*x^2*y^4+0.75*Dpar*bx^2*x^2*y^4+0.46875*Dperp*by^2*x*y^4+1.5625*Dperp*bx*by*x*y^4-1.5625*Dpar*bx*by*x*y^4+0.46875*Dpar*bx^2*x*y^4-0.0390625*Dperp*by^2*y^4+0.15625*Dperp*bx*by*y^4-0.15625*Dpar*bx*by*y^4-0.0390625*Dpar*bx^2*y^4-20.0*Dpar*by^2*x^5*y^3-20.0*Dperp*bx^2*x^5*y^3+5.0*Dpar*by^2*x^4*y^3+10.0*Dperp*bx*by*x^4*y^3-10.0*Dpar*bx*by*x^4*y^3+5.0*Dperp*bx^2*x^4*y^3+6.25*Dperp*by^2*x^3*y^3+6.25*Dpar*by^2*x^3*y^3-2.0*Dperp*bx*by*x^3*y^3+2.0*Dpar*bx*by*x^3*y^3+6.25*Dperp*bx^2*x^3*y^3+6.25*Dpar*bx^2*x^3*y^3-0.9375*Dperp*by^2*x^2*y^3-1.5625*Dpar*by^2*x^2*y^3-1.875*Dperp*bx*by*x^2*y^3+1.875*Dpar*bx*by*x^2*y^3-1.5625*Dperp*bx^2*x^2*y^3-0.9375*Dpar*bx^2*x^2*y^3-0.5859375*Dperp*by^2*x*y^3-0.3125*Dpar*by^2*x*y^3+0.3125*Dperp*bx*by*x*y^3-0.3125*Dpar*bx*by*x*y^3-0.3125*Dperp*bx^2*x*y^3-0.5859375*Dpar*bx^2*x*y^3+0.048828125*Dperp*by^2*y^3+0.078125*Dpar*by^2*y^3+0.03125*Dperp*bx*by*y^3-0.03125*Dpar*bx*by*y^3+0.078125*Dperp*bx^2*y^3+0.048828125*Dpar*bx^2*y^3-3.0*Dpar*by^2*x^5*y^2-3.0*Dperp*bx^2*x^5*y^2+0.75*Dpar*by^2*x^4*y^2-9.375*Dperp*bx*by*x^4*y^2+9.375*Dpar*bx*by*x^4*y^2+0.75*Dperp*bx^2*x^4*y^2+1.5625*Dperp*by^2*x^3*y^2+0.9375*Dpar*by^2*x^3*y^2+1.875*Dperp*bx*by*x^3*y^2-1.875*Dpar*bx*by*x^3*y^2+0.9375*Dperp*bx^2*x^3*y^2+1.5625*Dpar*bx^2*x^3*y^2-0.234375*Dperp*by^2*x^2*y^2-0.234375*Dpar*by^2*x^2*y^2+1.7578125*Dperp*bx*by*x^2*y^2-1.7578125*Dpar*bx*by*x^2*y^2-0.234375*Dperp*bx^2*x^2*y^2-0.234375*Dpar*bx^2*x^2*y^2-0.146484375*Dperp*by^2*x*y^2-0.046875*Dpar*by^2*x*y^2-0.29296875*Dperp*bx*by*x*y^2+0.29296875*Dpar*bx*by*x*y^2-0.046875*Dperp*bx^2*x*y^2-0.146484375*Dpar*bx^2*x*y^2+0.01220703125*Dperp*by^2*y^2+0.01171875*Dpar*by^2*y^2-0.029296875*Dperp*bx*by*y^2+0.029296875*Dpar*bx*by*y^2+0.01171875*Dperp*bx^2*y^2+0.01220703125*Dpar*bx^2*y^2+1.875*Dpar*by^2*x^5*y+1.875*Dperp*bx^2*x^5*y-0.46875*Dpar*by^2*x^4*y-1.5625*Dperp*bx*by*x^4*y+1.5625*Dpar*bx*by*x^4*y-0.46875*Dperp*bx^2*x^4*y-0.3125*Dperp*by^2*x^3*y-0.5859375*Dpar*by^2*x^3*y+0.3125*Dperp*bx*by*x^3*y-0.3125*Dpar*bx*by*x^3*y-0.5859375*Dperp*bx^2*x^3*y-0.3125*Dpar*bx^2*x^3*y+0.046875*Dperp*by^2*x^2*y+0.146484375*Dpar*by^2*x^2*y+0.29296875*Dperp*bx*by*x^2*y-0.29296875*Dpar*bx*by*x^2*y+0.146484375*Dperp*bx^2*x^2*y+0.046875*Dpar*bx^2*x^2*y+0.029296875*Dperp*by^2*x*y+0.029296875*Dpar*by^2*x*y-0.048828125*Dperp*bx*by*x*y+0.048828125*Dpar*bx*by*x*y+0.029296875*Dperp*bx^2*x*y+0.029296875*Dpar*bx^2*x*y-0.00244140625*Dperp*by^2*y-0.00732421875*Dpar*by^2*y-0.0048828125*Dperp*bx*by*y+0.0048828125*Dpar*bx*by*y-0.00732421875*Dperp*bx^2*y-0.00244140625*Dpar*bx^2*y+0.15625*Dpar*by^2*x^5+0.15625*Dperp*bx^2*x^5-0.0390625*Dpar*by^2*x^4+0.15625*Dperp*bx*by*x^4-0.15625*Dpar*bx*by*x^4-0.0390625*Dperp*bx^2*x^4-0.078125*Dperp*by^2*x^3-0.048828125*Dpar*by^2*x^3-0.03125*Dperp*bx*by*x^3+0.03125*Dpar*bx*by*x^3-0.048828125*Dperp*bx^2*x^3-0.078125*Dpar*bx^2*x^3+0.01171875*Dperp*by^2*x^2+0.01220703125*Dpar*by^2*x^2-0.029296875*Dperp*bx*by*x^2+0.029296875*Dpar*bx*by*x^2+0.01220703125*Dperp*bx^2*x^2+0.01171875*Dpar*bx^2*x^2+0.00732421875*Dperp*by^2*x+0.00244140625*Dpar*by^2*x+0.0048828125*Dperp*bx*by*x-0.0048828125*Dpar*bx*by*x+0.00244140625*Dperp*bx^2*x+0.00732421875*Dpar*bx^2*x-6.103515625e-4*Dperp*by^2-6.103515625e-4*Dpar*by^2+4.8828125e-4*Dperp*bx*by-4.8828125e-4*Dpar*bx*by-6.103515625e-4*Dperp*bx^2-6.103515625e-4*Dpar*bx^2
   end,

   -- optional exact solution (App will compute projection)
   sol = function (t, z)
      local x, y = z[1], z[2]
      return (x-1/2)*(x+1/2)*(x+1/4)*(x-1/4)^2*(y-1/2)*(y+1/2)*(y+1/4)^2*(y-1/4)
   end  
}
aPoisson()
