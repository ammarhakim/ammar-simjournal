-- Gkyl ------------------------------------------------------------------------
--
-- Original by Petr Cagas. Modified by AH. Not use in the actual
-- simulations
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local x0, y0 = 0.0, 0.0 -- location of O/X-point
local zet = 1e9 -- Dpar/Dperp
local Dpar = 1.0 -- parallel heat conduction
local Dperp = Dpar/zet -- perpendicular heat conduction

local numQuad = 8
local maxPolyOrder = 1
local numCellListX = {8}
local numCellListY = {8}
-- one of "bi-quad", "bi-quad-noxy", "bi-quad-xyc", "bi-cubic",
-- "bi-cubic-shift", "gauss", "xygauss", "blob", "2blobs", "4xyblobs",
-- "constant", "pulse"
local srcFuncType = "pulse"

-- for "pulse" set cell in which source is pulsed
local pList = { {2, 2}, {7, 2} }
local pVals = { 1.0, 0.0 }
local NX_G, NY_G = 0, 0 -- number of cells in grid (needed in "pulse")

-- boundary conditions
bcLower = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} }
bcUpper = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} }

local bx = function(x, y)
   return -(y-y0)/math.sqrt((x-x0)^2+(y-y0)^2)
end
local by = function(x, y)
   return (x-x0)/math.sqrt((x-x0)^2+(y-y0)^2)
end

-- diffusion coefficients
local DxxFn = function(t, z)
   local x, y = z[1], z[2]
   return Dpar*bx(x,y)^2 + Dperp*by(x,y)^2
end
local DyyFn = function(t, z)
   local x, y = z[1], z[2]
   return Dperp*bx(x,y)^2 + Dpar*by(x,y)^2
end
local DxyFn = function(t, z)
   local x, y = z[1], z[2]
   if srcFuncType == "bi-quad-noxy" then
      return 0
   elseif srcFuncType == "bi-quad-xyc" then
      return 0.1
   else
      return (Dpar-Dperp)*bx(x,y)*by(x,y)
   end
end

local solFunc = function(x,y) return 0 end
local srcFunc = function(x,y) return 0 end

if srcFuncType == "bi-quad" then
   solFunc = function (x,y)
      return (x-1/2)*(x+1/2)*(y-1/2)*(y+1/2)
   end
   
   srcFunc = function (x,y)
      return (8*Dperp*x^2*y^4)/(y^4+2*x^2*y^2+x^4)-(Dperp*y^4)/(y^4+2*x^2*y^2+x^4)+(8*Dperp*x^4*y^2)/(y^4+2*x^2*y^2+x^4)-(2*Dperp*x^2*y^2)/(y^4+2*x^2*y^2+x^4)-(Dperp*x^4)/(y^4+2*x^2*y^2+x^4)+(2*Dperp*y^2)/(2*y^2+2*x^2)+(2*Dperp*x^2)/(2*y^2+2*x^2)-(2*Dpar*y^4)/(y^2+x^2)-(24*Dperp*x^2*y^2)/(y^2+x^2)+(12*Dpar*x^2*y^2)/(y^2+x^2)+(Dperp*y^2)/(y^2+x^2)-(2*Dpar*x^4)/(y^2+x^2)+(Dperp*x^2)/(y^2+x^2) 
   end

elseif srcFuncType == "bi-quad-noxy" then

   solFunc = function (x,y)
      return (x-1/2)*(x+1/2)*(y-1/2)*(y+1/2)
   end
   
   srcFunc = function (x,y)
      return (4.0*Dperp*x^2*y^4)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dpar*x^2*y^4)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*y^4)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dperp*x^4*y^2)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dpar*x^4*y^2)/(y^4+2.0*x^2*y^2+x^4)-(2.0*Dpar*x^2*y^2)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*x^4)/(y^4+2.0*x^2*y^2+x^4)+(Dperp*y^2)/(2.0*y^2+2.0*x^2)+(Dpar*y^2)/(2.0*y^2+2.0*x^2)+(Dperp*x^2)/(2.0*y^2+2.0*x^2)+(Dpar*x^2)/(2.0*y^2+2.0*x^2)-(2.0*Dpar*y^4)/(y^2+x^2)-(12.0*Dperp*x^2*y^2)/(y^2+x^2)+(Dperp*y^2)/(y^2+x^2)-(2.0*Dpar*x^4)/(y^2+x^2)+(Dperp*x^2)/(y^2+x^2) 
   end

elseif srcFuncType == "bi-quad-xyc" then

   solFunc = function (x,y)
      return (x-1/2)*(x+1/2)*(y-1/2)*(y+1/2)
   end
   
   srcFunc = function (x,y)
      return (4.0*Dperp*x^2*y^4)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dpar*x^2*y^4)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*y^4)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dperp*x^4*y^2)/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dpar*x^4*y^2)/(y^4+2.0*x^2*y^2+x^4)-(2.0*Dpar*x^2*y^2)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*x^4)/(y^4+2.0*x^2*y^2+x^4)+(Dperp*y^2)/(2.0*y^2+2.0*x^2)+(Dpar*y^2)/(2.0*y^2+2.0*x^2)+(Dperp*x^2)/(2.0*y^2+2.0*x^2)+(Dpar*x^2)/(2.0*y^2+2.0*x^2)-(2.0*Dpar*y^4)/(y^2+x^2)-(12.0*Dperp*x^2*y^2)/(y^2+x^2)+(Dperp*y^2)/(y^2+x^2)-(2.0*Dpar*x^4)/(y^2+x^2)+(Dperp*x^2)/(y^2+x^2)-0.8*x*y 
   end   

elseif srcFuncType == "bi-cubic" then
   solFunc = function (x,y)
      return (x-1/2)*(x+1/2)*x*(y-1/2)*(y+1/2)*y
   end
   
   srcFunc = function (x,y)
      return (2*Dperp*x*y^3)/(8*y^4+16*x^2*y^2+8*x^4)+(2*Dperp*x^3*y)/(8*y^4+16*x^2*y^2+8*x^4)-(4*Dperp*x*y^5)/(2*y^4+4*x^2*y^2+2*x^4)-(8*Dperp*x^3*y^3)/(2*y^4+4*x^2*y^2+2*x^4)-(4*Dperp*x^5*y)/(2*y^4+4*x^2*y^2+2*x^4)+(12*Dperp*x^3*y^5)/(y^4+2*x^2*y^2+x^4)+(12*Dperp*x^5*y^3)/(y^4+2*x^2*y^2+x^4)-(2*Dperp*x*y)/(16*y^2+16*x^2)+(2*Dpar*x*y)/(16*y^2+16*x^2)-(3*Dperp*x*y)/(8*y^2+8*x^2)+(Dpar*x*y)/(8*y^2+8*x^2)+(4*Dperp*x*y^3)/(4*y^2+4*x^2)-(4*Dpar*x*y^3)/(4*y^2+4*x^2)+(4*Dperp*x^3*y)/(4*y^2+4*x^2)-(4*Dpar*x^3*y)/(4*y^2+4*x^2)+(10*Dperp*x*y^3)/(2*y^2+2*x^2)+(10*Dperp*x^3*y)/(2*y^2+2*x^2)-(6*Dpar*x*y^5)/(y^2+x^2)-(48*Dperp*x^3*y^3)/(y^2+x^2)+(24*Dpar*x^3*y^3)/(y^2+x^2)-(6*Dpar*x^5*y)/(y^2+x^2)
   end

elseif srcFuncType == "gauss" then
   solFunc = function (x,y)
      return math.exp(-(x^2+y^2)/(2*0.1^2))
   end
   srcFunc = function (x,y)
      return (-(200.0*Dperp*y^4*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4))-(399.9999999999999*Dperp*x^2*y^2*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)-(200.0*Dperp*x^4*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)-(9999.999999999996*Dperp*y^4*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(19999.99999999999*Dperp*x^2*y^2*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)+(399.9999999999999*Dperp*y^2*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(9999.999999999996*Dperp*x^4*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)+(399.9999999999999*Dperp*x^2*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)
   end

elseif srcFuncType == "xygauss" then
   solFunc = function (x,y)
      return x*y*math.exp(-(x^2+y^2)/(2*0.1^2))
   end
   srcFunc = function (x,y)
      return (-(200.0*Dperp*x*y^5*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4))-(399.9999999999999*Dperp*x^3*y^3*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dperp*x*y^3*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)-(200.0*Dperp*x^5*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)+(4.0*Dperp*x^3*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^4+2.0*x^2*y^2+x^4)-(9999.999999999996*Dperp*x*y^5*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(19999.99999999999*Dperp*x^3*y^3*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)+(799.9999999999998*Dperp*x*y^3*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(2.842170943040401e-14*Dpar*x*y^3*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(9999.999999999996*Dperp*x^5*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)+(799.9999999999999*Dperp*x^3*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)-(8.0*Dperp*x*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2)+(4.0*Dpar*x*y*2.718281828459045^((-49.99999999999999*y^2)-49.99999999999999*x^2))/(y^2+x^2) 
   end   

elseif srcFuncType == "bi-cubic-shift" then
   solFunc = function (x,y)
      return (x-1/2)*(x+1/2)*(x+1/6)*(y-1/2)*(y+1/2)*(y-1/4)
   end
   
   srcFunc = function (x,y)
      return (Dperp*y^3)/(48.0*y^4+96.0*x^2*y^2+48.0*x^4)+(Dperp*x^2*y)/(48.0*y^4+96.0*x^2*y^2+48.0*x^4)-(1.0*Dperp*x*y^2)/(32.0*y^4+64.0*x^2*y^2+32.0*x^4)-(1.0*Dperp*x^3)/(32.0*y^4+64.0*x^2*y^2+32.0*x^4)+(Dperp*y^4)/(24.0*y^4+48.0*x^2*y^2+24.0*x^4)+(2.0*Dperp*x^2*y^2)/(24.0*y^4+48.0*x^2*y^2+24.0*x^4)+(Dperp*x^4)/(24.0*y^4+48.0*x^2*y^2+24.0*x^4)-(1.0*Dperp*x^2*y^3)/(12.0*y^4+24.0*x^2*y^2+12.0*x^4)-(1.0*Dperp*x^4*y)/(12.0*y^4+24.0*x^2*y^2+12.0*x^4)+(Dperp*x*y^4)/(8.0*y^4+16.0*x^2*y^2+8.0*x^4)+(2.0*Dperp*x*y^3)/(8.0*y^4+16.0*x^2*y^2+8.0*x^4)+(4.0*Dperp*x^3*y^2)/(8.0*y^4+16.0*x^2*y^2+8.0*x^4)+(2.0*Dperp*x^3*y)/(8.0*y^4+16.0*x^2*y^2+8.0*x^4)+(3.0*Dperp*x^5)/(8.0*y^4+16.0*x^2*y^2+8.0*x^4)-(2.0*Dperp*x^2*y^4)/(6.0*y^4+12.0*x^2*y^2+6.0*x^4)-(1.0*Dperp*x^2*y^3)/(6.0*y^4+12.0*x^2*y^2+6.0*x^4)-(2.0*Dperp*x^4*y^2)/(6.0*y^4+12.0*x^2*y^2+6.0*x^4)-(1.0*Dperp*x^4*y)/(6.0*y^4+12.0*x^2*y^2+6.0*x^4)-(1.0*Dperp*y^5)/(4.0*y^4+8.0*x^2*y^2+4.0*x^4)+(Dperp*x*y^4)/(4.0*y^4+8.0*x^2*y^2+4.0*x^4)-(1.0*Dperp*x^2*y^3)/(4.0*y^4+8.0*x^2*y^2+4.0*x^4)+(Dperp*x^3*y^2)/(4.0*y^4+8.0*x^2*y^2+4.0*x^4)+(2.0*Dperp*x^2*y^5)/(3.0*y^4+6.0*x^2*y^2+3.0*x^4)+(2.0*Dperp*x^4*y^3)/(3.0*y^4+6.0*x^2*y^2+3.0*x^4)-(4.0*Dperp*x*y^5)/(2.0*y^4+4.0*x^2*y^2+2.0*x^4)-(3.0*Dperp*x^3*y^4)/(2.0*y^4+4.0*x^2*y^2+2.0*x^4)-(8.0*Dperp*x^3*y^3)/(2.0*y^4+4.0*x^2*y^2+2.0*x^4)-(3.0*Dperp*x^5*y^2)/(2.0*y^4+4.0*x^2*y^2+2.0*x^4)-(4.0*Dperp*x^5*y)/(2.0*y^4+4.0*x^2*y^2+2.0*x^4)+(12.0*Dperp*x^3*y^5)/(y^4+2.0*x^2*y^2+x^4)+(Dperp*x^2*y^5)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*x^3*y^4)/(y^4+2.0*x^2*y^2+x^4)+(12.0*Dperp*x^5*y^3)/(y^4+2.0*x^2*y^2+x^4)+(Dperp*x^4*y^3)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*x^5*y^2)/(y^4+2.0*x^2*y^2+x^4)-(1.0*Dperp*y)/(96.0*y^2+96.0*x^2)+(Dpar*y)/(96.0*y^2+96.0*x^2)+(Dperp*x)/(64.0*y^2+64.0*x^2)-(1.0*Dpar*x)/(64.0*y^2+64.0*x^2)-(2.0*Dperp*y^2)/(48.0*y^2+48.0*x^2)-(1.0*Dperp*y)/(48.0*y^2+48.0*x^2)-(2.0*Dperp*x^2)/(48.0*y^2+48.0*x^2)+(Dperp*x)/(32.0*y^2+32.0*x^2)-(1.0*Dperp*y^2)/(24.0*y^2+24.0*x^2)+(Dperp*x^2*y)/(24.0*y^2+24.0*x^2)-(1.0*Dpar*x^2*y)/(24.0*y^2+24.0*x^2)-(1.0*Dperp*x^2)/(24.0*y^2+24.0*x^2)-(1.0*Dperp*x*y^2)/(16.0*y^2+16.0*x^2)+(Dpar*x*y^2)/(16.0*y^2+16.0*x^2)-(2.0*Dperp*x*y)/(16.0*y^2+16.0*x^2)+(2.0*Dpar*x*y)/(16.0*y^2+16.0*x^2)-(3.0*Dperp*x^3)/(16.0*y^2+16.0*x^2)+(3.0*Dpar*x^3)/(16.0*y^2+16.0*x^2)+(Dpar*y^4)/(12.0*y^2+12.0*x^2)+(Dpar*y^3)/(12.0*y^2+12.0*x^2)+(4.0*Dperp*x^2*y^2)/(12.0*y^2+12.0*x^2)-(2.0*Dpar*x^2*y^2)/(12.0*y^2+12.0*x^2)+(3.0*Dperp*x^2*y)/(12.0*y^2+12.0*x^2)-(1.0*Dpar*x^2*y)/(12.0*y^2+12.0*x^2)+(Dpar*x^4)/(12.0*y^2+12.0*x^2)+(Dperp*y^3)/(8.0*y^2+8.0*x^2)-(1.0*Dpar*y^3)/(8.0*y^2+8.0*x^2)-(3.0*Dperp*x*y^2)/(8.0*y^2+8.0*x^2)-(2.0*Dpar*x*y^2)/(8.0*y^2+8.0*x^2)-(3.0*Dperp*x*y)/(8.0*y^2+8.0*x^2)+(Dpar*x*y)/(8.0*y^2+8.0*x^2)-(6.0*Dperp*x^3)/(8.0*y^2+8.0*x^2)-(1.0*Dpar*x^3)/(8.0*y^2+8.0*x^2)+(2.0*Dperp*x^2*y^2)/(6.0*y^2+6.0*x^2)+(2.0*Dperp*x^2*y)/(6.0*y^2+6.0*x^2)-(1.0*Dpar*x^2*y)/(6.0*y^2+6.0*x^2)+(4.0*Dperp*x*y^3)/(4.0*y^2+4.0*x^2)-(4.0*Dpar*x*y^3)/(4.0*y^2+4.0*x^2)+(2.0*Dperp*y^3)/(4.0*y^2+4.0*x^2)+(3.0*Dperp*x^3*y^2)/(4.0*y^2+4.0*x^2)-(3.0*Dpar*x^3*y^2)/(4.0*y^2+4.0*x^2)-(2.0*Dperp*x*y^2)/(4.0*y^2+4.0*x^2)+(Dpar*x*y^2)/(4.0*y^2+4.0*x^2)+(4.0*Dperp*x^3*y)/(4.0*y^2+4.0*x^2)-(4.0*Dpar*x^3*y)/(4.0*y^2+4.0*x^2)+(Dpar*x^2*y)/(4.0*y^2+4.0*x^2)-(1.0*Dpar*y^5)/(3.0*y^2+3.0*x^2)-(4.0*Dperp*x^2*y^3)/(3.0*y^2+3.0*x^2)+(Dpar*x^2*y^3)/(3.0*y^2+3.0*x^2)+(Dperp*x^2*y^2)/(3.0*y^2+3.0*x^2)-(1.0*Dpar*x^2*y^2)/(3.0*y^2+3.0*x^2)+(3.0*Dpar*x*y^4)/(2.0*y^2+2.0*x^2)-(1.0*Dperp*x^2*y^3)/(2.0*y^2+2.0*x^2)+(Dpar*x^2*y^3)/(2.0*y^2+2.0*x^2)+(10.0*Dperp*x*y^3)/(2.0*y^2+2.0*x^2)+(8.0*Dperp*x^3*y^2)/(2.0*y^2+2.0*x^2)-(1.0*Dpar*x^3*y^2)/(2.0*y^2+2.0*x^2)+(10.0*Dperp*x^3*y)/(2.0*y^2+2.0*x^2)+(Dpar*x^5)/(2.0*y^2+2.0*x^2)-(6.0*Dpar*x*y^5)/(y^2+x^2)-(48.0*Dperp*x^3*y^3)/(y^2+x^2)+(24.0*Dpar*x^3*y^3)/(y^2+x^2)-(4.0*Dperp*x^2*y^3)/(y^2+x^2)+(2.0*Dpar*x^2*y^3)/(y^2+x^2)+(4.0*Dperp*x^3*y^2)/(y^2+x^2)-(3.0*Dpar*x^3*y^2)/(y^2+x^2)-(6.0*Dpar*x^5*y)/(y^2+x^2)-(1.0*Dpar*x^4*y)/(y^2+x^2) 
   end

elseif srcFuncType == "blob" then
   solFunc = function (x,y)
      return 0.0
   end
   
   srcFunc = function (x,y)
      local xc, yc = 0.25, 0.0
      return math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
   end

elseif srcFuncType == "2blobs" then
   solFunc = function (x,y)
      return 0.0
   end
   
   srcFunc = function (x,y)
      local xc, yc = 0.25, 0.25
      local v1 =  math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
      xc, yc = -0.25, 0.25
      local v2 =  math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
      return v1-v2
   end   

elseif srcFuncType == "4xyblobs" then
   solFunc = function (x,y)
      return 0.0
   end
   
   srcFunc = function (x,y)
      local xc, yc = 0.25, 0.25
      local v = math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
      
      xc, yc = -0.25, -0.25
      v = v + math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))

      xc, yc = -0.25, 0.25
      v = v + math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))      

      xc, yc = 0.25, -0.25
      v = v + math.exp(-((x-xc)^2+(y-yc)^2)/(2*0.05^2))
      
      return x*y*v
   end

elseif srcFuncType == "constant" then
   solFunc = function (x,y)
      return 0.0
   end
   
   srcFunc = function (x,y)
      return 1.0
   end

elseif srcFuncType == "pulse" then
   solFunc = function (x,y)
      return 0.0
   end
   
   srcFunc = function (x,y)
      local dx, dy = 1/NX_G, 1/NY_G

      local val = 0.0
      for i = 1, #pList do
	 local px, py = pList[i][1], pList[i][2]
	 local xlo, xup = -0.5+(px-1)*dx, -0.5+px*dx
	 local ylo, yup = -0.5+(py-1)*dy, -0.5+py*dy
	 if (xlo <= x and x <= xup) and (ylo <=y and y <= yup) then
	    val = val + pVals[i]
	 end
      end
      return val
   end
end

local solFn = function(t, z) return solFunc(z[1],z[2]) end
local srcFn = function(t, z) return srcFunc(z[1],z[2]) end



for polyOrder = 1,maxPolyOrder do
   for ncIdx = 1, #numCellListX do
      local grid = Grid.RectCart {
         lower = {-0.5, -0.5},
         upper = {0.5, 0.5},
         cells = {numCellListX[ncIdx], numCellListY[ncIdx]},
         periodicDirs = {}
      }
      -- so pulse can access this
      NX_G = numCellListX[ncIdx]
      NY_G = numCellListY[ncIdx]
      
      local basis = {
         Basis.CartModalSerendipity {
            ndim = grid:ndim(),
            polyOrder = polyOrder
         }
      }
      if polyOrder > 1 then
         basis[2] = Basis.CartModalTensor {
            ndim = grid:ndim(),
            polyOrder = polyOrder
         }
      end
      for bIdx = 1,#basis do
         local function getField()
            return DataStruct.Field {
               onGrid = grid,
               numComponents = basis[bIdx]:numBasis(),
               ghost = {1, 1},
               metaData = {
                  polyOrder = basis[bIdx]:polyOrder(),
                  basisType = basis[bIdx]:id(),
               },
            }
         end
         
         local src = getField()
         local initSource = Updater.ProjectOnBasis {
            onGrid = grid,
            basis = basis[bIdx],
            numQuad = numQuad,
            evaluate = srcFn,
         }
         initSource:advance(0.0, {}, {src})

         local solExact = getField()
         local initSol = Updater.ProjectOnBasis {
            onGrid = grid,
            basis = basis[bIdx],
            numQuad = numQuad,
            evaluate = solFn,
         }
         initSol:advance(0.0, {}, {solExact})

         local Dxx = getField()
         local initDxx = Updater.ProjectOnBasis {
            onGrid = grid,
            basis = basis[bIdx],
            numQuad = numQuad,
            evaluate = DxxFn,
            projectOnGhosts = true,
         }
         initDxx:advance(0.0, {}, {Dxx})

         local Dyy = getField()
         local initDyy = Updater.ProjectOnBasis {
            onGrid = grid,
            basis = basis[bIdx],
            numQuad = numQuad,
            evaluate = DyyFn,
            projectOnGhosts = true,
         }
         initDyy:advance(0.0, {}, {Dyy})
         
         local Dxy = getField()
         local initDxy = Updater.ProjectOnBasis {
            onGrid = grid,
            basis = basis[bIdx],
            numQuad = numQuad,
            evaluate = DxyFn,
            projectOnGhosts = true,
         }
         initDxy:advance(0.0, {}, {Dxy})
         
         local solSim = getField()
	 local tmStart = Time.clock()
         local discontPoisson = Updater.DiscontGenPoisson {
            onGrid = grid,
            basis = basis[bIdx],
            Dxx = Dxx,
            Dyy = Dyy,
            Dxy = Dxy,
            bcLower = bcLower,
            bcUpper = bcUpper,
            writeMatrix = true,
         }
	 local tmEnd = Time.clock()
	 print(string.format("Setup took %g secs", tmEnd-tmStart))

	 tmStart = Time.clock()
         discontPoisson:advance(0.0, {src}, {solSim})
	 tmEnd = Time.clock()
	 print(string.format("Solve took %g secs", tmEnd-tmStart))

         src:write(string.format('%s_nc%d_p%d_src.bp',
                                    basis[bIdx]:id(),
                                    ncIdx,
                                    polyOrder),
                   0.0, 0)
         solExact:write(string.format('%s_nc%d_p%d_solExact.bp',
                                      basis[bIdx]:id(),
                                      ncIdx,
                                      polyOrder),
                        0.0, 0)
         solSim:write(string.format('%s_nc%d_p%d_solSim.bp',
                                    basis[bIdx]:id(),
                                    ncIdx,
                                    polyOrder),
                      0.0, 0)
         Dxx:write(string.format('%s_nc%d_p%d_Dxx.bp',
                                 basis[bIdx]:id(),
                                 ncIdx,
                                 polyOrder),
                   0.0, 0)
         Dyy:write(string.format('%s_nc%d_p%d_Dyy.bp',
                                 basis[bIdx]:id(),
                                 ncIdx,
                                 polyOrder),
                   0.0, 0)
         Dxy:write(string.format('%s_nc%d_p%d_Dxy.bp',
                                 basis[bIdx]:id(),
                                 ncIdx,
                                 polyOrder),
                   0.0, 0)
      end
   end
end
