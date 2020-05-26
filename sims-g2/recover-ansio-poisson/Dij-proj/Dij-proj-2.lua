local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local x0, y0 = 0.0, 0.0 -- location of O
local zet = 1e1 -- Dpar/Dperp
local Dpar = 1.0 -- parallel heat conduction
local Dperp = Dpar/zet -- perpendicular heat conduction

local bx = function(x, y)
   return -(y-y0)/math.sqrt((x-x0)^2+(y-y0)^2)
end
local by = function(x, y)
   return (x-x0)/math.sqrt((x-x0)^2+(y-y0)^2)
end

tbl = {
   lower = {-0.5, -0.5},
   upper = {0.5, 0.5},
   cells = {16, 16},

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
}

local grid = Grid.RectCart {
   lower = tbl.lower,
   upper = tbl.upper,
   cells = tbl.cells,
}

local basis = Basis.CartModalSerendipity {
   ndim = grid:ndim(),
   polyOrder = 1,
}
   
local function getField()
   return DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {1, 1},
      metaData = {
	 polyOrder = basis:polyOrder(),
	 basisType = basis:id(),
      },
   }
end

local Dxx = getField()
local initDxx = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = tbl.Dxx,
   projectOnGhosts = true,
}
initDxx:advance(0.0, {}, {Dxx})
   
local Dyy = getField()
local initDyy = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = tbl.Dyy,
   projectOnGhosts = true,
}
initDyy:advance(0.0, {}, {Dyy})
         
local Dxy = getField()
local initDxy = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = tbl.Dxy,
   projectOnGhosts = true,
}
initDxy:advance(0.0, {}, {Dxy})

Dxx:write("Dxx.bp")
Dyy:write("Dyy.bp")
Dxy:write("Dxy.bp")
