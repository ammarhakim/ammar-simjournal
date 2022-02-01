-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local polyOrder = 2
local lower = {0, 0}
local upper = {2*math.pi, 2*math.pi}
local cells = {32, 32}
local periodicDirs = {1,2}

local grid = Grid.RectCart {
   lower = {lower[1], lower[2]},
   upper = {upper[1], upper[2]},
   cells = {cells[1], cells[2]},
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity {
   ndim = grid:ndim(),
   polyOrder = polyOrder
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
local fIn = getField()
local fOut = getField()
local fExact = getField()

-- Initial conditions from:
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html
local initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x, y = z[1], z[2]
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
}
local initExact = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x, y = z[1], z[2]
      local a, b = 2, 5
      local c1, d0 = 0, 0
      local c0 = a/12 - 1/2
      local d1 = b/12 - 1/2
      local t1 = x^2/2 - a*x^4/12 + c0*x + c1
      local t2 = y^2/2 - b*y^4/12 + d0*y + d1
      return t1*t2
   end,
}

local discontPoisson = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { { }, { } },
   bcUpper = { { }, { } },
}

local tmStart = Time.clock()

initDist:advance(0.0, {}, {fIn})
discontPoisson:advance(0.0, {fIn}, {fOut})

local tmEnd = Time.clock()

print(string.format("Direct solver took %g sec", tmEnd-tmStart))

fIn:write('src.bp', 0.0, 0)
--fExact:write('fExact.bp', 0.0, 0)
fOut:write('f_1.bp', 0.0, 0)
