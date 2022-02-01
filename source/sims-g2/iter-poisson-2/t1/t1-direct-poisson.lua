-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

local polyOrder = 3
local lower = {0, 0}
local upper = {2*math.pi, 2*math.pi}
local cells = {4, 4}
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
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html#convergence-of-2d-solver
local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
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
}

local exactSol = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
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

local discontPoisson = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { { }, { } },
   bcUpper = { { }, { } },
}

local iterPoisson = Updater.IterPoisson {
   onGrid = grid,
   basis = basis,
}

initSource:advance(0.0, {}, {fIn})
discontPoisson:advance(0.0, {fIn}, {fOut})

exactSol:advance(0.0, {}, {fExact})

local gridVol = 1.0
for d = 1, grid:ndim() do
   gridVol = gridVol*(grid:upper(d)-grid:lower(d))
end

local fOutInt = iterPoisson:integrateField(fOut)/gridVol
local fExactInt = iterPoisson:integrateField(fExact)/gridVol

do
   local localRange = fExact:localRange()
   local indexer = fExact:genIndexer()
   local dfact = math.sqrt(2^grid:ndim())
   
   for idxs in localRange:colMajorIter() do
      local sr = fOut:get(indexer(idxs))
      sr[1] = sr[1] - dfact*fOutInt
   end
   for idxs in localRange:colMajorIter() do
      local sr = fExact:get(indexer(idxs))
      sr[1] = sr[1] - dfact*fExactInt
   end
end

fIn:write('fIn.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)
fExact:write("fExact.bp")

print(string.format("L2 error %g\n", iterPoisson:l2diff(fOut, fExact)))
