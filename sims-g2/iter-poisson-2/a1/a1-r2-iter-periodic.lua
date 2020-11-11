-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local polyOrder = 1
local lower = {0}
local upper = {2*math.pi}
local cells = {8}
local periodicDirs = {1}

local grid = Grid.RectCart {
   lower = lower,
   upper = upper,
   cells = cells,
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

-- Initial conditions from
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html#convergence-of-2d-solver
local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, xn)
      local x = xn[1]
      local am = {0, 5, -10} 
      local bm = {10, 5, 10}
      local t1, t2 = 0.0, 0.0
      local f = 0.0
      for m = 0,2 do
	 for n = 0,2 do
	    t1 = am[m+1]*math.cos(m*x)
	    t2 = bm[m+1]*math.sin(m*x)
	    f = f-m*m*(t1+t2)
	 end
      end
      return -f/50.0
   end,
}

local exactSol = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, xn)
      local x = xn[1]
      local t1, t2 = 0.0, 0.0
      local am = {0, 5, -10} 
      local bm = {10, 5, 10}
      local f = 0.0

      for m = 0,2 do
	 for n = 0,2 do
	    t1 = am[m+1]*math.cos(m*x)
	    t2 = bm[m+1]*math.sin(m*x)
	    f = f+(t1+t2)
	 end
      end
      return f/50.0
   end,
}

local iterPoisson = Updater.IterPoisson {
   onGrid = grid,
   basis = basis,

   -- there parameters will eventually be replaced by internal
   -- heuristics
   
   errEps = 1e-8, -- maximum residual error
   cflFrac = 1.5, -- CFL frac for internal iterations
   stepper = 'richard2',
   verbose = true,
}

initSource:advance(0.0, {}, {fIn})

local tmStart = Time.clock()
iterPoisson:advance(0.0, {fIn}, {fOut})
print(string.format("Simulation took %g", Time.clock()-tmStart))

fIn:write('fIn.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)

exactSol:advance(0.0, {}, {fExact})
fExact:write("fExact.bp")

iterPoisson:writeDiagnostics()
print(string.format("L2 error %g\n", iterPoisson:l2diff(fOut, fExact)))
