-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local Updater = require "Updater"

local polyOrder = 1
local lower = {0, 0, 0}
local upper = {2*math.pi, 2*math.pi, 2*math.pi}
local cells = {64, 64, 64}
local periodicDirs = {1, 2, 3}

local log = Logger { logToFile = false }
local decomp = DecompRegionCalc.CartProd { cuts = {2, 1, 1} }

local grid = Grid.RectCart {
   lower = lower,
   upper = upper,
   cells = cells,
   periodicDirs = periodicDirs,
   decomposition = decomp,
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
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
      local x, y, z = xn[1], xn[2], xn[3]
      local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local t1, t2 = 0.0, 0.0
      local f = 0.0
      for m = 0,2 do
	 for n = 0,2 do
	    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)*math.sin(3*z)
	    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)*math.sin(3*z)
	    f = f + -(m*m+n*n+9)*(t1+t2)
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
      local x, y, z = xn[1], xn[2], xn[3]
      local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local t1, t2 = 0.0, 0.0
      local f = 0.0

      for m = 0,2 do
	 for n = 0,2 do
	    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)*math.sin(3*z)
	    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)*math.sin(3*z)
	    f = f + (t1+t2)
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
   cflFrac = 8.0, -- CFL frac for internal iterations
   stepper = 'richard2',
   verbose = true,
}

initSource:advance(0.0, {}, {fIn})

local tmStart = Time.clock()
iterPoisson:advance(0.0, {fIn}, {fOut})
log(string.format("Simulation took %g\n", Time.clock()-tmStart))

fIn:write('fIn.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)

exactSol:advance(0.0, {}, {fExact})
fExact:write("fExact.bp")

iterPoisson:writeDiagnostics()
log(string.format("L2 error %g\n\n", iterPoisson:l2diff(fOut, fExact)))
