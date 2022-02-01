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
local lower = {0}
local upper = {1}
local cells = {16}
local periodicDirs = {1}
local nmax = 5 -- terms to retain in FT of step

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
local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
      local x = xn[1]
      local s = 0.0
      for j = 0, nmax do
	 s = s + 4*math.pi^2*(2*j+1)*math.sin(2*math.pi*(2*j+1)*x)
      end
      return -s
   end,
}

local exactSol = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
      local x = xn[1]
      local s = 0.0
      for j = 0, nmax do
	 s = s - math.sin(2*math.pi*(2*j+1)*x)/(2*j+1)
      end
      return s
   end,
}

local iterPoisson = Updater.IterPoisson {
   onGrid = grid,
   basis = basis,

   -- there parameters will eventually be replaced by internal
   -- heuristics
   
   errEps = 1e-8, -- maximum residual error
   factor = 8, -- factor over explicit scheme
   extraStages = 2, -- extra stages
   cflFrac = 0.6, -- CFL frac for internal iterations
   stepper = 'RKL1',
   extrapolateInterval = 2,
   
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
