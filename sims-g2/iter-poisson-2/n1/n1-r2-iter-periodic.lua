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
local lower = {-1.0, -1.0}
local upper = {1.0, 1.0}
local cells = {16, 16}
local periodicDirs = {1, 2}

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

local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = 2*polyOrder+1,
   evaluate = function(t, xn)
      local x, y = xn[1], xn[2]
      return math.exp(-10*(x^2+y^2))
   end,
}

local iterPoisson = Updater.IterPoisson {
   onGrid = grid,
   basis = basis,

   errEps = 1e-8, -- maximum residual error
   stepper = 'richard2',
   verbose = true,
}

initSource:advance(0.0, {}, {fIn})
fIn:write('fIn.bp', 0.0, 0)

local tmStart = Time.clock()
iterPoisson:advance(0.0, {fIn}, {fOut})
print(string.format("Simulation took %g", Time.clock()-tmStart))

fOut:write('fOut.bp', 0.0, 0)

iterPoisson:writeDiagnostics()
