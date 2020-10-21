-- load app code
local PoissonDiffusion = dofile("../code/sts-diffusion-mod.lua")

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Eq = require "Eq.ConstDiffusion"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local Updater = require "Updater"

local polyOrder = 2
local lower = {-1.0, -1.0}
local upper = {1.0, 1.0}
local cells = {32, 32}

poissonSlvr = PoissonDiffusion {
   polyOrder = polyOrder,
   lower = lower,
   upper = upper,
   cells = cells,
   errEps = 1e-8,
   cflFrac = 0.85,
   factor = 120,
   extraStages = 3,
   stepper = 'RKL1',
   extrapolateInterval = 2,
}

local grid = Grid.RectCart {
   lower = lower,
   upper = upper,
   cells = cells,   
   periodicDirs = {1, 2},
}

-- basis functions
local basis = Basis.CartModalSerendipity {
   ndim = grid:ndim(),
   polyOrder = polyOrder,
}

local function getField(numComponents)
   return DataStruct.Field {
      onGrid = grid,
      numComponents = numComponents,
      ghost = {1, 1},
      syncCorners = true,
      metaData = {
	 polyOrder = basis:polyOrder(),
	 basisType = basis:id(),
      },
   }
end

local srcFunc = function(t, xn)
   local x, y = xn[1], xn[2]
   return math.exp(-10*(x^2+y^2))
end

local initSrc = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = srcFunc,
}

local src = getField(basis:numBasis())
local fOut = getField(basis:numBasis())

-- compute RHS
initSrc:advance(0.0, {}, {src})

poissonSlvr(src, fOut)

src:write("src.bp", 0)
fOut:write("fOut.bp", 0)
