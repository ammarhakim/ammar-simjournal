local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local App = function(tbl)
   -- read stuff from input table
   local polyOrder = tbl.polyOrder
   local bcLower = tbl.bcLower
   local bcUpper = tbl.bcUpper
   
   local numQuad = 8
   
   local Dxx, Dyy = tbl.Dxx, tbl.Dyy
   local Dxy = tbl.Dxy

   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
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

   local src = getField()
   local initSource = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      numQuad = numQuad,
      evaluate = tbl.src,
   }
   initSource:advance(0.0, {}, {src})

   local solSim = getField()
   local solExact = getField()

   local solFn = function(t,z) return 0 end
   if tbl.sol then solFn = tbl.sol end
   
   local initSol = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      numQuad = numQuad,
      evaluate = solFn,
   }
   initSol:advance(0.0, {}, {solExact})

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

   local discontPoisson = Updater.DiscontGenPoisson {
      onGrid = grid,
      basis = basis,
      Dxx = Dxx, Dyy = Dyy, Dxy = Dxy,
      bcLower = bcLower,
      bcUpper = bcUpper,
      writeMatrix = true,
   }

   Dxx:write("Dxx.bp")
   Dyy:write("Dyy.bp")
   Dxy:write("Dxy.bp")
   
   return function ()
      local tmStart = Time.clock()
      discontPoisson:advance(0.0, {src}, {solSim})
      local tmEnd = Time.clock()
      print(string.format("Solve took %g secs", tmEnd-tmStart))

      src:write("src.bp")
      solSim:write("solSim.bp")
      solExact:write("solExact.bp")
   end
end

return App
