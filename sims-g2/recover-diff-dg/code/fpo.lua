local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

local App = function(tbl)
   -- read in stuff from input table
   local polyOrder = tbl.polyOrder
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local tEnd = tbl.tEnd
   local nFrames = tbl.numFrames
   
   local cells = {8, 8}
   local lower = {0.0, 0.0}
   local upper = {2*math.pi, 2*math.pi}
   local periodicDirs = {1, 2}

   local Dxx, Dyy = tbl.D.Dxx, tbl.D.Dyy
   local Dxy, Dyx = tbl.D.Dxy, tbl.D.Dxy -- Dxy = Dyx

   local cfl = 0.5*cflFrac/(2*polyOrder+1)

   --local updateKernels = dofile("../code/fpo-kernels-c.lua")
   local updateKernels = dofile("../code/fpo-kernels.lua")

   ----------------------
   -- Grids and Fields --
   ----------------------
   local grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = {1, 2},
   }

   -- basis functions
   local basis = Basis.CartModalSerendipity {
      ndim = grid:ndim(),
      polyOrder = polyOrder
   }

   -- fields
   local function getField()
      return DataStruct.Field {
	 onGrid = grid,
	 numComponents = basis:numBasis(),
	 ghost = {1, 1},
	 syncCorners = true,
	 metaData = {
	    polyOrder = basis:polyOrder(),
	    basisType = basis:id(),
	 },
      }
   end
   local f = getField()
   local f1 = getField()
   local f2 = getField()
   local fe = getField()
   local fNew = getField()

   local gpot = getField()

   --------------
   -- Updaters --
   --------------
   local initDist = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = tbl.init,
   }
   local initGpot = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      projectOnGhosts = true,
      evaluate = function (t, xn)
	 local x, y = xn[1], xn[2]
	 return 0.5*(Dxx*x^2 + Dyy*y^2) + Dxy*x*y
      end
   }

   local function applyBc(fld)
      fld:sync()
      -- need to manually sync corners for now
      local globalRange = fld:globalRange()
      local xlo, xup = globalRange:lower(1), globalRange:upper(1)
      local ylo, yup = globalRange:lower(2), globalRange:upper(2)

      local indexer = fld:indexer()
      local idxSkin, idxGhost

      -- lower-left
      idxSkin, idxGhost = fld:get(indexer(xup, yup)), fld:get(indexer(xlo-1, ylo-1))
      for k = 1, fld:numComponents() do
	 idxGhost[k] = idxSkin[k]
      end

      -- lower-right
      idxSkin, idxGhost = fld:get(indexer(xlo, yup)), fld:get(indexer(xup+1, ylo-1))
      for k = 1, fld:numComponents() do
	 idxGhost[k] = idxSkin[k]
      end

      -- upper-left
      idxSkin, idxGhost = fld:get(indexer(xup, ylo)), fld:get(indexer(xlo-1, yup+1))
      for k = 1, fld:numComponents() do
	 idxGhost[k] = idxSkin[k]
      end

      -- upper-right
      idxSkin, idxGhost = fld:get(indexer(xlo, ylo)), fld:get(indexer(xup+1, yup+1))
      for k = 1, fld:numComponents() do
	 idxGhost[k] = idxSkin[k]
      end
   end

   --------------------
   -- Initialization --
   --------------------
   initDist:advance(0.0, {}, {f})
   applyBc(f)
   f:write("f_0.bp", 0, 0)

   -- initialize g
   initGpot:advance(0.0, {}, {gpot})
   gpot:write("g.bp", 0, 0)

   local updateKernel = updateKernels[polyOrder]

   local function forwardEuler(dt, fIn, fOut)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local idxsL, idxsR = {}, {}
      local idxsT, idxsB = {}, {}
      local idxsTL, idxsTR = {}, {}
      local idxsBL, idxsBR = {}, {}   
      local dx, dy = grid:dx(1), grid:dx(2)
      local dxCells = {dx, dy}

      local kerOut = Lin.Vec(fIn:numComponents())

      for idxs in localRange:colMajorIter() do
	 idxsL[1], idxsL[2] = idxs[1]-1, idxs[2]
	 idxsR[1], idxsR[2] = idxs[1]+1, idxs[2]
	 idxsT[1], idxsT[2] = idxs[1], idxs[2]+1
	 idxsB[1], idxsB[2] = idxs[1], idxs[2]-1

	 idxsTL[1], idxsTL[2] = idxs[1]-1, idxs[2]+1
	 idxsTR[1], idxsTR[2] = idxs[1]+1, idxs[2]+1
	 idxsBL[1], idxsBL[2] = idxs[1]-1, idxs[2]-1
	 idxsBR[1], idxsBR[2] = idxs[1]+1, idxs[2]-1

	 local f = fIn:get(indexer(idxs))
	 local fL = fIn:get(indexer(idxsL))
	 local fR = fIn:get(indexer(idxsR))
	 local fT = fIn:get(indexer(idxsT))
	 local fB = fIn:get(indexer(idxsB))

	 local fTL = fIn:get(indexer(idxsTL))
	 local fTR = fIn:get(indexer(idxsTR))
	 local fBL = fIn:get(indexer(idxsBL))
	 local fBR = fIn:get(indexer(idxsBR))

	 local g = gpot:get(indexer(idxs))
	 local gL = gpot:get(indexer(idxsL))
	 local gR = gpot:get(indexer(idxsR))
	 local gT = gpot:get(indexer(idxsT))
	 local gB = gpot:get(indexer(idxsB))

	 local gTL = gpot:get(indexer(idxsTL))
	 local gTR = gpot:get(indexer(idxsTR))
	 local gBL = gpot:get(indexer(idxsBL))
	 local gBR = gpot:get(indexer(idxsBR))

	 -- compute increment
	 updateKernel(
	    dxCells,
	    fTL, fT, fTR, fL, f, fR, fBL, fB, fBR,
	    gTL, gT, gTR, gL, g, gR, gBL, gB, gBR,
	    kerOut)
	 -- update solution
	 local fO = fOut:get(indexer(idxs))
	 for k = 1, fIn:numComponents() do
	    fO[k] = f[k] + dt*kerOut[k]
	 end
      end
   end

   local function rk3(dt, fIn, fOut)
      -- Stage 1
      forwardEuler(dt, fIn, f1)
      applyBc(f1)

      -- Stage 2
      forwardEuler(dt, f1, fe)
      f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
      applyBc(f2)

      -- Stage 3
      forwardEuler(dt, f2, fe)
      fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
      applyBc(fOut)
   end

   -- run simulation with RK3
   return function ()
      local tmStart = Time.clock()
      
      local tCurr = 0.0
      local step = 1
      local dx, dy = grid:dx(1), grid:dx(2)
      local omegaCFL = Dxx/dx^2 + 2*math.abs(Dxy)/(dx*dy) + Dyy/dy^2
      local dt = cfl/omegaCFL

      local frameInt = tEnd/nFrames
      local nextFrame = 1
      local isDone = false

      while not isDone do
	 if (tCurr+dt >= tEnd) then
	    isDone = true
	    dt = tEnd-tCurr
	 end
	 print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))
	 rk3(dt, f, fNew)
	 f:copy(fNew)

	 tCurr = tCurr+dt
	 if tCurr >= nextFrame*frameInt or math.abs(tCurr-nextFrame*frameInt) < 1e-10 then
	    f:write(string.format("f_%d.bp", nextFrame), tCurr, nextFrame)
	    nextFrame = nextFrame+1
	 end
	 step = step+1
      end

      print(string.format("\nSimulation took %g sec", Time.clock()-tmStart))
   end
end

return App
