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
   local nFrames = 1
   local errEps = tbl.errEps and tbl.errEps or 1.0e-6
   local maxSteps = tbl.maxSteps and tbl.maxSteps or 1000
   
   local Dxx, Dyy = tbl.D.Dxx, tbl.D.Dyy
   local Dxy, Dyx = tbl.D.Dxy, tbl.D.Dxy -- Dxy = Dyx

   local cfl = 0.5*cflFrac/(2*polyOrder+1)

   local updateKernels
   if tbl.useFivePointStencil then
      updateKernels = dofile("../code/bad-diffusion-kernels.lua")
   else
      updateKernels = dofile("../code/diffusion-kernels.lua")
   end

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
   local basis = Basis.CartModalTensor {
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

   local src = getField()
   local exactSol = getField()

   -- function for source (optional)
   local srcFunc = function (t, xn) return 0 end
   if tbl.source  then srcFunc = tbl.source end
   -- function for exact solution (optional)
   local exactFunc = function (t, xn) return 0 end
   if tbl.exact  then exactFunc = tbl.exact end

   --------------
   -- Updaters --
   --------------
   local initDist = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = tbl.init,
   }
   local initSrc = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = srcFunc,
   }
   local initExactSol = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = exactFunc,
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

   -- compute integral of field
   local function integrateField(f1)
      local localRange = f1:localRange()
      local indexer = f1:genIndexer()
      local dx, dy = grid:dx(1), grid:dx(2)

      local intf = 0.0
      for idxs in localRange:colMajorIter() do
	 local f1Itr = f1:get(indexer(idxs))
	 intf = intf + f1Itr[1]/2
      end
      return intf*dx*dy
   end

   --------------------
   -- Initialization --
   --------------------
   initDist:advance(0.0, {}, {f})
   applyBc(f)
   f:write("f_0.bp", 0, 0)
   
   initSrc:advance(0.0, {}, {src})
   src:write("src.bp", 0, 0)

   local vol = (grid:upper(1)-grid:lower(1))*(grid:upper(2)-grid:lower(2))
   local srcInt = integrateField(src)/vol -- mean integrated source
   
   initExactSol:advance(0.0, {}, {exactSol})
   exactSol:write("exactSol.bp", 0, 0)

   local updateKernel = updateKernels[polyOrder]

   local diffCoeff = {}
   diffCoeff.Dxx = Dxx
   diffCoeff.Dxy = Dxy
   diffCoeff.Dyx = Dyx
   diffCoeff.Dyy = Dyy

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

	 local sr = src:get(indexer(idxs))

	 -- compute increment
	 updateKernel(diffCoeff, dxCells, fTL, fT, fTR, fL, f, fR, fBL, fB, fBR, kerOut)
	 local fO = fOut:get(indexer(idxs))

	 -- need to handle averages to ensure proper normalization
	 fO[1] = f[1] + dt*(kerOut[1] + sr[1] - 2*srcInt)
	 for k = 2, fIn:numComponents() do
	    fO[k] = f[k] + dt*(kerOut[k] + sr[k])
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

   -- compute the L2 norm of f1-f2
   local function l2diff(f1, f2)
      local localRange = f1:localRange()
      local indexer = f1:genIndexer()
      local dx, dy = grid:dx(1), grid:dx(2)

      local l2 = 0.0
      for idxs in localRange:colMajorIter() do
	 local f1Itr = f1:get(indexer(idxs))
	 local f2Itr = f2:get(indexer(idxs))

	 for k = 1, f1:numComponents() do
	    l2 = l2 + (f1Itr[k]-f2Itr[k])^2/4.0
	 end
      end
      return math.sqrt(l2*dx*dy)
   end

   -- run simulation with RK3
   return function ()
      local tmStart = Time.clock()
      
      local step = 1
      local dx, dy = grid:dx(1), grid:dx(2)
      local omegaCFL = Dxx/dx^2 + 2*math.abs(Dxy)/(dx*dy) + Dyy/dy^2
      local dt = cfl/omegaCFL
      local isDone = false

      while not isDone do
      	 rk3(dt, f, fNew)
	 
      	 local err = l2diff(f, fNew)
      	 print(string.format("Step %d, error %g", step, err))
      	 if err < errEps or step>maxSteps then
      	    isDone = true
      	 end
      	 f:copy(fNew)
      	 step = step + 1
      end
      if step>maxSteps then
	 print("WARNING: Solution has not converged! Increase 'maxSteps'")
      end
      
      f:write("f_1.bp", 1.0)

      print(string.format("\nSimulation took %g sec", Time.clock()-tmStart))
   end
end

return App
