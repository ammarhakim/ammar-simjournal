local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)
function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
function a(j) return 1-b(j) end
function w1(s) return 4/(s^2+s-2) end
function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end
function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
function nu(j) return -(j-1)/j*b(j)/b(j-2) end
function gbar(s,j) return -a(j-1)*mubar(s,j) end

function calcNumStages(dhdp, extraStages) 
   return math.ceil(math.sqrt(4*dhdp+9/4) - 1/2) + extraStages
end

local App = function(tbl)
   -- read in stuff from input table
   local polyOrder = tbl.polyOrder
   local cflFrac = 1.0
   local nFrames = 1
   local errEps = tbl.errEps and tbl.errEps or 1.0e-6
   local maxSteps = tbl.maxSteps and tbl.maxSteps or 1000
   local fact = tbl.factor and tbl.factor or 100
   local extraStages = tbl.extraStages and tbl.extraStages or 1

   -- initial factor to start iteration
   local initFact = tbl.initFactor and tbl.initFactor or fact
   local initFactNumSteps = tbl.initFactorNumSteps and tbl.initFactorNumSteps or 0
   
   local Dxx, Dyy = 1.0, 1.0
   local Dxy, Dyx = 0.0, 0.0

   local cfl = fact*0.5*cflFrac/(2*polyOrder+1)
   local cflInit = initFact*0.5*cflFrac/(2*polyOrder+1)

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
   local fNew = getField()
   local fDup = getField()
   local fNewDup = getField()
   local fDiff0 = getField()
   local fDiff = getField()
   local fJ = getField()
   local fJ1 = getField()
   local fJ2 = getField()

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

   -- compute the L2 norm
   local function l2norm(f1)
      local localRange = f1:localRange()
      local indexer = f1:genIndexer()
      local dx, dy = grid:dx(1), grid:dx(2)

      local l2 = 0.0
      for idxs in localRange:colMajorIter() do
	 local f1Itr = f1:get(indexer(idxs))
	 for k = 1, f1:numComponents() do
	    l2 = l2 + f1Itr[k]^2/4.0
	 end
      end
      return math.sqrt(l2*dx*dy)
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

   local srcL2 = l2norm(src) -- L2 norm of source

   initExactSol:advance(0.0, {}, {exactSol})
   exactSol:write("exactSol.bp", 0, 0)

   local updateKernel = updateKernels[polyOrder]

   local diffCoeff = {}
   diffCoeff.Dxx = Dxx
   diffCoeff.Dxy = Dxy
   diffCoeff.Dyx = Dyx
   diffCoeff.Dyy = Dyy

   local function calcRHS(fIn, fOut)
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
	 fO[1] = kerOut[1] + sr[1] - 2*srcInt
	 for k = 2, fIn:numComponents() do
	    fO[k] = kerOut[k] + sr[k]
	 end
      end
   end

   print(string.format("Number of stages = %d", calcNumStages(fact, extraStages)))
   
   local function sts(dt, fIn, fOut, fact)
      local numStages = calcNumStages(fact, extraStages)

      -- we need this in each stage
      calcRHS(fIn, fDiff0)

      -- stage 1
      fJ2:copy(fIn)
      fJ1:combine(1.0, fIn, mubar(numStages,1)*dt, fDiff0)
      applyBc(fJ1)

      -- rest of stages
      for j = 2, numStages do
	 calcRHS(fJ1, fDiff)
	 fJ:combine(mu(j), fJ1, nu(j), fJ2, 1-mu(j)-nu(j), fIn,
		    mubar(numStages,j)*dt, fDiff, gbar(numStages,j)*dt, fDiff0)
	 applyBc(fJ)
	 
	 -- reset fields for next stage
	 fJ2:copy(fJ1)
	 fJ1:copy(fJ)
      end
      fOut:copy(fJ)
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

   -- run simulation with STS
   return function ()
      local tmStart = Time.clock()

      local fn = io.open(GKYL_OUT_PREFIX .. "_err.txt", "w")
      
      local step = 1
      local dx, dy = grid:dx(1), grid:dx(2)
      local omegaCFL = Dxx/dx^2 + 2*math.abs(Dxy)/(dx*dy) + Dyy/dy^2
      local isDone = false
      local dt

      while not isDone do
	 if step < initFactNumSteps then
	    dt = cflInit/omegaCFL
	    -- we may want to take a few steps with smaller factor
	    sts(dt, f, fNew, initFact)
	 else
	    dt = cfl/omegaCFL
	    sts(dt, f, fNew, fact)
	 end
	 
      	 local err = l2diff(f, fNew)
	 local resNorm = err/dt/srcL2
      	 print(string.format("Step %d, dt = %g. Error = %g (Res. norm = %g)", step, dt, err, resNorm))
	 fn:write(string.format("%d %g\n", step, err))
      	 if err < errEps or step>=maxSteps then
      	    isDone = true
      	 end
      	 f:copy(fNew)
      	 step = step + 1
      end
      if step>maxSteps then
	 print("WARNING: Solution has not converged! Increase 'maxSteps'")
      end
      
      f:write("f_1.bp", 1.0)
      io.close(fn)

      print(string.format("\nSimulation took %g sec", Time.clock()-tmStart))
   end
end

return App
