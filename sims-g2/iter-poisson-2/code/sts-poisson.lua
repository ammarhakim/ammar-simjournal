local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Eq = require "Eq.ConstDiffusion"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local Updater = require "Updater"

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)
local function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
local function a(j) return 1-b(j) end
local function w1(s) return 4/(s^2+s-2) end
local function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end
local function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
local function nu(j) return -(j-1)/j*b(j)/b(j-2) end
local function gbar(s,j) return -a(j-1)*mubar(s,j) end

local function calcNumStagesRKL2(dhdp, extraStages) 
   return math.ceil(math.sqrt(4*dhdp+9/4) - 1/2) + extraStages
end

-- For RKL1 scheme
local function muRKL1(j) return (2*j-1)/j end
local function nuRKL1(j) return (1-j)/j end
local function mubarRKL1(s,j) return (2*j-1)/j*2/(s^2+s) end

local function calcNumStagesRKL1(dhdp, extraStages) 
   return math.ceil(1/2*(math.sqrt(1+8*dhdp)-1)) + extraStages
end

local App = function(tbl)
   -- read in stuff from input table
   local polyOrder = tbl.polyOrder
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local nFrames = 1
   local errEps = tbl.errEps and tbl.errEps or 1.0e-6
   local maxSteps = tbl.maxSteps and tbl.maxSteps or 1000
   local fact = tbl.factor and tbl.factor or 100
   local extraStages = tbl.extraStages and tbl.extraStages or 1
   
   local extrapolateInterval = tbl.extrapolateInterval and tbl.extrapolateInterval or maxSteps+1

   -- one of 'RKL2' or 'RKL1'
   local stepper = tbl.stepper and tbl.stepper or 'RKL2'

   -- initial factor to start iteration
   local initFact = tbl.initFactor and tbl.initFactor or fact
   local initFactNumSteps = tbl.initFactorNumSteps and tbl.initFactorNumSteps or 0

   local hasExact = tbl.exact and true or false -- flag to indicate if exact sol is specified
   
   local Dxx, Dyy = 1.0, 1.0

   local cfl = fact*0.5*cflFrac/(2*polyOrder+1)
   local cflInit = initFact*0.5*cflFrac/(2*polyOrder+1)

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
   local f = getField(basis:numBasis())
   local fNew = getField(basis:numBasis())
   local fDiff0 = getField(basis:numBasis())
   local fDiff = getField(basis:numBasis())
   local fJ = getField(basis:numBasis())
   local fJ1 = getField(basis:numBasis())
   local fJ2 = getField(basis:numBasis())

   -- for extrapolation
   local fE1 = getField(basis:numBasis())
   local fE2 = getField(basis:numBasis())

   -- source
   local src = getField(basis:numBasis())

   -- exact solution
   local fExact = getField(basis:numBasis())

   -- for time-stepping (never used)
   local cflRateByCell = getField(1)

   -- error history
   local errHist = DataStruct.DynVector { numComponents = 1 }
   -- extrapolation factors
   local extraHist = DataStruct.DynVector { numComponents = 1 }

   -- function for source (optional)
   local srcFunc = function (t, xn) return 0 end
   if tbl.source  then srcFunc = tbl.source end
   -- function for exact solution (optional)
   local exactFunc = function (t, xn) return 0 end
   if tbl.exact then exactFunc = tbl.exact end

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

   -- constant diffusion equation object.
   local constDiffusionCalc = Eq {
      coefficient = {Dxx, Dyy},
      basis = basis,
   }

   -- ipdater to solve the (parabolic) diffusion equation.
   local constDiffusionSlvr = Updater.HyperDisCont {
      onGrid = grid,
      basis = basis,
      cfl = 0.5*cflFrac/(2*polyOrder+1),
      equation = constDiffusionCalc,
   }

   local function applyBc(fld)
      fld:sync()
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

   if hasExact then
      initExactSol:advance(0.0, {}, {fExact})
      fExact:write("fExact.bp", 0.0)
   end

   local vol = (grid:upper(1)-grid:lower(1))*(grid:upper(2)-grid:lower(2))
   local srcInt = integrateField(src)/vol -- mean integrated source

   -- we need to adjust sources when all directions are periodic
   do
      local localRange = src:localRange()
      local indexer = src:genIndexer()
      
      for idxs in localRange:colMajorIter() do
	 local sr = src:get(indexer(idxs))
	 sr[1] = sr[1] - 2*srcInt
      end
   end

   local srcIntZero = integrateField(src)
   local srcL2 = l2norm(src) -- L2 norm of source

   local function calcRHS(fIn, fOut)
      local dt, tCurr = 0.1, 1.0 -- these are ignored
      constDiffusionSlvr:setDtAndCflRate(dt, cflRateByCell)
      constDiffusionSlvr:advance(tCurr, {fIn}, {fOut})
      fOut:accumulate(1.0, src)
   end

   if stepper == 'RKL2' then
      print(string.format("Number of stages = %d", calcNumStagesRKL2(fact, extraStages)))
   else
      print(string.format("Number of stages = %d", calcNumStagesRKL1(fact, extraStages)))
   end
   
   local function stsRKL2(dt, fIn, fOut, fact)
      local numStages = calcNumStagesRKL2(fact, extraStages)
      local nsFh = io.open(GKYL_OUT_PREFIX .. "_numStages", "w")
      nsFh:write(string.format("%d", numStages))
      nsFh:close()

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

   local function stsRKL1(dt, fIn, fOut, fact)
      local mu, nu, mubar = muRKL1, nuRKL1, mubarRKL1
      
      local numStages = calcNumStagesRKL1(fact, extraStages)
      local nsFh = io.open(GKYL_OUT_PREFIX .. "_numStages", "w")
      nsFh:write(string.format("%d", numStages))      
      nsFh:close()

      -- we need this in each stage
      calcRHS(fIn, fDiff0)

      -- stage 1
      fJ2:copy(fIn)
      fJ1:combine(1.0, fIn, mubar(numStages,1)*dt, fDiff0)
      applyBc(fJ1)

      -- rest of stages
      for j = 2, numStages do
	 calcRHS(fJ1, fDiff)
	 fJ:combine(mu(j), fJ1, nu(j), fJ2, mubar(numStages,j)*dt, fDiff)
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

      local step = 1
      local dx, dy = grid:dx(1), grid:dx(2)
      local omegaCFL = Dxx/dx^2 + Dyy/dy^2
      local isDone = false

      local numPrevStored = 0 -- flag to check if we are ready to extrapolate
      local errE1, errE2 = 1e10, 1e10 -- errors for use in extrapolation

      local numStages = calcNumStagesRKL2(fact, extraStages)
      if stepper == 'RKL1' then
	 numStages = calcNumStagesRKL1(fact, extraStages)
      end

      while not isDone do
	 local dt = cfl/omegaCFL
	 if stepper == 'RKL2' then
	    stsRKL2(dt, f, fNew, fact)
	 else
	    stsRKL1(dt, f, fNew, fact)
	 end
	 
      	 local err = l2diff(f, fNew)
	 local resNorm = err/dt/srcL2
      	 print(string.format("Step %d, dt = %g. Error = %g (Res. norm = %g)", step, dt, err, resNorm))

      	 if err < errEps or step>=maxSteps then
      	    isDone = true
      	 end
      	 f:copy(fNew)

	 -- check if we should store the solution for use in
	 -- extrapolation
	 if step % extrapolateInterval == 0 then
	    fE1:copy(fE2)
	    fE2:copy(f)
	    errE1 = errE2
	    errE2 = err
	    numPrevStored = numPrevStored + 1

	    if numPrevStored > 1 then -- need two values to extrapolate
	       local eps = errE2/errE1
	       extraHist:appendData(numPrevStored-1, { eps } )
	       f:combine(1.0, fE2, eps, fE2, -eps, fE1)
	    end
	 end

	 errHist:appendData(numStages*step, { err })	 
	 
      	 step = step + 1
      end
      if step>maxSteps then
	 print("WARNING: Solution has not converged! Increase 'maxSteps'")
      end

      extraHist:write("extraHist.bp", 1.0)
      errHist:write("errHist.bp", 1.0)
      f:write("f_1.bp", 1.0)

      if hasExact then
	 -- compute L2 norm of error and write to file
	 local l2Err = l2diff(fExact, f)
	 local l2Fh = io.open(GKYL_OUT_PREFIX .. "_l2Error", "w")
	 l2Fh:write(string.format("%g", l2Err))
	 l2Fh:close() 
      end
      
      print(
	 string.format(
	    "\nSimulation took %g sec, %d stages", Time.clock()-tmStart, (step-1)*numStages
	 )
      )
   end
end

return App
