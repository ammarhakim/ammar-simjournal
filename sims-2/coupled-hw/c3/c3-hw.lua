-- Input file for Hasegawa-Wakatani

----------------------------------
-- Problem dependent parameters --
----------------------------------

polyOrder = 2
coupleCoeff = 0.3 -- adiabacity coefficient
cfl = 0.5/(2*polyOrder+1)

LX, LY = 40, 40 -- domain size
NX, NY = 32, 32 -- number of cells
numOverlappingCells = 4 -- number of overlapping cells

-- (extend left domain to overlap with right domain)
NX_LEFT = NX/2+numOverlappingCells
NX_RIGHT = NX/2
dx = LX/NX -- cell-size

-- sanity
assert(NX_LEFT+NX_RIGHT-numOverlappingCells == NX, "Overlapping grids are not matched")

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------

grid = Grid.RectCart2D {
   lower = {-LX/2, -LY/2},
   upper = {LX/2, LY/2},
   cells = {NX, NY},
   periodicDirs = {0, 1},
}
basis = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid,
   polyOrder = polyOrder,
}

gridLeft = Grid.RectCart2D {
   lower = {-LX/2, -LY/2},
   upper = {0.0+numOverlappingCells*dx, LY/2},
   cells = {NX_LEFT, NY},
   periodicDirs = {1},
}
basisLeft = NodalFiniteElement2D.SerendipityElement {
   onGrid = gridLeft,
   polyOrder = polyOrder,
}

gridRight = Grid.RectCart2D {
   lower = {0.0, -LY/2},
   upper = {LX/2, LY/2},
   cells = {NX_RIGHT, NY},
   periodicDirs = {1},
}
basisRight = NodalFiniteElement2D.SerendipityElement {
   onGrid = gridRight,
   polyOrder = polyOrder,
}

-- vorticity on full domain
chiFull = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- number density on full domain
numDensFull = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- vorticity (left domain)
chiLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
-- extra fields for performing RK updates
chiNewLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
chi1Left = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
chiDupLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}

-- vorticity (right domain)
chiRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
   writeGhost = {0, 0},
}
-- extra fields for performing RK updates
chiNewRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
chi1Right = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
chiDupRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}

-- number density fluctuations (left domain)
numDensLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
-- extra fields for performing RK updates
numDensNewLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
numDens1Left = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
numDensDupLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}

-- number density fluctuations (right domain)
numDensRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
-- extra fields for performing RK updates
numDensNewRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
numDens1Right = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
numDensDupRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}

-- background number density (remains fixed, left domain)
numDensBackLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
   writeGhost = {1, 1},
}
-- total number density
numDensTotalLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}

-- background number density (remains fixed, right domain)
numDensBackRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
   writeGhost = {1, 1},   
}
-- total number density
numDensTotalRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}

-- Poisson source term
poissonSrc = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

--  discontinuous field for potential (for use in source term)
phiDG = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
phiDGDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- potential on left/right grids
phiDGLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}
phiDGDupLeft = DataStruct.Field2D {
   onGrid = gridLeft,
   numComponents = basisLeft:numNodes(),
   ghost = {1, 1},
}

phiDGRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}
phiDGDupRight = DataStruct.Field2D {
   onGrid = gridRight,
   numComponents = basisRight:numNodes(),
   ghost = {1, 1},
}

-----------------------
-- Utility functions --
-----------------------

-- generic function to run an updater
function runUpdater(updater, currTime, timeStep, inpFlds, outFlds)
   updater:setCurrTime(currTime)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(currTime+timeStep)
end

--------------------------------
-- INITIAL CONDITION UPDATERS --
--------------------------------

function oneVortex(x,y,z)
   local r1 = x^2 + y^2
   local s = 2.0
   return math.exp(-r1/s^2)
end

function chiVortex(x,y,z)
   local r1 = x^2 + y^2
   local s = 2.0
   local beta = 0.1
   return 4*(r1-s^2)*math.exp(-r1/s^2)/s^4
end

-- create updater to initialize vorticity
function makeInitChi(grid, basis)
   return Updater.EvalOnNodes2D {
      onGrid = grid,
      basis = basis,
      shareCommonNodes = false, -- In DG, common nodes are not shared
      evaluate = function (x,y,z,t)
	 return -chiVortex(x,y,z)
      end
   }
end

-- initialize vorticity on left/right domains
initChiLeft = makeInitChi(gridLeft, basisLeft)
runUpdater(initChiLeft, 0.0, 0.0, {}, {chiLeft})

initChiRight = makeInitChi(gridRight, basisRight)
runUpdater(initChiRight, 0.0, 0.0, {}, {chiRight})

-- create updater to initialize background density
function makeInitNumDensBack(grid, basis)
   return Updater.EvalOnNodes2D {
      onGrid = grid,
      basis = basis,
      shareCommonNodes = false, -- In DG, common nodes are not shared
      evaluate = function (x,y,z,t)
	 return x
      end
   }
end

-- initialize background density on left/right domains
initNumDensBackLeft = makeInitNumDensBack(gridLeft, basisLeft)
runUpdater(initNumDensBackLeft, 0.0, 0.0, {}, {numDensBackLeft})

initNumDensBackRight = makeInitNumDensBack(gridRight, basisRight)
runUpdater(initNumDensBackRight, 0.0, 0.0, {}, {numDensBackRight})

numDensBackLeft:write("numDensBackLeft.h5", 0.0)
numDensBackRight:write("numDensBackRight.h5", 0.0)

-- create updater to number density
function makeInitNumDens(grid, basis)
   return Updater.EvalOnNodes2D {
      onGrid = grid,
      basis = basis,
      shareCommonNodes = false, -- In DG, common nodes are not shared
      evaluate = function (x,y,z,t)
	 return oneVortex(x,y,z)
      end
   }
end

-- initialize density on left/right domains
initNumDensLeft = makeInitNumDens(gridLeft, basisLeft)
runUpdater(initNumDensLeft, 0.0, 0.0, {}, {numDensLeft})

initNumDensRight = makeInitNumDens(gridRight, basisRight)
runUpdater(initNumDensRight, 0.0, 0.0, {}, {numDensRight})

----------------------
-- EQUATION SOLVERS --
----------------------

-- updater for Poisson bracket
function makePbSlvr(grid, basis)
   return Updater.PoissonBracket {
      onGrid = grid,
      basis = basis,
      cfl = cfl,
      fluxType = "upwind",
      hamilNodesShared = false,
      onlyIncrement = true,
   }
end
pbSlvrLeft = makePbSlvr(gridLeft, basisLeft)
pbSlvrRight = makePbSlvr(gridRight, basisRight)

-- updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson2D {
   onGrid = grid,
   basis = basis,
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- solution is  discontinous   
   periodicDirs = {0, 1}
}

-- set boundaries by copy stuff over to ghost cells
setCouplingBCs = Updater.OverlappingFieldCopy2D {
   onGrid = grid,
   dir = 0, -- direction of overlap
   numOverlappingCells = numOverlappingCells, -- number of overlapping cells
   copyPeriodicDirs = true,
   noiseLevel = 0.1, -- noise due to errors
}
-- average subdomains to compute global field
averageToGlobal = Updater.OverlappingFieldAverage2D {
   onGrid = grid,
   dir = 0, -- direction of overlap
   numOverlappingCells = numOverlappingCells, -- number of overlapping cells
}
-- copy from global to sub-domains
copyToSubdomains = Updater.OverlappingFieldSplit2D {
   onGrid = grid,
   dir = 0, -- direction of overlap
   numOverlappingCells = numOverlappingCells, -- number of overlapping cells
}

-------------------------
-- Boundary Conditions --
-------------------------

-- function to apply copy boundary conditions
function applyBc(fld)
   fld:sync()
end

function applyCouplingBc(fldLeft, fldRight)
   fldLeft:sync(); fldRight:sync()
   -- set interior ghost cells by copying from overlapping regions
   runUpdater(setCouplingBCs, 0.0, 0.0, {}, {fldLeft, fldRight})
end

----------------------
-- SOLVER UTILITIES --
----------------------

-- generic function to run an updater
function runUpdater(updater, currTime, timeStep, inpFlds, outFlds)
   updater:setCurrTime(currTime)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(currTime+timeStep)
end

poissonSolveTime = 0.0 -- time spent in Poisson solve
-- function to solve Poisson equation
function solvePoissonEqn(srcFld, outFld)
   local t1 = os.time() -- begin timer
   -- set poissonSrc <- -srcFld
   poissonSrc:combine(-1.0, srcFld)
   applyBc(poissonSrc)
   local s, dt, msg = runUpdater(poissonSlvr, 0.0, 0.0, {poissonSrc}, {outFld})

   if (s == false) then
      Lucee.logError(string.format("Poisson solver failed to converge (%s)", msg))
   end
   applyBc(outFld)

   poissonSolveTime = poissonSolveTime + os.difftime(os.time(), t1)
end

-- function to solve vorticity equations
function solveVorticityEqn(pbSlvr, t, dt, vortIn, vortOut, phiIn)
   local pbStatus, pbDt = runUpdater(pbSlvr, t, dt, {vortIn, phiIn}, {vortOut})
   vortOut:scale(dt)
   vortOut:accumulate(1.0, vortIn)
   return pbStatus, pbDt
end

-- function to solve number density equations
function solveNumDensEqn(pbSlvr, t, dt, numDensTotal, numDensBack, numDensIn, numDensOut, phiIn)
   -- accumulate background number density
   numDensTotal:combine(1.0, numDensBack, 1.0, numDensIn)
   local pbStatus, pbDt = runUpdater(pbSlvr, t, dt, {numDensTotal, phiIn}, {numDensOut})
   numDensOut:scale(dt)
   numDensOut:accumulate(1.0, numDensIn)
   return pbStatus, pbDt
end

-- solve Poisson equation to determine initial potential
applyCouplingBc(chiLeft, chiRight)
applyCouplingBc(numDensLeft, numDensRight)
-- copy chi into a global field for use in Poisson solver
runUpdater(averageToGlobal, 0.0, 0.0, {chiLeft, chiRight}, {chiFull})
solvePoissonEqn(chiFull, phiDG)
-- now split potential between left/right subdomains
runUpdater(copyToSubdomains, 0.0, 0.0, {phiDG}, {phiDGLeft, phiDGRight})
applyCouplingBc(phiDGLeft, phiDGRight)

-- function to compute diagnostics
function calcDiagnostics(tc, dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local myStatus, myDtSuggested

   ---------------------------------------------- RK

   -- RK stage 1 (left)
   solveNumDensEqn(pbSlvrLeft, tCurr, myDt, numDensTotalLeft, numDensBackLeft, numDensLeft, numDens1Left, phiDGLeft)
   local myStatusL, myDtSuggestedL = solveVorticityEqn(pbSlvrLeft, tCurr, myDt, chiLeft, chi1Left, phiDGLeft)

   -- accumulate source term into updated solutions
   numDens1Left:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDensLeft)
   chi1Left:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDensLeft)

   -- RK stage 1 (right)
   solveNumDensEqn(pbSlvrRight, tCurr, myDt, numDensTotalRight, numDensBackRight, numDensRight, numDens1Right, phiDGRight)
   local myStatusR, myDtSuggestedR = solveVorticityEqn(pbSlvrRight, tCurr, myDt, chiRight, chi1Right, phiDGRight)

   -- accumulate source term into updated solutions
   numDens1Right:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDensRight)
   chi1Right:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDensRight)

   myStatus, myDtSuggested = myStatusL and myStatusR, math.min(myDtSuggestedL, myDtSuggestedR)
   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end   
   
   -- apply BCs
   applyCouplingBc(chi1Left, chi1Right)
   applyCouplingBc(numDens1Left, numDens1Right)
   -- copy chi into a global field for use in Poisson solver
   runUpdater(averageToGlobal, 0.0, 0.0, {chi1Left, chi1Right}, {chiFull})
   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chiFull, phiDG)
   -- now split potential between left/right subdomains
   runUpdater(copyToSubdomains, 0.0, 0.0, {phiDG}, {phiDGLeft, phiDGRight})
   applyCouplingBc(phiDGLeft, phiDGRight)

   ---------------------------------------------- RK   
   
   -- RK stage 2 (left)
   solveNumDensEqn(pbSlvrLeft, tCurr, myDt, numDensTotalLeft, numDensBackLeft, numDens1Left, numDensNewLeft, phiDGLeft)
   local myStatusL, myDtSuggestedL = solveVorticityEqn(pbSlvrLeft, tCurr, myDt, chi1Left, chiNewLeft, phiDGLeft)

   -- accumulate source term into updated solutions
   numDensNewLeft:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDens1Left)
   chiNewLeft:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDens1Left)

   chi1Left:combine(3.0/4.0, chiLeft, 1.0/4.0, chiNewLeft)
   numDens1Left:combine(3.0/4.0, numDensLeft, 1.0/4.0, numDensNewLeft)   

   -- RK stage 2 (right)
   solveNumDensEqn(pbSlvrRight, tCurr, myDt, numDensTotalRight, numDensBackRight, numDens1Right, numDensNewRight, phiDGRight)
   local myStatusR, myDtSuggestedR = solveVorticityEqn(pbSlvrRight, tCurr, myDt, chi1Right, chiNewRight, phiDGRight)
   
   -- accumulate source term into updated solutions
   numDensNewRight:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDens1Right)
   chiNewRight:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDens1Right)

   chi1Right:combine(3.0/4.0, chiRight, 1.0/4.0, chiNewRight)
   numDens1Right:combine(3.0/4.0, numDensRight, 1.0/4.0, numDensNewRight)   
   
   myStatus, myDtSuggested = myStatusL and myStatusR, math.min(myDtSuggestedL, myDtSuggestedR)
   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   -- apply BCs
   applyCouplingBc(chi1Left, chi1Right)
   applyCouplingBc(numDens1Left, numDens1Right)
   -- copy chi into a global field for use in Poisson solver
   runUpdater(averageToGlobal, 0.0, 0.0, {chi1Left, chi1Right}, {chiFull})
   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chiFull, phiDG)
   -- now split potential between left/right subdomains
   runUpdater(copyToSubdomains, 0.0, 0.0, {phiDG}, {phiDGLeft, phiDGRight})
   applyCouplingBc(phiDGLeft, phiDGRight)

   ---------------------------------------------- RK 

   -- RK stage 3 (left)
   solveNumDensEqn(pbSlvrLeft, tCurr, myDt, numDensTotalLeft, numDensBackLeft, numDens1Left, numDensNewLeft, phiDGLeft)
   local myStatusL, myDtSuggestedL = solveVorticityEqn(pbSlvrLeft, tCurr, myDt, chi1Left, chiNewLeft, phiDGLeft)

   -- accumulate source term into updated solutions
   numDensNewLeft:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDens1Left)
   chiNewLeft:accumulate(-coupleCoeff*myDt, phiDGLeft, -coupleCoeff*myDt, numDens1Left)

   chi1Left:combine(1.0/3.0, chiLeft, 2.0/3.0, chiNewLeft)
   numDens1Left:combine(1.0/3.0, numDensLeft, 2.0/3.0, numDensNewLeft)

   -- RK stage 3 (right)
   solveNumDensEqn(pbSlvrRight, tCurr, myDt, numDensTotalRight, numDensBackRight, numDens1Right, numDensNewRight, phiDGRight)
   local myStatus, myDtSuggested = solveVorticityEqn(pbSlvrRight, tCurr, myDt, chi1Right, chiNewRight, phiDGRight)

   -- accumulate source term into updated solutions
   numDensNewRight:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDens1Right)
   chiNewRight:accumulate(-coupleCoeff*myDt, phiDGRight, -coupleCoeff*myDt, numDens1Right)

   chi1Right:combine(1.0/3.0, chiRight, 2.0/3.0, chiNewRight)
   numDens1Right:combine(1.0/3.0, numDensRight, 2.0/3.0, numDensNewRight)      

   myStatus, myDtSuggested = myStatusL and myStatusR, math.min(myDtSuggestedL, myDtSuggestedR)
   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   -- apply BCs
   applyCouplingBc(chi1Left, chi1Right)
   applyCouplingBc(numDens1Left, numDens1Right)
   -- copy chi into a global field for use in Poisson solver
   runUpdater(averageToGlobal, 0.0, 0.0, {chi1Left, chi1Right}, {chiFull})
   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chiFull, phiDG)
   -- now split potential between left/right subdomains
   runUpdater(copyToSubdomains, 0.0, 0.0, {phiDG}, {phiDGLeft, phiDGRight})
   applyCouplingBc(phiDGLeft, phiDGRight)

   -- copy over solution
   chiLeft:copy(chi1Left)
   chiRight:copy(chi1Right)   
   numDensLeft:copy(numDens1Left)
   numDensRight:copy(numDens1Right)

   return myStatus, myDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local lastGood = 0.0

   -- main loop
   while tCurr<=tEnd do
      -- copy stuff in case we need to retake step
      chiDupLeft:copy(chiLeft)
      chiDupRight:copy(chiRight)
      
      numDensDupLeft:copy(numDensLeft)
      numDensDupRight:copy(numDensRight)

      phiDGDupLeft:copy(phiDGLeft)
      phiDGDupRight:copy(phiDGRight)
      
      phiDGDup:copy(phiDG)
      
      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- take a time-step
      local advStatus, advDtSuggested = rk3(tCurr, myDt)
      lastGood = advDtSuggested

      if (advStatus == false) then
	 -- time-step too large
	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, advDtSuggested))
	 -- copy in case current solutions were messed up
	 chiLeft:copy(chiDupLeft)
	 chiRight:copy(chiDupRight)
	 
	 numDensLeft:copy(numDensDupLeft)
	 numDensRight:copy(numDensDupRight)

	 phiDGLeft:copy(phiDGDupLeft)
	 phiDGRight:copy(phiDGDupRight)

	 phiDG:copy(phiDGDup)
	 myDt = advDtSuggested
      else
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)

	 tCurr = tCurr + myDt
	 myDt = advDtSuggested
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end

      myDt = math.min(myDt, 0.01)
   end

   return lastGood
end

-- write out data
function writeFields(frame, tm)
   phiDG:write( string.format("phi_%d.h5", frame), tm)
   
   chiLeft:write( string.format("chiLeft_%d.h5", frame), tm )
   chiRight:write( string.format("chiRight_%d.h5", frame), tm)

   chiFull:write( string.format("chiFull_%d.h5", frame), tm)
   
   numDensLeft:write( string.format("numDensLeft_%d.h5", frame), tm )
   numDensRight:write( string.format("numDensRight_%d.h5", frame), tm )
end

-- write out initial conditions
writeFields(0, 0.0)

-- parameters to control time-stepping
tStart = 0.0
tEnd = 100.0
dtSuggested = 0.01 -- initial time-step to use (will be adjusted)
nFrames = 10
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   local retDtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   dtSuggested = retDtSuggested
   -- write out data
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

Lucee.logInfo (string.format("Poisson solves took %g seconds", poissonSolveTime))


