-- Input file for Hasegawa-Wakatani

----------------------------------
-- Problem dependent parameters --
----------------------------------

polyOrder = 2
coupleCoeff = 0.3 -- adiabacity coefficient
cfl = 0.5/(2*polyOrder+1)

LX, LY = 40, 40 -- domain size
NX, NY = 32, 32 -- number of cells

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

-- vorticity
chi = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- extra fields for performing RK updates
chiNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
chi1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
chiDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- number density fluctuations
numDens = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- extra fields for performing RK updates
numDensNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
numDens1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
numDensDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- background number density (remains fixed)
numDensBack = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- total number density
numDensTotal = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
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
initChi = Updater.EvalOnNodes2D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   evaluate = function (x,y,z,t)
      return -chiVortex(x,y,z)
   end
}
initChi:setOut( {chi} )
initChi:advance(0.0)

-- create updater to initialize background density
initNumDensBack = Updater.EvalOnNodes2D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   evaluate = function (x,y,z,t)
      return x
   end
}
initNumDensBack:setOut( {numDensBack} )
initNumDensBack:advance(0.0)

-- create updater to number density
initNumDens = Updater.EvalOnNodes2D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   evaluate = function (x,y,z,t)
      return oneVortex(x,y,z)
   end
}
initNumDens:setOut( {numDens} )
initNumDens:advance(0.0) -- time is irrelevant

----------------------
-- EQUATION SOLVERS --
----------------------

-- updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   basis = basis,
   cfl = cfl,
   fluxType = "upwind",
   hamilNodesShared = false,
   onlyIncrement = true,
}
-- updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson2D {
   onGrid = grid,
   basis = basis,
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- solution is  discontinous   
   periodicDirs = {0, 1}
}

-- updater to initialize chi
copyCToD = Updater.CopyContToDisCont2D {
   onGrid = grid,
   basis = basis,
}

-------------------------
-- Boundary Conditions --
-------------------------

-- function to apply copy boundary conditions
function applyBc(fld)
   fld:sync()
end

-- apply BCs
applyBc(chi)
applyBc(numDens)
applyBc(phiDG)

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
   local s, dt, msg = runUpdater(poissonSlvr, 0.0, 0.0, {poissonSrc}, {outFld})

   if (s == false) then
      Lucee.logError(string.format("Poisson solver failed to converge (%s)", msg))
   end
   applyBc(outFld)

   poissonSolveTime = poissonSolveTime + os.difftime(os.time(), t1)
end

-- function to solve vorticity equations
function solveVorticityEqn(t, dt, vortIn, vortOut, phiIn)
   local pbStatus, pbDt = runUpdater(pbSlvr, t, dt, {vortIn, phiIn}, {vortOut})
   vortOut:scale(dt)
   vortOut:accumulate(1.0, vortIn)
   return pbStatus, pbDt
end

-- function to solve number density equations
function solveNumDensEqn(t, dt, numDensIn, numDensOut, phiIn)
   -- accumulate background number density
   numDensTotal:combine(1.0, numDensBack, 1.0, numDensIn)
   local pbStatus, pbDt = runUpdater(pbSlvr, t, dt, {numDensTotal, phiIn}, {numDensOut})
   numDensOut:scale(dt)
   numDensOut:accumulate(1.0, numDensIn)
   return pbStatus, pbDt
end

-- function to copy potential to DG field
function copyPotential(tCurr, dt, cgIn, dgOut)
   copyCToD:setCurrTime(tCurr)
   copyCToD:setIn( {cgIn} )
   copyCToD:setOut( {dgOut} )
   copyCToD:advance(tCurr+dt)
end

-- solve Poisson equation to determine initial potential
solvePoissonEqn(chi, phiDG)

-- function to compute diagnostics
function calcDiagnostics(tc, dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   solveNumDensEqn(tCurr, myDt, numDens, numDens1, phiDG)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi, chi1, phiDG)

   -- accumulate source term into updated solutions
   numDens1:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens)
   chi1:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   -- apply BCs
   applyBc(chi1)
   applyBc(numDens1)
   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi1, phiDG)

   -- RK stage 2
   solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phiDG)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phiDG)

   -- accumulate source term into updated solutions
   numDensNew:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens1)
   chiNew:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens1)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)
   numDens1:combine(3.0/4.0, numDens, 1.0/4.0, numDensNew)
   -- apply BCs
   applyBc(chi1)
   applyBc(numDens1)
   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi1, phiDG)

   -- RK stage 3
   solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phiDG)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phiDG)

   -- accumulate source term into updated solutions
   numDensNew:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens1)
   chiNew:accumulate(-coupleCoeff*myDt, phiDG, -coupleCoeff*myDt, numDens1)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   numDens1:combine(1.0/3.0, numDens, 2.0/3.0, numDensNew)
   -- apply BCs
   applyBc(chi1)
   applyBc(numDens1)

   -- copy over solution
   chi:copy(chi1)
   numDens:copy(numDens1)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi, phiDG)

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
      -- copy chi over
      chiDup:copy(chi)
      phiDGDup:copy(phiDG)
      numDensDup:copy(numDens)

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
	 chi:copy(chiDup)
	 phiDG:copy(phiDGDup)
	 numDens:copy(numDensDup)
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

   end

   return lastGood
end

-- write out data
function writeFields(frame)
   phiDG:write( string.format("phi_%d.h5", frame) )
   chi:write( string.format("chi_%d.h5", frame) )
   numDens:write( string.format("numDens_%d.h5", frame) )
end

-- write out initial conditions
writeFields(0)

-- parameters to control time-stepping
tStart = 0.0
tEnd = 100.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 10
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   retDtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   dtSuggested = retDtSuggested
   -- write out data
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

Lucee.logInfo (string.format("Poisson solves took %g seconds", poissonSolveTime))


