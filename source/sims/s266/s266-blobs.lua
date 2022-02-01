--[[

In this file we solve the coupled two-field reduced MHD equations
given by

dn/dt + { phi,n } = -alpha*n + D del_per^2 n
dC/dt + { phi,C } = alpha*C - beta*{x,n} + mu del_per^2 C

where n(x,y,t) is number density, C(x,y,t) is vorticity and phi is the
potential determined from

del_perp^2 phi = C

--]]

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.05

-- simulation parameters
alpha = 0.0
beta = 1.0
Ra = 1e4 -- Rayleigh number
D = math.sqrt(1/Ra)
mu = math.sqrt(1/Ra)
n0 = 0.1 -- background density

LX, LY = 60, 40
-- number of cells
NX, NY = 384, 256

-- time to run simulation
tEnd = 25.0
-- number of frames
nFrames = 10

-- grid on which equations are to be solved
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {LX, LY},
   cells = {NX, NY},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2.
   polyOrder = polyOrder,
}

-- number of CG nodes per cell
numCgNodesPerCell = basis:numExclusiveNodes()
-- number of DG nodes per cell
numDgNodesPerCell = basis:numNodes()

-- vorticity
chi = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
-- clear out contents
chi:clear(0.0)

-- extra fields for performing RK updates
chiNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
chi1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
chiDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
chiD2 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- number density fluctuations
numDens = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
-- clear out contents
numDens:clear(0.0)

-- extra fields for performing RK updates
numDensNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
numDens1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
numDensDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
numDensD2 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- Poisson source term
poissonSrc = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- field to store number density advection
numAdvectFld = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
-- clear out contents
numAdvectFld:clear(0.0)

-- updater to solve diffusion equation
diffSolver = Updater.HyperDiffusion2D {
   onGrid = grid,
   basis = basis,
   diffusionCoeff = mu,
   cfl = cfl,
   onlyIncrement = true, -- just compute increments
}

-- create updater to initialize vorticity
initNumDens = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local xc, yc = LX/6.0, LY/2.0
		 local r2 = (x-xc)^2 + (y-yc)^2
		 return n0 + math.exp(-r2/2)
	      end
}
initNumDens:setOut( {numDens} )
-- initialize number density
initNumDens:advance(0.0) -- time is irrelevant

-- A generic function to run an updater.
--
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

-- A HACK
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val}
   end
end
function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   end
end

bcPhi = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}
bcNumDens = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, n0),
}

bcLowerPhi = Updater.Bc2D {
   onGrid = grid,
   boundaryConditions = {bcPhi},
   dir = 0,
   edge = "lower",
}
bcUpperPhi = Updater.Bc2D {
   onGrid = grid,
   boundaryConditions = {bcPhi},
   dir = 0,
   edge = "upper",
}
bcLowerNumDens = Updater.Bc2D {
   onGrid = grid,
   boundaryConditions = {bcNumDens},
   dir = 0,
   edge = "lower",
}
bcUpperNumDens = Updater.Bc2D {
   onGrid = grid,
   boundaryConditions = {bcNumDens},
   dir = 0,
   edge = "upper",
}

-- function to apply boundary conditions
function applyBcPhiChi(fld)
   for i,bc in ipairs({bcLowerPhi, bcUpperPhi}) do
      runUpdater(bc, 0.0, 0.0, {}, {fld})
   end
   fld:applyPeriodicBc(1)
end
function applyBcNumDens(fld)
   for i,bc in ipairs({bcLowerNumDens, bcUpperNumDens}) do
      runUpdater(bc, 0.0, 0.0, {}, {fld})
   end
   fld:applyPeriodicBc(1)
end

-- apply BCs
function applyBc(phi, chi, numDens)
   applyBcPhiChi(phi)
   applyBcPhiChi(chi)
   applyBcNumDens(numDens)
end

-- potential
phi = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex",   
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
phiDup = DataStruct.Field2D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

--  discontinuous field for potential (for use in source term)
phiDG = DataStruct.Field2D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- create updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   -- only compute increments
   onlyIncrement = true,
}

-- define equation to solve
advectionEqn = HyperEquation.Advection {
   -- advection velocity
   speeds = {0.0, 1.0, 0.0}
}
-- updater to solve advection term in vorticity equation
advectSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = advectionEqn,
   -- CFL number
   cfl = cfl,
   -- only compute increments
   onlyIncrement = true,
}

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- flag to indicate if nodes in src field are shared
   sourceNodesShared = false, -- default true
   -- periodic directions
   periodicDirs = {1},
   -- left boundary
   bcLeft = { T = "D", V = 0.0 },
   -- right boundary
   bcRight = { T = "D", V = 0.0 },
}

-- create updater to initialize chi
copyCToD = Updater.CopyContToDisCont2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
}

poissonSolveTime = 0.0 -- time spent in Poisson solve

-- to store center-of-mass
centerOfMass = DataStruct.DynVector { numComponents = 2, }

-- to compute total number of particles in domain
centerOfMassCalc = Updater.CenterOfMass2D {
   onGrid = grid,
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
   -- background value to subtract
   background = n0,
}

-- function to solve Poisson equation
function solvePoissonEqn(srcFld, outFld)
   local t1 = os.time() -- begin timer

   -- set poissonSrc <- -srcFld
   poissonSrc:combine(-1.0, srcFld)

   poissonSlvr:setIn( {poissonSrc} )
   poissonSlvr:setOut( {outFld} )
   -- solve for potential (time is irrelevant here)
   local s, dt, msg = poissonSlvr:advance(0.0)

   if (s == false) then
      Lucee.logError(string.format("Poisson solver failed to converge (%s)", msg))
   end

   poissonSolveTime = poissonSolveTime + os.difftime(os.time(), t1)
end

-- function to solve vorticity equations
function solveVorticityEqn(t, dt, vortIn, vortOut, phiIn, numDensIn)
   pbSlvr:setCurrTime(t)
   pbSlvr:setIn( {vortIn, phiIn} )
   pbSlvr:setOut( {vortOut} )
   -- solve for number density
   local pbStatus, pbDt = pbSlvr:advance(t+dt)
   -- compute number density advection term
   advectSlvr:setCurrTime(t)
   advectSlvr:setIn( {numDensIn} )
   advectSlvr:setOut( {numAdvectFld} )
   local adStatus, adDt = advectSlvr:advance(t+dt)

   -- vortOut <- vortIn + dt*vortOut
   vortOut:scale(dt)
   vortOut:accumulate(1.0, vortIn, beta*dt, numAdvectFld)

   return pbStatus and adStatus, math.min(pbDt, adDt)
end

-- function to solve number density equations
function solveNumDensEqn(t, dt, numDensIn, numDensOut, phiIn)
   pbSlvr:setCurrTime(t)
   pbSlvr:setIn( {numDensIn, phiIn} )
   pbSlvr:setOut( {numDensOut} )
   -- solve for number density
   local pbStatus, pbDt = pbSlvr:advance(t+dt)
   -- numDensOut <- numDensIn + dt*numDensOut
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

-- solve advection equation
function solveDiffusion(curr, dt, qIn, qOut)
   diffSolver:setIn( {qIn} )
   diffSolver:setOut( {qOut} )
   diffSolver:setCurrTime(curr)
   return diffSolver:advance(curr+dt)
end

-- solve Poisson equation to determine initial potential
solvePoissonEqn(chi, phi)
copyPotential(0.0, 0.0, phi, phiDG)

-- function to compute diagnostics
function calcDiagnostics(numDensIn, curr, dt)
   -- for now not doing anything
   centerOfMassCalc:setIn( {numDensIn} )
   centerOfMassCalc:setOut( {centerOfMass} )
   centerOfMassCalc:setCurrTime(curr)
   centerOfMassCalc:advance(curr+dt)
end

-- apply Bcs to initial conditions
applyBc(phi, chi, numDens)

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local myStatus, myDtSuggested
   local myNumDensDiffStatus, myNumDensDiffDtSuggested
   local myChiDiffStatus, myChiDiffDtSuggested
   
   -- RK stage 1
   myStatus, myDtSuggested = solveNumDensEqn(tCurr, myDt, numDens, numDens1, phi)
   myNumDensDiffStatus, myNumDensDiffDtSuggested = solveDiffusion(tCurr, myDt, numDens, numDensD2)

   myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi, chi1, phi, numDens)
   myChiDiffStatus, myChiDtSuggested = solveDiffusion(tCurr, myDt, chi, chiD2)

   numDens1:accumulate(-alpha*myDt, numDens, myDt, numDensD2)
   copyPotential(tCurr, myDt, phi, phiDG)
   chi1:accumulate(alpha*myDt, phiDG, myDt, chiD2)

   if ((myStatus == false) or (myNumDensDiffStatus == false) or (myChiDiffStatus == false)) then
      return false, math.min(myDtSuggested, myNumDensDiffDtSuggested, myChiDtSuggested)
   end
   applyBc(phi, chi1, numDens1)
   solvePoissonEqn(chi1, phi)

   -- RK stage 2
   myStatus, myDtSuggested = solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phi)
   myNumDensDiffStatus, myNumDensDiffDtSuggested = solveDiffusion(tCurr, myDt, numDens1, numDensD2)

   myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phi, numDens1)
   myChiDiffStatus, myChiDtSuggested = solveDiffusion(tCurr, myDt, chi1, chiD2)

   numDensNew:accumulate(-alpha*myDt, numDens1, myDt, numDensD2)
   copyPotential(tCurr, myDt, phi, phiDG)
   chiNew:accumulate(alpha*myDt, phiDG, myDt, chiD2)

   if ((myStatus == false) or (myNumDensDiffStatus == false) or (myChiDiffStatus == false)) then
      return false, math.min(myDtSuggested, myNumDensDiffDtSuggested, myChiDtSuggested)
   end
   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)
   numDens1:combine(3.0/4.0, numDens, 1.0/4.0, numDensNew)
   applyBc(phi, chi1, numDens1)
   solvePoissonEqn(chi1, phi)

   -- RK stage 3
   myStatus, myDtSuggested = solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phi)
   myNumDensDiffStatus, myNumDensDiffDtSuggested = solveDiffusion(tCurr, myDt, numDens1, numDensD2)

   myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phi, numDens1)
   myChiDiffStatus, myChiDtSuggested = solveDiffusion(tCurr, myDt, chi1, chiD2)

   numDensNew:accumulate(-alpha*myDt, numDens1, myDt, numDensD2)
   copyPotential(tCurr, myDt, phi, phiDG)
   chiNew:accumulate(alpha*myDt, phiDG, myDt, chiD2)

   if ((myStatus == false) or (myNumDensDiffStatus == false) or (myChiDiffStatus == false)) then
      return false, math.min(myDtSuggested, myNumDensDiffDtSuggested, myChiDtSuggested)
   end
   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   numDens1:combine(1.0/3.0, numDens, 2.0/3.0, numDensNew)
   applyBc(phi, chi1, numDens1)

   chi:copy(chi1)
   numDens:copy(numDens1)

   solvePoissonEqn(chi, phi)
   copyPotential(tCurr, myDt, phi, phiDG)

   return true, math.min(myDtSuggested, myNumDensDiffDtSuggested, myChiDtSuggested)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local lastGood = 0.0

   -- main loop
   while tCurr<=tEnd do
      chiDup:copy(chi)
      phiDup:copy(phi)
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
	 phi:copy(phiDup)
	 numDens:copy(numDensDup)
	 myDt = advDtSuggested
      else
	 calcDiagnostics(numDens, tCurr, myDt)
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
function writeFields(frame, tm)
   phiDG:write( string.format("phi_%d.h5", frame), tm )
   chi:write( string.format("chi_%d.h5", frame), tm )
   numDens:write( string.format("numDens_%d.h5", frame), tm )
   centerOfMass:write( string.format("centerOfMass_%d.h5", frame), tm )
end

calcDiagnostics(numDens, 0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)

tStart = 0.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
-- parameters to control time-stepping
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   retDtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   dtSuggested = retDtSuggested
   -- write out data
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

Lucee.logInfo (string.format("Poisson solves took %g seconds", poissonSolveTime))


