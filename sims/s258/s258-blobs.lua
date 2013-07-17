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
cfl = 0.5/(2*polyOrder-1)

-- simulation parameters
alpha = 0.0
beta = 1.0
D = 0.0
mu = 0.0

LX, LY = 60, 40
-- number of cells
NX, NY = 192, 128

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

function oneVortex(x,y,z)
   local r1 = x^2 + y^2
   local s = 2.0
   return math.exp(-r1/s^2)
end

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
		 return 0.1 + math.exp(-r2/2)
	      end
}
initNumDens:setOut( {numDens} )
-- initialize number density
initNumDens:advance(0.0) -- time is irrelevant

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
end

-- apply BCs
applyBc(numDens)

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
   periodicDirs = {0, 1}
}

-- create updater to initialize chi
copyCToD = Updater.CopyContToDisCont2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
}

poissonSolveTime = 0.0 -- time spent in Poisson solve

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

-- solve Poisson equation to determine initial potential
solvePoissonEqn(chi, phi)
copyPotential(0.0, 0.0, phi, phiDG)

-- function to compute diagnostics
function calcDiagnostics(tc, dt)
   -- for now not doing anything
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   solveNumDensEqn(tCurr, myDt, numDens, numDens1, phi)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi, chi1, phi, numDens)

   copyPotential(tCurr, myDt, phi, phiDG)
   numDens1:accumulate(-alpha*myDt, numDens)
   chi1:accumulate(alpha*myDt, phiDG)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   applyBc(chi1)
   applyBc(numDens1)
   solvePoissonEqn(chi1, phi)

   -- RK stage 2
   solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phi)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phi, numDens1)

   copyPotential(tCurr, myDt, phi, phiDG)
   numDensNew:accumulate(-alpha*myDt, numDens1)
   chiNew:accumulate(alpha*myDt, phiDG)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)
   numDens1:combine(3.0/4.0, numDens, 1.0/4.0, numDensNew)
   applyBc(chi1)
   applyBc(numDens1)
   solvePoissonEqn(chi1, phi)

   -- RK stage 3
   solveNumDensEqn(tCurr, myDt, numDens1, numDensNew, phi)
   local myStatus, myDtSuggested = solveVorticityEqn(tCurr, myDt, chi1, chiNew, phi,numDens1)

   copyPotential(tCurr, myDt, phi, phiDG)
   numDensNew:accumulate(-alpha*myDt, numDens1)
   chiNew:accumulate(alpha*myDt, phiDG)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   numDens1:combine(1.0/3.0, numDens, 2.0/3.0, numDensNew)
   applyBc(chi1)
   applyBc(numDens1)

   chi:copy(chi1)
   numDens:copy(numDens1)

   solvePoissonEqn(chi, phi)
   copyPotential(tCurr, myDt, phi, phiDG)

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
      chiDup:copy(chi)
      phiDup:copy(phi)
      numDensDup:copy(numDens)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- take a time-step
      local advStatus, advDtSuggested = rk3(tCurr, myDt)
      lastGood = advDtSuggested

      if (advStatus == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, advDtSuggested))
	 -- copy in case current solutions were messed up
	 chi:copy(chiDup)
	 phi:copy(phiDup)
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
tEnd = 40.0
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


