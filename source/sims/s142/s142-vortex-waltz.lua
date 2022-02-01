-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.05

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
end

-- Determine number of global nodes per cell for use in creating DG
-- fields.
if (polyOrder == 1) then
   numDgNodesPerCell = 4
elseif (polyOrder == 2) then
   numDgNodesPerCell = 8
end

L = 10.0
Nx = 64
Ny = 64

-- grid on which equations are to be solved
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {L, L},
   cells = {Nx, Ny},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 4 and 8.
   polyOrder = polyOrder,
}

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

-- Poisson source term
poissonSrc = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

function vortices(x,y,z)
   local x1, y1 = 3.5, 5.0
   local x2, y2 = 6.5, 5.0
   local r1 = (x-x1)^2 + (y-y1)^2
   local r2 = (x-x2)^2 + (y-y2)^2
   return math.exp(-r1/0.8^2) + math.exp(-r2/0.8^2)
end

-- create updater to initialize vorticity
initChi = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return vortices(x,y,z)
	      end
}
initChi:setOut( {chi} )
-- initialize potential
initChi:advance(0.0) -- time is irrelevant

-- function to apply copy boundary conditions field
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
end

-- apply copy BCs to vorticity
applyBc(chi)

-- write initial value of vorticity
chi:write("chi_0.h5")

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
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "central",
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

-- solve Poisson equation to determine initial potential
solvePoissonEqn(chi, phi)

-- write out initial potential
phi:write("phi_0.h5")

-- solve Poisson equation to determine Potential
solvePoissonEqn(chi, phi)

-- total energy diagnostic
totalEnergy = DataStruct.DynVector {
   -- number of components in diagnostic
   numComponents = 1,
}

-- updater to compute total energy
energyCalc = Updater.EnergyFromStreamFunction {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
}
-- set input/output (this never changes, so it once)
energyCalc:setIn( {phi} )
energyCalc:setOut( {totalEnergy} )

-- compute initial energy in system
energyCalc:advance(0)

-- total enstrophy diagnostic
totalEnstrophy = DataStruct.DynVector {
   -- number of components in diagnostic
   numComponents = 1,
}

-- updater to compute total energy
enstrophyCalc = Updater.TotalEnstrophy {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
}
-- set input/output (this never changes, so it once)
enstrophyCalc:setIn( {chi} )
enstrophyCalc:setOut( {totalEnstrophy} )

-- compute initial enstrophy of system
enstrophyCalc:advance(0)

-- function to compute diagnostics
function calcDiagnostics(tc, dt)
   energyCalc:setCurrTime(tc)
   energyCalc:advance(tc+dt)

   enstrophyCalc:setCurrTime(tc)
   enstrophyCalc:advance(tc+dt)

end

-- function to take a time-step using RK2 time-stepping scheme
function rk2(tCurr, myDt)
   -- RK stage 1 (chi1 <- chi + L(chi))
   pbSlvr:setCurrTime(tCurr)
   pbSlvr:setIn( {chi, phi} )
   pbSlvr:setOut( {chi1} )
   local myStatus, myDtSuggested = pbSlvr:advance(tCurr+myDt)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   -- apply BCs
   applyBc(chi1)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi1, phi)

   -- RK stage 2 (chiNew <- chi1 + L(chi1))
   pbSlvr:setCurrTime(tCurr)
   pbSlvr:setIn( {chi1, phi} )
   pbSlvr:setOut( {chiNew} )
   local myStatus, myDtSuggested = pbSlvr:advance(tCurr+myDt)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   -- do final update of chi1 <-- 0.5*(chi + chiNew)
   chi1:combine(0.5, chi, 0.5, chiNew)
   -- copy over solution
   chi:copy(chi1)

   -- apply BCs
   applyBc(chi)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi, phi)

   return myStatus, myDtSuggested
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1 (chi1 <- chi + L(chi))
   pbSlvr:setCurrTime(tCurr)
   pbSlvr:setIn( {chi, phi} )
   pbSlvr:setOut( {chi1} )
   local myStatus, myDtSuggested = pbSlvr:advance(tCurr+myDt)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   -- apply BCs
   applyBc(chi1)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi1, phi)

   -- RK stage 2
   pbSlvr:setCurrTime(tCurr)
   pbSlvr:setIn( {chi1, phi} )
   pbSlvr:setOut( {chiNew} )
   local myStatus, myDtSuggested = pbSlvr:advance(tCurr+myDt)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)

   -- apply BCs
   applyBc(chi1)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi1, phi)

   -- RK stage 3
   pbSlvr:setCurrTime(tCurr)
   pbSlvr:setIn( {chi1, phi} )
   pbSlvr:setOut( {chiNew} )
   local myStatus, myDtSuggested = pbSlvr:advance(tCurr+myDt)

   -- check if step failed and return immediately if it did
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   -- apply BCs
   applyBc(chi1)

   -- copy over solution
   chi:copy(chi1)

   -- solve Poisson equation to determine Potential
   solvePoissonEqn(chi, phi)

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
      phiDup:copy(phi)

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
   phi:write( string.format("phi_%d.h5", frame) )
   chi:write( string.format("chi_%d.h5", frame) )
   totalEnergy:write( string.format("totalEnergy_%d.h5", frame) )
   totalEnstrophy:write( string.format("totalEnstrophy_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

Lucee.logInfo (string.format("Poisson solves took %g seconds", poissonSolveTime))

