-- DG solver for Maxwell equations

log = Lucee.logInfo

gasGamma = 1.4 -- gas adiabatic index
-- resolution and time-stepping
NX = 400
polyOrder = 1 -- DG polynomial order
cfl = 0.9/(2*polyOrder+1)
tStart = 0.0
tEnd = 0.0039
nFrames = 1

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {0.5},
   cells = {NX},
   decomposition = decomp,
}

-- create FEM nodal basis
basis = NodalFiniteElement1D.LagrangeTensor {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- number of nodes per cell for DG
numDgNodesPerCell = basis:numNodes()

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}

-- for RK time-stepping
q1 = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- extra copy for use in "safe stage" 
qSafe = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial condition to apply
function init(x,y,z)
   local sloc = 0.4
   -- right state
   local rho, u, pr = 6.591493, 2.2654207, 3.1544874
   if (x<sloc) then
      -- left state
      rho, u, pr = 0.1261192, 8.9047029, 782.92899
   end
   return rho, rho*u, 0.0, 0.0, pr/(gasGamma-1) + 0.5*rho*u^2
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   local bcList = {}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------

-- function run specified updater
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

-- updater to apply initial conditions
initField = Updater.ProjectOnNodalBasis1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return init(x,y,z)
	      end
}

-- define equation to solve
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma
}

-- updater to solve hyperbolic equations
eulerSlvr = Updater.NodalDgHyper1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = eulerEqn,
   -- CFL number
   cfl = cfl,
}

-- updaters to fix positivity
eulerFlatten = Updater.NodalPositiveFilter1D {
   onGrid = grid,
   basis = basis,
   equation = eulerEqn,
   operation = "flatten"
}
eulerFilter = Updater.NodalPositiveFilter1D {
   onGrid = grid,
   basis = basis,
   equation = eulerEqn,
   operation = "filter"
}

-- write initial conditions
q:write("q_0.h5")

-- solve euler equation
function solveEuler(curr, dt, qIn, qOut)
   qSafe:copy(qIn) -- make a copy so we don't mess up qIn

   -- solve Euler equation
   local status, dtSuggested = runUpdater(eulerSlvr, curr, dt, {qSafe}, {qOut})
   -- check if positivity was violated
   local isPositive, dtInf = runUpdater(eulerFlatten, curr, dt, {qOut}, {qSafe})
   if (isPositive == false) then
      log ("** Positivity violated, retaking stage ...")
      status, dtSuggested = runUpdater(eulerSlvr, curr, dt, {qSafe}, {qOut})
   end
   -- filter solution
   runUpdater(eulerFilter, curr, dt, {qOut}, {qOut})

   return status, dtSuggested
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = solveEuler(tCurr, myDt, q, q1)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   applyBc(q1, tCurr, myDt)
   -- RK stage 2
   local myStatus, myDtSuggested = solveEuler(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1, tCurr, myDt)
   -- RK stage 3
   local myStatus, myDtSuggested = solveEuler(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   q1:combine(1.0/3.0, q, 2.0/3.0, qNew)
   applyBc(q1, tCurr, myDt)
   qNew:copy(q1)
   return myStatus, myDtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- dynvector to store EM energy

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )
end

----------------------------
-- TIME-STEPPING FUNCTION --
----------------------------
-- function to advance solution from tStart to tEnd
function runSimulation(tStart, tEnd, nFrames, initDt)
   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- the grand loop 
   while true do
      -- copy q and qNew in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
        myDt = tEnd-tCurr
      end

      -- advance fluids and fields
      log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
	 q:copy(qDup)
      else
	 -- check if a nan occured
	 if (qNew:hasNan()) then
	    log (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
	    -- write out data jut before quiting
	    q:write( string.format("q_nan.h5", frame) )
	    qNew:write( string.format("qNew_nan.h5", frame) )
	    break
	 end
	 -- copy updated solution back
	 q:copy(qNew)
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)
	 -- write out data
	 if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then
	    log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
	    writeFields(frame, tCurr+myDt)
	    frame = frame + 1
	    nextIOt = nextIOt + tFrame
	    step = 0
	 end

	 tCurr = tCurr + myDt
	 myDt = dtSuggested
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end 
   end -- end of time-step loop
   return dtSuggested
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- initialize
initField:setOut( {q} )
initField:advance(0.0) -- time is irrelevant

-- initial conditions
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 1.0
runSimulation(tStart, tEnd, nFrames, initDt)
