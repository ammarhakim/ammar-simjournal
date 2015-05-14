-- DG solver for Maxwell equations

log = Lucee.logInfo

-- domain size
X = 80.0 -- [m]
Y = 40.0 -- [m]

-- resolution and time-stepping
NX = 80
NY = 40
polyOrder = 1 -- DG polynomial order
cfl = 0.5/(2*polyOrder+1)
tStart = 0.0
tEnd = 150e-9
nFrames = 2

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {X, Y},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {0, 1},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- number of nodes per cell for DG
numDgNodesPerCell = basis:numNodes()

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8*numDgNodesPerCell,
   ghost = {2, 2},
}

-- for RK time-stepping
q1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8*numDgNodesPerCell,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8*numDgNodesPerCell,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8*numDgNodesPerCell,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8*numDgNodesPerCell,
   ghost = {2, 2},
}

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial condition to apply
function init(x,y,z)
   local m, n = 8, 4
   local a = m*Lucee.Pi/X
   local b = n*Lucee.Pi/Y
   local Ez = 1.0*math.sin(a*x)*math.sin(b*y)
   return 0.0, 0.0, Ez, 0.0, 0.0, 0.0, 0.0, 0.0
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
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- updater to apply initial conditions
initField = Updater.EvalOnNodes2D {
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
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = Lucee.SpeedOfLight,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 1.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 1.0,
}

-- updater to solve hyperbolic equations
maxwellSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = maxwellEqn,
   -- CFL number
   cfl = cfl,
}

-- write initial conditions
q:write("q_0.h5")

-- solve maxwell equation
function solveMaxwell(curr, dt, qIn, qOut)
   maxwellSlvr:setCurrTime(curr)
   maxwellSlvr:setIn( {qIn} )
   maxwellSlvr:setOut( {qOut} )
   return maxwellSlvr:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = solveMaxwell(tCurr, myDt, q, q1)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   applyBc(q1, tCurr, myDt)
   -- RK stage 2
   local myStatus, myDtSuggested = solveMaxwell(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1, tCurr, myDt)
   -- RK stage 3
   local myStatus, myDtSuggested = solveMaxwell(tCurr, myDt, q1, qNew)
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
emEnergy = DataStruct.DynVector { numComponents = 1 }

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
   emEnergy:write( string.format("emEnergy_%d.h5", frame) )
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
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)
	 -- copy updated solution back
	 q:copy(qNew)
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

initField:setOut( {q} )
-- initialize
initField:advance(0.0) -- time is irrelevant
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 1.0
runSimulation(tStart, tEnd, nFrames, initDt)
