-- KH test. See
-- http://www.astro.virginia.edu/VITA/ATHENA/kh.html

log = Lucee.logInfo

-- physical parameters
gasGamma = 1.4
gravity = 0.0

-- domain size
XL = -0.5
XU = 0.5
YL = -0.5
YU = 0.5
rhoL = 0.01 -- gradient length scale for density
velL = 0.01 -- gradient length scale for velocity

Lx = XU-XL
Ly = YU-YL

-- resolution and time-stepping
polyOrder = 1
NX = 500
NY = 500
cfl = 0.45/(2*polyOrder+1)
tStart = 0.0
tEnd = 20.0
nFrames = 200

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {XL, YL},
   upper = {XU, YU},
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
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}

-- for RK time-stepping
q1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- extra copy for use in "safe stage" 
qSafe = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}
-- to store solution increment
qIncr = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5*numDgNodesPerCell,
   ghost = {2, 2},
}

-----------------------
-- INITIAL CONDITION --
-----------------------

-- Construct a smooth transition using tanh profiles
-- pr: Right state. pl: Left state. L: Gradient length scale. x0: Location
--
function tanhx(pr, pl, L, x0, X)
   local tanh = math.tanh
   return 0.5*pr*(1+tanh(2*(X-x0)/L)) + 0.5*pl*(1-tanh(2*(X-x0)/L))
end

-- initial conditions
function init(x,y,z)
   local nmodes = 10
   local rho1 = 1.0
   local vx1 = 0.5

   local rho2 = 2.0
   local vx2 = -0.5

   local pi = Lucee.Pi
   local rho, vx = 0.0, 0.0

   if (y<0) then
      rho = tanhx(rho2, rho1, rhoL, -0.25, y)
      vx = tanhx(vx2, vx1, velL, -0.25, y)
   else
      rho = tanhx(rho1, rho2, rhoL, 0.25, y)
      vx = tanhx(vx1, vx2, velL, 0.25, y)
   end

   vx = vx + 0.01*math.sin(nmodes*Lucee.Pi*x/Lx)
   local vy = 0.01*math.sin(nmodes*Lucee.Pi*x/Lx)
   local pr = 2.5

   return rho, rho*vx, rho*vy, 0.0, pr/(gasGamma-1) + 0.5*rho*(vx^2+vy^2)
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
initField = Updater.ProjectOnNodalBasis2D {
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
eulerSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = eulerEqn,
   -- CFL number
   cfl = cfl,
}

-- updaters to fix positivity
eulerFlatten = Updater.NodalPositiveFilter2D {
   onGrid = grid,
   basis = basis,
   equation = eulerEqn,
   operation = "flatten"
}
eulerFilter = Updater.NodalPositiveFilter2D {
   onGrid = grid,
   basis = basis,
   equation = eulerEqn,
   operation = "filter"
}

function solveEuler(curr, dt, qIn, qOut)
   -- solve Euler equation, computing increments   
   local myStatus, myDtSuggested = runUpdater(eulerSlvr, curr, dt, {qIn}, {qOut})
   return myStatus, myDtSuggested
end

-- solve euler equation
function safeSolveEuler(curr, dt, qIn, qOut)
   qSafe:copy(qIn) -- make a copy so we don't mess up qIn

   -- solve Euler equation
   local status, dtSuggested = solveEuler(curr, dt, qSafe, qOut)
   -- check if positivity was violated
   local isPositive, dtInf = runUpdater(eulerFlatten, curr, dt, {qOut}, {qSafe})
   if (isPositive == false) then
      log ("** Positivity violated, retaking stage ...")
      applyBc(qSafe, curr, dt)
      status, dtSuggested = solveEuler(curr, dt, qSafe, qOut)
   end
   -- filter solution
   runUpdater(eulerFilter, curr, dt, {qOut}, {qOut})

   return status, dtSuggested
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = safeSolveEuler(tCurr, myDt, q, q1)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   applyBc(q1, tCurr, myDt)
   -- RK stage 2
   local myStatus, myDtSuggested = safeSolveEuler(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1, tCurr, myDt)
   -- RK stage 3
   local myStatus, myDtSuggested = safeSolveEuler(tCurr, myDt, q1, qNew)
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
