-- DG solver for Maxwell equations

log = Lucee.logInfo

gasGamma = 5.0/3.0 -- gas adiabatic index
-- resolution and time-stepping
NX, NY = 50, 50
polyOrder = 1 -- DG polynomial order
cfl = 0.45/(2*polyOrder+1)
tStart = 0.0
tEnd = 2.0
nFrames = 1

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {1, 1},
   cells = {NX, NY},
   decomposition = decomp,
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

-----------------------
-- INITIAL CONDITION --
-----------------------
function initXX(x,y,z)
   -- Noh problem
   local rho = 1.0
   local pr = 1.0
   local u = 0.0
   local v = 0.0
   
   return rho, rho*u, rho*v, 0.0, 0.5*rho*(u^2+v^2) + pr/(gasGamma-1)
end

-- initial condition to apply
function init(x,y,z)
   -- Noh problem
   local rho = 1.0
   local pr = 1e-6
   local vrad = -1.0
   local cost = x/math.sqrt(x^2+y^2)
   local sint = y/math.sqrt(x^2+y^2)
   local u = vrad*cost
   local v = vrad*sint
   
   return rho, rho*u, rho*v, 0.0, 0.5*rho*(u^2+v^2) + pr/(gasGamma-1)
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- solid walls at symmetry lines
bcFluidCopy = BoundaryCondition.NodalDgCopy2D { 
   components = {0, 4},
   basis = basis,
}
bcFluidWall = BoundaryCondition.NodalDgZeroNormal2D { 
   components = {1, 2, 3},
   basis = basis,
}
-- exact solution in outflow region
bcNohFunc = BoundaryCondition.NodalDgFunction2D {
   components = {0, 1, 2, 3, 4},
   basis = basis,
   bc = function(x,y,z,t)
	   local rho = 1.0 + t/math.sqrt(x^2+y^2)
	   local pr = 1e-6
	   local vrad = -1.0
	   local cost = x/math.sqrt(x^2+y^2)
	   local sint = y/math.sqrt(x^2+y^2)
	   local u = vrad*cost
	   local v = vrad*sint

	   return rho, rho*u, rho*v, 0.0, 0.5*rho*(u^2+v^2) + pr/(gasGamma-1)
	end,
}

-- create boundary condition object to apply wall BCs
function createWallBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {
	 bcFluidCopy, bcFluidWall,
      },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- create boundary condition object to apply exact BCs on right/top edges
function createFlowBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = { bcNohFunc },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- walls on left and bottom
bcLeft = createWallBc(0, "lower")
bcBottom = createWallBc(1, "lower")
-- time-dependent exact solution on top and right
bcRight = createFlowBc(0, "upper")
bcTop = createFlowBc(1, "upper")

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   local bcList = {bcLeft, bcRight, bcTop, bcBottom}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   --fld:applyCopyBc(0, "upper")
   --fld:applyCopyBc(1, "upper")
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
      applyBc(qSafe, curr, dt)
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
