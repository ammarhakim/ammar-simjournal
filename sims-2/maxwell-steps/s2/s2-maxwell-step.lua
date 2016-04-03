-- EM waves reflecting off a sphere

Pi = Lucee.Pi
log = Lucee.logInfo

lightSpeed = 1.0

-- wavelength and frequency of EM wave
wavelength = 0.2
freq = lightSpeed/wavelength

Lx = 20.0
Ly = 20.0
rad = 2.0

-- resolution and time-stepping
NX = 250
NY = 250
cfl = 0.9
tStart = 0.0
tEnd = 10*Lx/lightSpeed 
nFrames = 20

-- compute coordinate of left 
dx = Lx/NX
xFirstEdge = -Lx/2+dx

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 8,
   ghost = {2, 2},
}

function circle(x, y, x0, y0, rad)
   return (x-x0)^2+(y-y0)^2<rad^2
end
-- in/out field representing embedded object
inOut = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}
inOut:set(
   function (x,y,z)
      local xc, yc = 0.0, 0.0
      return  circle(x, y, 0.0, 0.0, rad) and -1.0 or 1.0
   end
)
inOut:sync()
-- write field
inOut:write("inOut.h5")

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   -- no fields initially
   return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end


------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

bcElcFld = BoundaryCondition.ZeroTangent { components = {0, 1, 2} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {3, 4, 5} }
bcPot = BoundaryCondition.Copy { components = {6, 7}, fact = {-1, 1} }
--FIXME: fact in bcPot

-- updater for embedded BC (solid wall)
embeddedBcUpdater = Updater.StairSteppedBc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcElcFld, bcMgnFld, bcPot},
   -- in/out field
   inOutField = inOut,
}

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt, dir)
   for i,bc in ipairs({}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- copy BCs on top, bottom and back
   fld:applyCopyBc(0, "upper")
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")

   -- apply BCs on embedded boundary
   embeddedBcUpdater:setDir(dir)
   embeddedBcUpdater:setOut( {fld} )
   embeddedBcUpdater:advance(tCurr+myDt)   
   
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = 1.0,
   mgnErrorSpeedFactor = 1.0,
}

maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,   
}
maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,   
}

-- Current source from an "antenna"
currentSrc = PointSource.Function {
   -- contributes to E_z component
   outComponents = {2},
   -- source term to apply
   source = function (x,y,z,t)
      local J0 = 1.0
      if (x<xFirstEdge) then
	 return -J0*math.sin(freq*t)
      else
	 return 0.0
      end
   end,
}

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.GridOdePointIntegrator2D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {currentSrc},
}

-- function to update source terms
function updateSource(fld, tCurr, t)
   sourceSlvr:setOut( {fld} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = False
   -- X-direction updates
   for i,slvr in ipairs({maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if (myStatus == false) then
      return myStatus, myDtSuggested, false
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr, 1)

   -- Y-direction updates
   for i,slvr in ipairs({maxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested, false
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr, 0)
   
   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(qNewtCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, 0)   
   
   return status, dtSuggested, useLaxSolver
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
function runSimulation(tStart, tEnd, nFrames, initDt)

   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested
   local useLaxSolver = false

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
      if (useLaxSolver) then
	 -- call Lax solver if positivity violated
	 log (string.format(" Taking step %5d at time %6g with dt %g (using Lax solvers)", step, tCurr, myDt))
	 status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
	 useLaxSolver = false
      else
	 log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
	 status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
	 -- time-step too large
	 log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
	 q:copy(qDup)
      elseif (useLaxSolver == true) then
	 -- negative density/pressure occured
	 log (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
	 q:copy(qDup)
	 qNew:copy(qNewDup)
      else
	 -- check if a nan occured
	 if (qNew:hasNan()) then
	    log (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
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
-- setup initial condition
q:set(init)
q:sync()
qNew:copy(q)

maxSlvrDir0:setIn( {q} )
maxSlvrDir0:setOut( {qX} )

maxSlvrDir1:setIn( {qX} )
maxSlvrDir1:setOut( {qNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)


