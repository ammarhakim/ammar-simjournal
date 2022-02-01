
log = Lucee.logInfo

-- physical parameters
gasGamma = 1.4

Lx = 1.0
Ly = 1.0

-- resolution and time-stepping
NX = 400
NY = 400
cfl = 0.9
tStart = 0.0
tEnd = 0.25
nFrames = 5

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {Lx, Ly},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}

-- aliases to various sub-systems
fluid = q:alias(0, 5)
fluidX = qX:alias(0, 5)
fluidNew = qNew:alias(0, 5)

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   -- Case 4 from Liska and Wendroff 2003 (see Table 4.3)
   -- these tables store pressure, rho, u and v
   local upLeft = {0.35, 0.5065, 0.8939, 0.0}
   local upRight = {1.1, 1.1, 0.0, 0.0}
   local loLeft = {1.1, 1.1, 0.8939, 0.8939}
   local loRight = {0.35, 0.5065, 0.0, 0.8939}

   local rho, u, v, pr
   if (y>0.5) then
      -- upper half
      if (x<0.5) then
	 pr, rho, u, v = upLeft[1], upLeft[2], upLeft[3], upLeft[4]
      else
	 pr, rho, u, v = upRight[1], upRight[2], upRight[3], upRight[4]
      end
   else
      -- lower half
      if (x<0.5) then
	 pr, rho, u, v = loLeft[1], loLeft[2], loLeft[3], loLeft[4]
      else
	 pr, rho, u, v = loRight[1], loRight[2], loRight[3], loRight[4]
      end
   end
   
   return rho, rho*u, rho*v, 0.0, 0.5*rho*(u^2+v^2) + pr/(gasGamma-1)
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
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(0, "upper")
   fld:applyCopyBc(1, "upper")
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
eulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}

-- ds solvers for regular Euler equations along X
fluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- ds solvers for regular Euler equations along Y
fluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- ds solvers for Lax Euler equations along X
fluidLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax Euler equations along Y
fluidLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false
   -- X-direction updates
   for i,slvr in ipairs({fluidSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if (myStatus == false) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   if (eulerEqn:checkInvariantDomain(fluidX) == false) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({fluidSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if (eulerEqn:checkInvariantDomain(fluidNew) == false) then
       useLaxSolver = true
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({fluidLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({fluidLaxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

-- function to take one time-step with Lax Euler solver
function solveTwoFluidLaxSystem(tCurr, t)
   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- dynvector to store fluid energy
fluidEnergy = DataStruct.DynVector { numComponents = 1 }
fluidEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (rho, rhou, rhov, rhow, er)
		  return er
	       end,
}
fluidEnergyCalc:setIn( {fluid} )
fluidEnergyCalc:setOut( {fluidEnergy} )

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({fluidEnergyCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )
   fluidEnergy:write( string.format("fluidEnergy_%d.h5", frame) )
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

-- set input/output arrays for various solvers
fluidSlvrDir0:setIn( {fluid} )
fluidSlvrDir0:setOut( {fluidX} )

fluidSlvrDir1:setIn( {fluidX} )
fluidSlvrDir1:setOut( {fluidNew} )

fluidLaxSlvrDir0:setIn( {fluid} )
fluidLaxSlvrDir0:setOut( {fluidX} )

fluidLaxSlvrDir1:setIn( {fluidX} )
fluidLaxSlvrDir1:setOut( {fluidNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)


