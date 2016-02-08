-- Plasma-beach problem

-- Electrostatic shock problem

log = Lucee.logInfo

-- problem parameters
Te_Ti = 9.0 -- ratio of electron to ion temperaute
machNum = 0.25 -- Mach number computed from ion thermal speed
n0 = 1.0 -- initial number density

-- physical parameters
gasGamma = 5.0/3.0
elcTemp = 1.0e-4
elcMass = 1.0
elcCharge = -1.0

ionTemp = elcTemp/Te_Ti
ionMass = 1836.2
ionCharge = 1.0

-- permittivity of free space
epsilon0 = 1.0

-- thermal speeds
cs = math.sqrt(elcTemp/ionMass)
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- electron and ion drift speeds
elcDrift = machNum*cs
ionDrift = elcDrift -- no net current

-- domain size and simulation time
LX = 500*lambdaD
tEnd = 1000.0/wpe
nFrames = 100

-- Resolution, time-stepping etc.
NX = 1000
cfl = 0.9

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Sound speed=%g", cs))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Electron/Ion drift speed=%g", elcDrift))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc1D.CartGeneral {}
-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {LX},
   cells = {NX},
   decomposition = decomp,
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- aliases to various sub-systems
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   local sloc = 0.5*LX
   local elcU = elcDrift
   local ionU = ionDrift
   if (x>sloc) then
      elcU = -elcDrift
      ionU = -ionDrift
   end
   local rhoElc = elcMass*n0
   local rhoIon = ionMass*n0
   local erElc = n0*elcTemp/(gasGamma-1) + 0.5*rhoElc*elcU^2
   local erIon = n0*ionTemp/(gasGamma-1) + 0.5*rhoIon*ionU^2
   return rhoElc, rhoElc*elcU, 0.0, 0.0, erElc, rhoIon, rhoIon*ionU, 0.0, 0.0, erIon, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
elcEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
ionEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
elcEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
ionEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = 1.0,
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = 0.0,
}

-- ds solvers for regular Euler equations along X
elcFluidSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax Euler equations along X
elcLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- updater for source terms
sourceSlvr = Updater.ImplicitFiveMomentSrc1D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = 1.0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, t)
   -- two-fluid sources
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = False
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({elcLaxSlvrDir0, ionLaxSlvrDir0, maxLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

-- function to take one time-step with Lax Euler solver
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- compute diagnostic
function calcDiagnostics(tCurr, t)
   for i,diag in ipairs({}) do
      diag:setCurrTime(tCurr)
      diag:advance(t)
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

-- set input/output arrays for various solvers
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidNew} )
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidNew} )
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldNew} )

elcLaxSlvrDir0:setIn( {elcFluid} )
elcLaxSlvrDir0:setOut( {elcFluidNew} )
ionLaxSlvrDir0:setIn( {ionFluid} )
ionLaxSlvrDir0:setOut( {ionFluidNew} )
maxLaxSlvrDir0:setIn( {emField} )
maxLaxSlvrDir0:setOut( {emFieldNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

tStart = 0.0
initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)

-- print some timing information
tFluids = elcFluidSlvrDir0:totalAdvanceTime() + ionFluidSlvrDir0:totalAdvanceTime()
tMax = maxSlvrDir0:totalAdvanceTime()
tSource = sourceSlvr:totalAdvanceTime()

Lucee.logInfo(string.format("Total time in fluid solvers = %g", tFluids))
Lucee.logInfo(string.format("Total time in Maxwell solver = %g", tMax))
Lucee.logInfo(string.format("Total time in source solver = %g", tSource))



