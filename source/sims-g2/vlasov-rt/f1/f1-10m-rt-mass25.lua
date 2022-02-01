-- Program to solve Two-Fluid equations
Pi = Lucee.Pi
log = Lucee.logInfo

-- decomposition object to use
decomp = DecompRegionCalc2D.CartGeneral {}

-- global parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25.0
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)
mgnErrorSpeedFactor = 1.0

-- problem parameters
Atwood = 0.17
n1 = 1.0 -- heavy-side number density
betaIon = 0.071
gHat = 0.09 -- gHat = g/Omega_i v_A
theta = 0.0 -- inclination angle of field to Z-direction
Te_Ti = 0.1 -- Te/Ti
B0 = 0.1 -- Left wall magnetic field (this may not be correct)
maxPert = 0.01 -- Perturbation amplitude

-- secondary quantities
n2 = n1*(1-Atwood)/(1+Atwood) -- light-side number density
T_ion = (betaIon/n1)*B0^2/(2*mu0) -- ion temperature (uniform)
T_elc = Te_Ti*T_ion -- electron temperature (uniform)

wpe = math.sqrt(n1*elcCharge^2/(epsilon0*elcMass))
wpi = math.sqrt(n1*ionCharge^2/(epsilon0*ionMass))
de = lightSpeed/wpe
di = lightSpeed/wpi
Valf = B0/math.sqrt(mu0*n1*ionMass)
OmegaCe0 = elcCharge*B0/elcMass
OmegaCi0 = ionCharge*B0/ionMass

-- average wave number for pressure damping
elcCollAverageWaveNumber = 1/de
ionCollAverageWaveNumber = 1/de

grav = gHat*OmegaCi0*Valf -- gravitational acceleration

-- domain size
Lx = 3.0*di
Ly = 3.75*di

-- resolution and time-stepping
NX = 128
NY = 160
cfl = 0.5
tStart = 0.0
tEnd = 60/OmegaCi0
nFrames = 60

log(string.format("Atwood=%g", Atwood))
log(string.format("n1=%g, n2=%g", n1, n2))
log(string.format("di=%g", di))
log(string.format("wpe/OmegaCe=%g", math.abs(wpe/OmegaCe0)))
log(string.format("OmegaCi=%g", OmegaCi0))
log(string.format("Lx=%g, Ly=%g", Lx, Ly))
log(string.format("betaIon=%g", betaIon))
log(string.format("Vthe/c=%g", math.sqrt(2*T_elc/elcMass)/lightSpeed))
log(string.format("Valf/c=%g", Valf/lightSpeed))
log(string.format("gravitational acceleration =%g", grav))
log(string.format("tEnd=%g,  nFrames=%d",tEnd,nFrames))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------

-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {Lx, Ly},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {1},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}

-- aliases to various sub-systems
elcFluid = q:alias(0, 10)
ionFluid = q:alias(10, 20)
emField = q:alias(20, 28)

elcFluidX = qX:alias(0, 10)
ionFluidX = qX:alias(10, 20)
emFieldX = qX:alias(20, 28)

elcFluidNew = qNew:alias(0, 10)
ionFluidNew = qNew:alias(10, 20)
emFieldNew = qNew:alias(20, 28)

-----------------------
-- INITIAL CONDITION --
-----------------------

-- set random seed based on processor number to avoid imprinting
-- domain decomp on ICs.
r = Lucee.getRank()
math.randomseed(r+os.time())
math.random()


-- function to apply initial conditions
function initElc(x,y,z)
   local xloc = 0.5*Lx
   local numDens, prElc, Bz
   local B1 = B0
   local maxPert = maxPert -- add this pertubation
   local ph = 0.0
   local xregion = Lx/50.0

   if (x<xloc) then
      -- heavy side
      numDens = n1
      prElc = n1*T_elc
   else
      -- light side
      numDens = n2
      prElc = n2*T_elc
   end

   local rhoe = elcMass*numDens
   local vx = maxPert*Valf*(math.random()-0.5)

   return rhoe, rhoe*vx, 0.0, 0.0, prElc+rhoe*vx*vx, 0.0, 0.0, prElc, 0.0, prElc
end

function initIon(x,y,z)
   local xloc = 0.5*Lx
   local numDens, prIon, Bz
   local B1 = B0
   local maxPert = maxPert -- add this pertubation
   local ph = 0.0
   local xregion = Lx/50.0

   if (x<xloc) then
      -- heavy side
      numDens = n1
      prIon = n1*T_ion
   else
      -- light side
      numDens = n2
      prIon = n2*T_ion
   end

   local rhoi = ionMass*numDens
   local vx = maxPert*Valf*(math.random()-0.5)

   return rhoi, rhoi*vx, 0.0, 0.0, prIon+rhoi*vx*vx, 0.0, 0.0, prIon, 0.0, prIon
end

function initEmField(x,y,z)
   local xloc = 0.5*Lx
   local T = (T_elc + T_ion)
   local B1 = B0
   local maxPert = maxPert -- add this pertubation
   
   if (x<xloc) then
      -- heavy side
      numDens = n1
      Bz = math.sqrt(B1^2 + 2*mu0*(n1*T - numDens*T + ionMass*grav*n1*x))
   else
      -- light side
      numDens = n2
      Bz = math.sqrt(B1^2 + 2*mu0*(n1*T - numDens*T + ionMass*grav*n1*xloc + ionMass*grav*n2*(x-xloc)))
   end

   return 0.0, 0.0, 0.0, 0.0, 0.0, Bz, 0.0, 0.0
end

-- define various equations to solve
elcEqn = HyperEquation.TenMoment {
}
ionEqn = HyperEquation.TenMoment {
}
elcLaxEqn = HyperEquation.TenMoment {
   numericalFlux = "lax",
}
ionLaxEqn = HyperEquation.TenMoment {
   numericalFlux = "lax",
}

maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = lightSpeed,
   -- correction speeds
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-- ds solvers for regular equations along X
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for regular equations along Y
elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- ds solvers for Lax equations along X
elcLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax equations along Y
elcLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- updater for source terms
sourceSlvr = Updater.ImplicitTenMomentSrc2D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
   gravity = grav, -- gravitational acceleration
   dir = 0, -- direction of acceleration
}

-- Collisional source updaters
elcCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = elcCollAverageWaveNumber,
}
ionCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = ionCollAverageWaveNumber,
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, tEnd)
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(tEnd)

   -- electron collisional relaxation
   elcCollSrcSlvr:setOut( {elcIn} )
   elcCollSrcSlvr:setCurrTime(tCurr)
   elcCollSrcSlvr:advance(tEnd)

   -- ion collisional relaxation
   ionCollSrcSlvr:setOut( {ionIn} )
   ionCollSrcSlvr:setCurrTime(tCurr)
   ionCollSrcSlvr:advance(tEnd)   
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields
bcElcCopy = BoundaryCondition.Copy { components = {0} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcElcPrCopyY = BoundaryCondition.Copy { components = {4, 7, 8, 9} }
bcElcPrFlipY = BoundaryCondition.Copy { components = {5, 6}, fact = {-1, -1} }

bcIonCopy = BoundaryCondition.Copy { components = {10} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {11, 12, 13} }
bcIonPrCopyY = BoundaryCondition.Copy { components = {14, 17, 18, 19} }
bcIonPrFlipY = BoundaryCondition.Copy { components = {15, 16}, fact = {-1, -1} }

bcElcFld = BoundaryCondition.ZeroTangent { components = {20, 21, 22} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {23, 24, 25} }
bcPot = BoundaryCondition.Copy { components = {26, 27} }

-- create boundary condition object
function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {
	 bcElcCopy, bcElcWall, 
	 bcElcPrCopyY, bcElcPrFlipY,
	 bcIonCopy, bcIonWall,
	 bcIonPrCopyY, bcIonPrFlipY,
	 bcElcFld, bcMgnFld, bcPot,
      },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- create updaters to apply boundary conditions
bcLeft = createBc(0, "lower")
bcRight = createBc(0, "upper")

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   for i,bc in ipairs({bcLeft, bcRight}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- sync ghost cells
   fld:sync()
end

-- apply BCs to initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)


-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEqn:checkInvariantDomain(elcFluidX) == false)
    or (qX:hasNan()) or (ionEqn:checkInvariantDomain(ionFluidX) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir1, ionFluidSlvrDir1, maxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEqn:checkInvariantDomain(elcFluidNew) == false)
    or (qNew:hasNan()) or (ionEqn:checkInvariantDomain(ionFluidNew) == false)) then
      useLaxSolver = true
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with solver
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

   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elcLaxSlvrDir1, ionLaxSlvrDir1, maxLaxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

-- function to take one time-step with Lax solver
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

   if Lucee.IsRestarting then
      fileName = "q_" .. Lucee.RestartFrame .. ".h5"
      tCurr = q:read(fileName)
      if not tCurr then
         tCurr = tStart + (tEnd - tStart) * Lucee.RestartFrame / numFrames
      end
      frame = Lucee.RestartFrame + 1
      nextIOt = tCurr + tFrame
      log('\nRestarting from frame %d tCurr = %g\n', Lucee.RestartFrame, tCurr)
   else
      elcFluid:set(initElc)
      ionFluid:set(initIon)
      emField:set(initEmField)
   end
   applyBc(q, 0.0, 0.0)
   qNew:copy(q)
   qNew:sync()
   
   if not Lucee.IsRestarting then
      calcDiagnostics(0.0, 0.0)
      writeFields(0, 0.0)
   end

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
elcFluid:set(initElc)
ionFluid:set(initIon)
emField:set(initEmField)
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- set input/output arrays for various solvers
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidX} )
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidX} )
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldX} )

elcFluidSlvrDir1:setIn( {elcFluidX} )
elcFluidSlvrDir1:setOut( {elcFluidNew} )
ionFluidSlvrDir1:setIn( {ionFluidX} )
ionFluidSlvrDir1:setOut( {ionFluidNew} )
maxSlvrDir1:setIn( {emFieldX} )
maxSlvrDir1:setOut( {emFieldNew} )

elcLaxSlvrDir0:setIn( {elcFluid} )
elcLaxSlvrDir0:setOut( {elcFluidX} )
ionLaxSlvrDir0:setIn( {ionFluid} )
ionLaxSlvrDir0:setOut( {ionFluidX} )
maxLaxSlvrDir0:setIn( {emField} )
maxLaxSlvrDir0:setOut( {emFieldX} )

elcLaxSlvrDir1:setIn( {elcFluidX} )
elcLaxSlvrDir1:setOut( {elcFluidNew} )
ionLaxSlvrDir1:setIn( {ionFluidX} )
ionLaxSlvrDir1:setOut( {ionFluidNew} )
maxLaxSlvrDir1:setIn( {emFieldX} )
maxLaxSlvrDir1:setOut( {emFieldNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)
