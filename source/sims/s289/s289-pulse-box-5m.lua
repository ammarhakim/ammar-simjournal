-- Program to solve Two-Fluid equations

-- decomposition object to use
decomp = DecompRegionCalc2D.CartGeneral {}

-- global parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
lightSpeed = 1.0
epsilon0 = 1.0
mgnErrorSpeedFactor = 1.0

Lx = 2.0
Ly = 2.0
cfl = 0.49
bGuideFactor = 0.0

nSpecies = 2

NX = 101
NY = 101

-- computational domain
grid = Grid.RectCart2D {
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {0, 1},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- create duplicate copy in case we need to take step again
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- create aliases to various sub-system
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-- function to apply initial conditions
function init(x,y,z)
   local me = elcMass
   local mi = ionMass
   local gasGamma1 = gasGamma-1
   local n0 = 1.0
   local rhoe = n0*elcMass
   local rhoi = n0*ionMass
   local pre, pri = 1.0e-3, 1.0e-3
   local E0 = 1.0e-1
   local B0 = 1.0e-1
   local Ez = E0*math.exp(-25*(x^2+y^2))
   local Bz = B0*math.exp(-25*(x^2+y^2))

   return rhoe, 0.0, 0.0, 0.0, pre/gasGamma1, rhoi, 0.0, 0.0, 0.0, pri/gasGamma1, 0.0, 0.0, Ez, 0.0, 0.0, Bz, 0.0, 0.0
end
qTwoFluids = q:alias(0, 5*nSpecies+8)
-- set initial conditions for fields and fluids
qTwoFluids:set(init)


-- get ghost cells correct
q:sync()
-- copy initial conditions over
qNew:copy(q)

-- define various equations to solve
elcEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
ionEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
elcEulerLaxEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
   -- use Lax fluxes
   numericalFlux = "lax",   
}
ionEulerLaxEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
   -- use Lax fluxes
   numericalFlux = "lax",
}
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = lightSpeed,
   -- correction speeds
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-- updater for electron equations
elcEulerSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
elcEulerSlvr:setIn( {elcFluid} )
elcEulerSlvr:setOut( {elcFluidNew} )

-- updater for ion equations
ionEulerSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
ionEulerSlvr:setIn( {ionFluid} )
ionEulerSlvr:setOut( {ionFluidNew} )

-- updater for Maxwell equations
maxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
maxSlvr:setIn( {emField} )
maxSlvr:setOut( {emFieldNew} )

-- (Lax equation solver are used to fix negative pressure/density)
-- updater for electron equations
elcEulerLaxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
elcEulerLaxSlvr:setIn( {elcFluid} )
elcEulerLaxSlvr:setOut( {elcFluidNew} )

-- updater for ion equations
ionEulerLaxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
ionEulerLaxSlvr:setIn( {ionFluid} )
ionEulerLaxSlvr:setOut( {ionFluidNew} )

-- updater for two-fluid sources
sourceSlvrImpl = Updater.ImplicitFiveMomentSrc2D {
   -- grid on which to run updater
   onGrid = grid,
   -- number of fluids
   numFluids = 2,
   -- species charge
   charge = {elcCharge, ionCharge},
   -- species mass
   mass = {elcMass, ionMass},
   -- premittivity of free space
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   -- has static magnetic field
   hasStaticField = false,
}


-- function to apply boundary conditions
function applyBc(fld, t)
   -- sync ghost cells, including periodic BCs
   fld:sync()
end

-- function to update source terms
function calcSourceImpl(elcIn, ionIn, emIn, tCurr, t)
   sourceSlvrImpl:setOut( {elcIn, ionIn, emIn} )
   sourceSlvrImpl:setCurrTime(tCurr)
   sourceSlvrImpl:advance(t)
end

-- function to take one time-step
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   calcSourceImpl(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q)

   -- advance electrons
   elcEulerSlvr:setCurrTime(tCurr)
   local elcStatus, elcDtSuggested = elcEulerSlvr:advance(t)
   -- advance ions
   ionEulerSlvr:setCurrTime(tCurr)
   local ionStatus, ionDtSuggested = ionEulerSlvr:advance(t)
   -- advance fields
   maxSlvr:setCurrTime(tCurr)
   local maxStatus, maxDtSuggested = maxSlvr:advance(t)

   -- check if any updater failed
   local status = false
   local dtSuggested = math.min(elcDtSuggested, ionDtSuggested, maxDtSuggested)
   if (elcStatus and ionStatus and maxStatus) then
      status = true
   end

   -- update source terms
   calcSourceImpl(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew)

   return status, dtSuggested
end

-- function to take one time-step
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   calcSourceImpl(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q)

   -- advance electrons
   elcEulerLaxSlvr:setCurrTime(tCurr)
   local elcStatus, elcDtSuggested = elcEulerLaxSlvr:advance(t)
   -- advance ions
   ionEulerLaxSlvr:setCurrTime(tCurr)
   local ionStatus, ionDtSuggested = ionEulerLaxSlvr:advance(t)
   -- advance fields
   maxSlvr:setCurrTime(tCurr)
   local maxStatus, maxDtSuggested = maxSlvr:advance(t)

   -- check if any updater failed
   local status = false
   local dtSuggested = math.min(elcDtSuggested, ionDtSuggested, maxDtSuggested)
   if (elcStatus and ionStatus and maxStatus) then
      status = true
   end

   -- update source terms
   calcSourceImpl(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew)

   return status, dtSuggested
end

-- compute diagnostic
function calcDiagnostics(tCurr, t)
   -- nothing to do
end

-- advance solution from tStart to tEnd, using optimal time-steps.
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local tfStatus, tfDtSuggested
   local useLaxSolver = false
   while true do
      qDup:copy(q)
      qNewDup:copy(qNew)
      if (tCurr+myDt > tEnd) then -- hit tEnd exactly
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format(" Taking step %d at time %g with dt %g", step, tCurr, myDt))
      if (useLaxSolver) then
	 -- (call Lax solver if positivity violated)
	 tfStatus, tfDtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
	 useLaxSolver = false
      else
	 tfStatus, tfDtSuggested = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (tfStatus == false) then
	 Lucee.logInfo (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = tfDtSuggested
	 qNew:copy(qNewDup)
	 q:copy(qDup)
      elseif ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
	   or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
	 Lucee.logInfo (string.format("** Negative pressure or density at %g! Will retake step with Lax fluxes", tCurr+myDt))
	 q:copy(qDup)
	 qNew:copy(qNewDup)
	 useLaxSolver = true
      else
	 if (qNew:hasNan()) then
	    Lucee.logInfo (string.format(" ** Nan occured at %g! Stopping simulation", tCurr))
	    break
	 end

	 calcDiagnostics(tCurr, tCurr+myDt)
	 q:copy(qNew)

	 myDt = tfDtSuggested
	 tCurr = tCurr + myDt
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end
   end
   
   return tfDtSuggested
end

function writeFrame(frame, tCurr)
   qNew:write(string.format("q_%d.h5", frame), tCurr )
end

-- compute diagnostics and write out initial conditions
calcDiagnostics(0.0, 0.0)
writeFrame(0, 0.0)

dtSuggested = 1.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 10.0

nFrames = 100
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
-- main loop
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   writeFrame(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
