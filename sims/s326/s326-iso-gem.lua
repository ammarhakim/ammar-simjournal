-- Program to solve Two-Fluid equations

--[[

Standard GEM challenge problem.

--]]

-- decomposition object to use
decomp = DecompRegionCalc2D.CartGeneral {}

-- global parameters
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
lightSpeed = 1.0
epsilon0 = 1.0
mgnErrorSpeedFactor = 1.0

Lx = 8*Lucee.Pi
Ly = 4*Lucee.Pi
B0 = 0.1
n0 = 1.0
lambda = 0.5
cfl = 0.9
wci = ionCharge*B0/ionMass -- ion cyclotron frequency

Te = B0^2/12.0 -- electron temperature
Ti = 5*Te -- ion temperature

nSpecies = 2

NX = 501
NY = 251

tStart = 0.0 -- start time
tEnd = 40.0/wci -- end time
nFrames = 40 -- number of frames

-- computational domain
grid = Grid.RectCart2D {
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {0},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}
-- solution after X update
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}
-- create duplicate copy in case we need to take step again
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}

-- create aliases to various sub-system
elcFluid = q:alias(0, 4)
ionFluid = q:alias(4, 8)
emField = q:alias(8, 16)

elcFluidX = qX:alias(0, 4)
ionFluidX = qX:alias(4, 8)
emFieldX = qX:alias(8, 16)

elcFluidNew = qNew:alias(0, 4)
ionFluidNew = qNew:alias(4, 8)
emFieldNew = qNew:alias(8, 16)

-- alias for various fields for diagnostics
byAlias = qNew:alias(12, 13)
ezAlias = qNew:alias(10, 11)
neAlias = qNew:alias(0, 1)
uzeAlias = qNew:alias(3, 4)
uziAlias = qNew:alias(7, 8)

-- function to apply initial conditions
function init(x,y,z)
   local me = elcMass
   local mi = ionMass
   local qe = elcCharge
   local qi = ionCharge
   local psi0 = 0.1*B0

   local pi = Lucee.Pi
   local twopi = 2*pi

   -- unperturbed field has only Bx component
   local Bx = B0*math.tanh(y/lambda)
   -- add perturbation so net field is divergence free
   local Bx = Bx - psi0*(pi/Ly)*math.cos(twopi*x/Lx)*math.sin(pi*y/Ly)
   local By = psi0*(twopi/Lx)*math.sin(twopi*x/Lx)*math.cos(pi*y/Ly)

   -- electron momentum is computed from plasma current that supports field
   local ezmom = - (1.0/6.0)*B0*(1/lambda)*(1/math.cosh(y/lambda))^2*(me/qe)
   local rhoe = n0*me*(1/math.cosh(y/lambda))^2 + 0.2*n0*me
   
   local izmom = - (5.0/6.0)*B0*(1/lambda)*(1/math.cosh(y/lambda))^2*(mi/qi)
   local rhoi = n0*mi*(1/math.cosh(y/lambda)^2) + 0.2*n0*mi

   return rhoe, 0.0, 0.0, ezmom, rhoi, 0.0, 0.0, izmom, 0.0, 0.0, 0.0, Bx, By, 0.0, 0.0, 0.0
end
-- set initial conditions for fields and fluids
q:set(init)

-- get ghost cells correct
q:sync()
-- copy initial conditions over
qNew:copy(q)

-- define various equations to solve
elcEulerEqn = HyperEquation.IsothermalEuler {
   -- thermal velocity
   thermalVelocity = math.sqrt(Te/elcMass),
}
ionEulerEqn = HyperEquation.IsothermalEuler {
   -- thermal velocity
   thermalVelocity = math.sqrt(Ti/ionMass),
}
-- (Lax equations are used to fix negative pressure/density)
elcEulerLaxEqn = HyperEquation.IsothermalEuler {
   -- thermal velocity
   thermalVelocity = math.sqrt(Te/elcMass),
   -- use Lax fluxes
   numericalFlux = "lax",   
}
ionEulerLaxEqn = HyperEquation.IsothermalEuler {
   -- thermal velocity
   thermalVelocity = math.sqrt(Ti/ionMass),
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
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidX} )

-- updater for ion equations
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidX} )

-- updater for Maxwell equations
maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldX} )

-- updater for electron equations
elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
elcFluidSlvrDir1:setIn( {elcFluidX} )
elcFluidSlvrDir1:setOut( {elcFluidNew} )

-- updater for ion equations
ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
ionFluidSlvrDir1:setIn( {ionFluidX} )
ionFluidSlvrDir1:setOut( {ionFluidNew} )

-- updater for Maxwell equations
maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
maxSlvrDir1:setIn( {emFieldX} )
maxSlvrDir1:setOut( {emFieldNew} )

-- (Lax equation solver are used to fix negative pressure/density)
-- updater for electron equations
elcLaxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
elcLaxSlvr:setIn( {elcFluid} )
elcLaxSlvr:setOut( {elcFluidNew} )

-- updater for ion equations
ionLaxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- set input/output arrays (these do not change so set it once)
ionLaxSlvr:setIn( {ionFluid} )
ionLaxSlvr:setOut( {ionFluidNew} )

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

-- updater for two-fluid sources
sourceSlvr = Updater.ImplicitFiveMomentSrc2D {
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
   -- we have no pressure equation
   hasPressure = false,
}

-- boundary applicator objects for fluids and fields
bcElcCopy = BoundaryCondition.Copy { components = {0} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcIonCopy = BoundaryCondition.Copy { components = {4} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {5, 6, 7} }
bcElcFld = BoundaryCondition.ZeroTangent { components = {8, 9, 10} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {11, 12, 13} }
bcPot = BoundaryCondition.Copy { components = {14, 15} }

-- top and bottom BC updater
bcBottom = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {
      bcElcCopy, bcElcWall, 
      bcIonCopy, bcIonWall,
      bcElcFld, bcMgnFld, bcPot
   },
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "lower",
}
bcBottom:setOut( {qNew} )

bcTop = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {
      bcElcCopy, bcElcWall, 
      bcIonCopy, bcIonWall,
      bcElcFld, bcMgnFld, bcPot
   },
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "upper",
}
bcTop:setOut( {qNew} )

-- function to apply boundary conditions
function applyBc(fld, t)
   -- walls
   bcBottom:setOut( {fld} )
   bcBottom:advance(t)

   bcTop:setOut( {fld} )
   bcTop:advance(t)
   -- periodic
   fld:sync()
end

-- apply BCs to initial conditions
applyBc(q)
applyBc(qNew)

function updateSource(inpElc, inpIon, inpEm, tCurr, tEnd)
   sourceSlvr:setOut( {inpElc, inpIon, inpEm} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(tEnd)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e6*math.abs(t-tCurr)
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidX) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidX) == false)) then
      -- if positivity violated, return immediatelty
      return false, myDtSuggested
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX)

   -- Y-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir1, ionFluidSlvrDir1, maxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

-- function to take one time-step
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source term
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndField(tCurr, t)

   if (status) then
      -- update source terms
      updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
      applyBc(qNew)
   end

   return status, dtSuggested
end

-- function to take one time-step
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source term
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q)

   -- advance electrons
   elcLaxSlvr:setCurrTime(tCurr)
   local elcStatus, elcDtSuggested = elcLaxSlvr:advance(t)
   -- advance ions
   ionLaxSlvr:setCurrTime(tCurr)
   local ionStatus, ionDtSuggested = ionLaxSlvr:advance(t)
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
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew)

   return status, dtSuggested
end

-- dynvector to store integrated flux
byFlux = DataStruct.DynVector { numComponents = 1 }
byFluxCalc = Updater.IntegrateFieldAlongLine2D {
   onGrid = grid,
   -- start cell
   startCell = {0, NY/2},
   -- direction to integrate in
   dir = 0,
   -- number of cells to integrate
   numCells = NX,
   -- integrand
   integrand = function (by)
		  return math.abs(by)
	       end,
}
byFluxCalc:setIn( {byAlias} )
byFluxCalc:setOut( {byFlux} )

-- dynvector to store Ez at X-point
xpointEz = DataStruct.DynVector { numComponents = 1 }
xpointEzRec = Updater.RecordFieldInCell2D {
   onGrid = grid,
   -- index of cell to record
   cellIndex = {(NX-1)/2, (NY-1)/2},
}
xpointEzRec:setIn( {ezAlias} )
xpointEzRec:setOut( {xpointEz} )

-- dynvector to store number density at X-point
xpointNe = DataStruct.DynVector { numComponents = 1 }
xpointNeRec = Updater.RecordFieldInCell2D {
   onGrid = grid,
   -- index of cell to record
   cellIndex = {(NX-1)/2, (NY-1)/2},
}
xpointNeRec:setIn( {neAlias} )
xpointNeRec:setOut( {xpointNe} )

-- dynvector to store electron uz at X-point
xpointUze = DataStruct.DynVector { numComponents = 1 }
xpointUzeRec = Updater.RecordFieldInCell2D {
   onGrid = grid,
   -- index of cell to record
   cellIndex = {(NX-1)/2, (NY-1)/2},
}
xpointUzeRec:setIn( {uzeAlias} )
xpointUzeRec:setOut( {xpointUze} )

-- dynvector to store ion uz at X-point
xpointUzi = DataStruct.DynVector { numComponents = 1 }
xpointUziRec = Updater.RecordFieldInCell2D {
   onGrid = grid,
   -- index of cell to record
   cellIndex = {(NX-1)/2, (NY-1)/2},
}
xpointUziRec:setIn( {uziAlias} )
xpointUziRec:setOut( {xpointUzi} )

-- compute diagnostic
function calcDiagnostics(tCurr, t)
   for i,diag in ipairs({byFluxCalc, xpointEzRec, xpointNeRec, xpointUzeRec, xpointUziRec}) do
      diag:setCurrTime(tCurr)
      diag:advance(t)
   end
end

function writeFields(frame, tCurr)
   qNew:write(string.format("q_%d.h5", frame), tCurr )
   byFlux:write( string.format("byFlux_%d.h5", frame) )
   xpointEz:write(string.format("xpointEz_%d.h5", frame) )
   xpointNe:write(string.format("xpointNe_%d.h5", frame) )
   xpointUze:write(string.format("xpointUze_%d.h5", frame) )
   xpointUzi:write(string.format("xpointUzi_%d.h5", frame) )
end

-- check if to trigger data output
function triggerTurnstile(gate, t)
   if (t >= gate.trigger) then
      gate.trigger = gate.trigger + gate.increment
      gate.open = true
   else
      gate.open = false
   end
   return gate
end

-- advance solution from tStart to tEnd, using optimal time-steps.
function runSimulation(tStart, tEnd, nFrames, initDt)

   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local myGate = { trigger = tFrame, open = false, increment = tFrame }
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested
   local useLaxSolver = false
   while true do
      qDup:copy(q)
      qNewDup:copy(qNew)

      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format(" Taking step %d at time %g with dt %g", step, tCurr, myDt))
      -- advance fluids and fields
      if (useLaxSolver) then
	 -- (call Lax solver if positivity violated)
	 status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
	 useLaxSolver = false
      else
	 status, dtSuggested = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
	 -- time-step too large
	 Lucee.logInfo (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
	 q:copy(qDup)
      elseif ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
	   or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
	 -- negative density/pressure occured
	 Lucee.logInfo (string.format("** Negative pressure or density at %g! Will retake step with Lax fluxes", tCurr+myDt))
	 q:copy(qDup)
	 qNew:copy(qNewDup)
	 useLaxSolver = true
      else
	 -- check if a nan occured
	 if (qNew:hasNan()) then
	    Lucee.logInfo (string.format(" ** Nan occured at %g! Stopping simulation", tCurr))
	    break
	 end
	 -- compute diagnostics
	 calcDiagnostics(tCurr, tCurr+myDt)
	 -- copy updated solution back
	 q:copy(qNew)

	 myGate = triggerTurnstile(myGate, tCurr+myDt)
	 if (myGate.open) then
	    -- write out data
	    Lucee.logInfo (string.format("Writing data at time %g ...\n", tCurr+myDt))
	    writeFields(frame, tCurr+myDt)
	    frame = frame + 1
	    step = 0
	 end

	 tCurr = tCurr + myDt
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end
   end
   
   return dtSuggested
end

-- write initial conditions
writeFields(0, 0.0)

-- write initial condition
writeFields(0, 0.0)
-- run simulation
dtSuggested = 1.0 -- initial time-step, will be adjusted
runSimulation(tStart, tEnd, nFrames, dtSuggested)
