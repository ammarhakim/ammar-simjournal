-- Program to solve Two-Fluid equations

--[[

According to Daughton et. al. 2006 PoP paper one has
 
rhoi/L = 1.0
mi/me = 25
Ti/Te = 5
wpe/wce = 3
nb/n0 = 0.3

--]]

-- decomposition object to use
decomp = DecompRegionCalc2D.CartGeneral {}

elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
lightSpeed = 1.0
epsilon0 = 1.0
mgnErrorSpeedFactor = 1.0
elcCollisionFreq = 1.0

Lx = 50.0
Ly = 25.0
B0 = 1/15.0
n0 = 1.0
nb = 0.3*n0
lambda = math.sqrt(10/12)
cfl = 0.9
bGuideFactor = 0.0
wci = ionCharge*B0/ionMass -- ion cyclotron frequency
elcPlasmaFreq = math.sqrt(n0*elcCharge^2/(epsilon0*elcMass)) -- plasma frequency
elcSkinDepth = lightSpeed/elcPlasmaFreq

nSpecies = 2

NX = 1001
NY = 501

-- A generic function to run an updater.
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

-- computational domain
grid = Grid.RectCart2D {
   lower = {-Lx/2, -Ly/2},
   upper = {Lx/2, Ly/2},
   cells = {NX, NY},
   decomposition = decomp,
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- solution after X update
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}
-- create duplicate copy in case we need to take step again
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
   ghost = {2, 2},
}

-- create aliases to various sub-system
elcFluid = q:alias(0, 10)
ionFluid = q:alias(10, 20)
emField = q:alias(20, 28)

elcFluidX = qX:alias(0, 10)
ionFluidX = qX:alias(10, 20)
emFieldX = qX:alias(20, 28)

elcFluidNew = qNew:alias(0, 10)
ionFluidNew = qNew:alias(10, 20)
emFieldNew = qNew:alias(20, 28)

-- alias for various fields for diagnostics
byAlias = qNew:alias(24, 25)
ezAlias = qNew:alias(22, 23)
neAlias = qNew:alias(0, 1)
uzeAlias = qNew:alias(3, 4)
niAlias = qNew:alias(10, 11)
uziAlias = qNew:alias(13, 14)

-- function to apply initial conditions
function initElc(x,y,z)
   local me = elcMass
   local qe = elcCharge

  local numDens = n0*(1/math.cosh(y/lambda))^2 + nb

   -- electron momentum is computed from plasma current that supports field
   local ezmom = -(1.0/6.0)*B0*(1/lambda)*(1/math.cosh(y/lambda))^2*(me/qe)
   local rhoe = numDens*me
   local pre = numDens*B0^2/12.0
   local pzz = pre + ezmom*ezmom/rhoe
   
   return rhoe, 0.0, 0.0, ezmom, pre, 0.0, 0.0, pre, 0.0, pzz
end
-- set electron initial conditions
elcFluid:set(initElc)

function initIon(x,y,z)
   local mi = ionMass
   local qi = ionCharge

   local izmom = -(5.0/6.0)*B0*(1/lambda)*(1/math.cosh(y/lambda))^2*(mi/qi)
   local numDens = n0*(1/math.cosh(y/lambda))^2 + nb
   local rhoi = numDens*mi
   local pre = numDens*B0^2/12.0
   local pri = 5*pre
   local pzz = pri + izmom*izmom/rhoi

   return rhoi, 0.0, 0.0, izmom, pri, 0.0, 0.0, pri, 0.0, pzz
end
-- set ion initial conditions
ionFluid:set(initIon)

function initEmField(x,y,z)
   local psi0 = 0.1*B0
   local pi = Lucee.Pi
   local twopi = 2*pi

   -- unperturbed field has only Bx component
   local Bx = B0*math.tanh(y/lambda)
   -- add perturbation so net field is divergence free
   local Bx = Bx - psi0*(pi/Ly)*math.cos(twopi*x/Lx)*math.sin(pi*y/Ly)
   local By = psi0*(twopi/Lx)*math.sin(twopi*x/Lx)*math.cos(pi*y/Ly)

   return 0.0, 0.0, 0.0, Bx, By, 0.0, 0.0, 0.0
end
-- set ion initial conditions
emField:set(initEmField)

-- make sure ghosts are set correctly
q:sync()
-- copy initial conditions over
qNew:copy(q)

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

-- updater for electron equations
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidX} )

elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEqn,
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
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidX} )

ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqn,
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
   equation = elcLaxEqn,
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
   equation = ionLaxEqn,
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

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.ImplicitTenMomentSrc2D {
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
}

-- Collisional source updaters
elcCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = 1/elcSkinDepth,
}
ionCollSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = 1/elcSkinDepth,
}

-- boundary applicator objects for fluids and fields
bcElcCopy = BoundaryCondition.Copy { components = {0, 4} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcIonCopy = BoundaryCondition.Copy { components = {5, 9} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {6, 7, 8} }
bcElcFld = BoundaryCondition.ZeroTangent { components = {10, 11, 12} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {13, 14, 15} }
bcPot = BoundaryCondition.Copy { components = {16, 17} }

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
   -- in X
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   -- in Y
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
   -- sync ghost cells
   fld:sync()
end

function updateSource(inpElc, inpIon, inpEm, tCurr, tEnd)
   -- EM sources
   sourceSlvr:setOut( {inpElc, inpIon, inpEm} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(tEnd)

   -- electron collisional relaxation
   elcCollSrcSlvr:setOut( {inpElc} )
   elcCollSrcSlvr:setCurrTime(tCurr)
   elcCollSrcSlvr:advance(tEnd)

   -- ion collisional relaxation
   ionCollSrcSlvr:setOut( {inpIon} )
   ionCollSrcSlvr:setCurrTime(tCurr)
   ionCollSrcSlvr:advance(tEnd)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
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

-- apply BCs to initial conditions
applyBc(q)
applyBc(qNew)

-- function to take one time-step
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source term
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew)

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
-- dynvectors stuff at X-point
xpointNe = DataStruct.DynVector { numComponents = 1 }
xpointUze = DataStruct.DynVector { numComponents = 1 }
xpointNi = DataStruct.DynVector { numComponents = 1 }
xpointUzi = DataStruct.DynVector { numComponents = 1 }
xpointEz = DataStruct.DynVector { numComponents = 1 }

-- updater for Maxwell equations
byFluxCalc = Updater.IntegrateFieldAlongLine2D {
   onGrid = grid,
   -- start cell
   startCell = {0, (NY-1)/2},
   -- direction to integrate in
   dir = 0,
   -- number of cells to integrate
   numCells = NX,
   -- integrand
   integrand = function (by)
		  return math.abs(by)
	       end,
}
-- updater to record stuff at X-point
xpointRec = Updater.RecordFieldInCell2D {
   onGrid = grid,
   -- index of cell to record
   cellIndex = {(NX-1)/2, (NY-1)/2},
}

-- compute diagnostic
function calcDiagnostics(tCurr, t)
   local dt = t-tCurr
   runUpdater(byFluxCalc, tCurr, dt, {byAlias}, {byFlux})
   runUpdater(xpointRec, tCurr, dt, {ezAlias}, {xpointEz})
   runUpdater(xpointRec, tCurr, dt, {neAlias}, {xpointNe})
   runUpdater(xpointRec, tCurr, dt, {uzeAlias}, {xpointUze})
   runUpdater(xpointRec, tCurr, dt, {niAlias}, {xpointNi})
   runUpdater(xpointRec, tCurr, dt, {uziAlias}, {xpointUzi})
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

      if (tCurr+myDt > tEnd) then -- adjust to hit tEnd exactly
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format(" Taking step %d at time %g with dt %g", step, tCurr, myDt))
      -- advance fluids and fields
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
      elseif ((elcEqn:checkInvariantDomain(elcFluidNew) == false)
	   or (ionEqn:checkInvariantDomain(ionFluidNew) == false)) then
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
   -- write out data
   qNew:write(string.format("q_%d.h5", frame), tCurr )
   byFlux:write( string.format("byFlux_%d.h5", frame) )
   xpointEz:write( string.format("xpointEz_%d.h5", frame) )
   xpointNe:write( string.format("xpointNe_%d.h5", frame) )
   xpointUze:write( string.format("xpointUze_%d.h5", frame) )
   xpointNi:write( string.format("xpointNi_%d.h5", frame) )
   xpointUzi:write( string.format("xpointUzi_%d.h5", frame) )
end

-- compute diagnostics and write out initial conditions
calcDiagnostics(0.0, 0.0)
writeFrame(0, 0.0)

dtSuggested = 1.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 120.0/wci

nFrames = 60
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
