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

-- global parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
lightSpeed = 1.0
epsilon0 = 1.0
mgnErrorSpeedFactor = 1.0

Lx = 100.0
Ly = 50.0
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

NX = 1000
NY = 500

-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {Lx, Ly},
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
-- solution after X update
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 28,
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

elcFluidX = qX:alias(0, 5)
ionFluidX = qX:alias(5, 10)
emFieldX = qX:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-- alias for By
byAlias = qNew:alias(14, 15)

-- function to apply initial conditions
function init(x,y,z)
   local me = elcMass
   local mi = ionMass
   local qe = elcCharge
   local qi = ionCharge
   local gasGamma1 = gasGamma-1
   local psi0 = 0.1*B0

   local pi = Lucee.Pi
   local twopi = 2*pi

   local Lx4 = Lx/4
   local Ly4 = Ly/4

   local numDens = n0*(1/math.cosh((y-Ly4)/lambda))^2 + n0*(1/math.cosh((y-3*Ly4)/lambda))^2 + nb

   -- electron momentum is computed from plasma current that supports field
   local ezmom = -B0*(1/lambda)*(1/math.cosh((y-Ly4)/lambda)^2 - 1/math.cosh((y-3*Ly4)/lambda)^2)*(me/qe)
   local rhoe = numDens*me
   local pre = numDens*B0^2/12.0
   local pri = 5*pre

   -- electron total energy is thermal plus kinetic
   local ere = pre/gasGamma1 + 0.5*ezmom*ezmom/rhoe
   -- ion total energy is thermal: ions do not carry any current
   local eri = pri/gasGamma1

   -- unperturbed field has only Bx component
   local Bx = B0*(-1+math.tanh((y-Ly4)/lambda)-math.tanh((y-3*Ly4)/lambda))
   local dBx1 = -psi0*(pi/Ly)*cos(2*pi*(x-Lx4)/Lx)*sin(pi*(y-Ly4)/Ly)
   local dBy1 = psi0*(2*pi/Lx)*sin(2*pi*(x-Lx4)/Lx)*cos(pi*(y-Ly4)/Ly)
   
   local dBx2 = -psi0*(pi/Ly)*cos(2*pi*(x+Lx4)/Lx)*sin(pi*(y+Ly4)/Ly)
   local dBy2 = psi0*(2*pi/Lx)*sin(2*pi*(x+Lx4)/Lx)*cos(pi*(y+Ly4)/Ly)
   
   -- add perturbation so net field is divergence free
   local Bx = Bx+dBx1+dBx2
   local By = dBy1+dBy2

   return rhoe, 0.0, 0.0, ezmom, ere, rhoi, 0.0, 0.0, 0.0, eri, 0.0, 0.0, 0.0, Bx, By, 0.0, 0.0, 0.0
end
-- set initial conditions for fields and fluids
q:set(init)

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
   -- sync ghost cells (this automatically applies periodic BCs as
   -- grid object has explicit specification of periodic directions)
   fld:sync()
end

-- function to update source terms
function calcSourceImpl(elcIn, ionIn, emIn, tCurr, t)
   sourceSlvrImpl:setIn( {staticEB} )
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

-- dynvector to store integrated flux
byFlux = DataStruct.DynVector { numComponents = 1 }

-- updater for Maxwell equations
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

-- compute diagnostic
function calcDiagnostics(tCurr, t)
   byFluxCalc:setCurrTime(tCurr)
   byFluxCalc:advance(t)
end

-- advance solution from tStart to tEnd, using optimal time-steps.
function advanceFrame(tStart, tEnd, initDt)

   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local tfStatus, tfDtSuggested
   local useLaxSolver = false
   while true do
      -- copy in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
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
	 -- time-step too large
	 Lucee.logInfo (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, tfDtSuggested))
	 myDt = tfDtSuggested
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
end

-- compute diagnostics and write out initial conditions
calcDiagnostics(0.0, 0.0)
writeFrame(0, 0.0)

dtSuggested = 1.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 60.0/wci

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
