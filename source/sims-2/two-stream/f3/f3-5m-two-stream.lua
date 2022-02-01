-- Gkeyll
-- Two-fluid 5m simulation
-- Two-stream instability

----------------------------------------------------------------------
-- Simulation parameters ---------------------------------------------
pi = Lucee.Pi
-- Constatnts
gasGamma = 5/3 --2.0
chargeElc = -1.0
massElc = 1.0
epsilon0 = 1/25.0
mu0 = 1.0
lightSpeed = 1.0/math.sqrt(epsilon0*mu0)

k0 = 0.5 -- wave-number
vDrift = 1.0 -- drift velocity
perturbation = 1.0e-6 -- distribution function perturbation

-- Initial conditions
nElc10 = 0.5
nElc20 = 0.5
uxElc10 = -vDrift
uyElc10 = 0.0
uzElc10 = 0.0
uxElc20 = vDrift
uyElc20 = 0.0
uzElc20 = 0.0
vthElc10 = 0.1
vthElc20 = 0.1
-- IC automatically calculated
TElc10 = massElc*vthElc10^2
TElc20 = massElc*vthElc20^2

n0 = nElc10+nElc20
lambdaD = math.sqrt(epsilon0*TElc10/(n0*chargeElc^2))

-- Domain and time
numCells = {200, 1} -- Nx, Ny
lowerBoundary = {-lambdaD*10*Lucee.Pi/k0, 0.0} -- xLow, yLow
upperBoundary = {lambdaD*10*Lucee.Pi/k0, 1.0} -- xUp, yUp
cfl = 0.9
tStart = 0.0
tEnd = 20.0
numFrames = 100

-- Other setup
numFluids = 2
Lucee.IsRestarting = false
Lucee.RestartFrame = -1
limiter = "van-leer"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 0

----------------------------------------------------------------------
-- Initial conditions ------------------------------------------------
function init(x, y, z)
   local alpha = perturbation
   
   local rho1 = massElc*nElc10*(1+alpha*math.cos(k0*x))
   local mx1 = rho1*uxElc10
   local my1 = rho1*uyElc10
   local mz1 = rho1*uzElc10
   local p1 = nElc10*TElc10/(gasGamma-1) +
      0.5*(mx1*mx1 + my1*my1 + mz1*mz1)/rho1


   local rho2 = massElc*nElc20*(1+alpha*math.cos(k0*x))
   local mx2 = rho2*uxElc20
   local my2 = rho2*uyElc20
   local mz2 = rho2*uzElc20
   local p2 = nElc20*TElc20/(gasGamma-1) +
      0.5*(mx2*mx2 + my2*my2 + mz2*mz2)/rho2

   local Ex = 0.0 --math.abs(chargeElc)/epsilon0*(nElc10+nElc20)*alpha*math.sin(k0*x)/k0
   local Ey = 0.0
   local Ez = 0.0
   local Bx = 0.0
   local By = 0.0
   local Bz = 0.0

   return rho1, mx1, my1, mz1, p1, rho2, mx2, my2, mz2, p2, Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
end

----------------------------------------------------------------------
-- Convenience and Debugging -----------------------------------------
log = function(...) Lucee.logInfo(string.format(...)) end

function HERE()
   local info = debug.getinfo(2)
   str = string.format("HERE: %d %s: %s", info.currentline,
		       info.source, tostring(info.name))
   Lucee.logInfo(string.format(str))
end

----------------------------------------------------------------------
-- Setup verification ------------------------------------------------
dtL = (upperBoundary[1]-lowerBoundary[1])/numCells[1]*cfl/lightSpeed

log("Setup verification:")
log(" * Speed of light = %g", lightSpeed)
log(" * Vthe1/c = %g", vthElc10/lightSpeed)
log(" * Vthe2/c = %g", vthElc20/lightSpeed)
log(" * Time-step from light-speed = %g", dtL)
log(" * tEnd = %g,  numFrames = %d", tEnd, numFrames)

----------------------------------------------------------------------
-- Computational domain, Data fields ---------------------------------
decomp = DecompRegionCalc2D.CartGeneral {}
grid = Grid.RectCart2D {
   lower = lowerBoundary,
   upper = upperBoundary,
   cells = numCells,
   decomposition = decomp,
   periodicDirs = {0, 1},
}

createData = function(numComponents)
   return DataStruct.Field2D {
      onGrid = grid,
      numComponents = numComponents,
      ghost = {2, 2},
   }
end
q = createData(18)
qTemp = createData(18)
qNew = createData(18)
qDup = createData(18)

-- aliases to various sub-systems
getFields = function(fld)
   return fld:alias(0, 5), fld:alias(5, 10), fld:alias(10, 18)
end
elcFluid1, elcFluid2, emField = getFields(q)
elcFluid1Temp, elcFluid2Temp, emFieldTemp = getFields(qTemp)
elcFluid1New, elcFluid2New, emFieldNew = getFields(qNew)

qIn = {q, qTemp}
qOut = {qTemp, qNew}

elcFluid1Out = {elcFluid1Temp, elcFluid1New}
elcFluid2Out = {elcFluid2Temp, elcFluid2New}

----------------------------------------------------------------------
-- Boundary conditions ----------------------------------------------- 
function applyBc(fld, t, dt)
   -- copy BCs in Y
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
   fld:sync()
end


----------------------------------------------------------------------
-- Equation solvers --------------------------------------------------
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
eulerEqnLax = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

createSlvr = function(eqn, input, output, dir, limiter)
   local slvr = Updater.WavePropagation2D {
      onGrid = grid,
      equation = eqn,
      limiter = limiter,
      cfl = cfl,
      cflm = 1.1*cfl,
      updateDirections = {dir}
   }
   slvr:setIn( {input} )
   slvr:setOut( {output} )
   return slvr
end
elc1SlvrDir0 = createSlvr(eulerEqn, elcFluid1, elcFluid1Temp, 0, limiter)
elc2SlvrDir0 = createSlvr(eulerEqn, elcFluid2, elcFluid2Temp, 0, limiter)
emSlvrDir0 = createSlvr(maxwellEqn, emField, emFieldTemp, 0, limiter)

elc1SlvrDir1 = createSlvr(eulerEqn, elcFluid1Temp, elcFluid1New, 1, limiter)
elc2SlvrDir1 = createSlvr(eulerEqn, elcFluid2Temp, elcFluid2New, 1, limiter)
emSlvrDir1 = createSlvr(maxwellEqn, emFieldTemp, emFieldNew, 1, limiter)

elc1SlvrDir0Lax = createSlvr(eulerEqnLax, elcFluid1, elcFluid1Temp, 0, "zero")
elc2SlvrDir0Lax = createSlvr(eulerEqnLax, elcFluid2, elcFluid2Temp, 0, "zero")
emSlvrDir0Lax = createSlvr(maxwellEqn, emField, emFieldTemp, 0, "zero")

elc1SlvrDir1Lax = createSlvr(eulerEqnLax, elcFluid1Temp, elcFluid1New, 1,"zero")
elc2SlvrDir1Lax = createSlvr(eulerEqnLax, elcFluid2Temp, elcFluid2New, 1,"zero")
emSlvrDir1Lax = createSlvr(maxwellEqn, emFieldTemp, emFieldNew, 1, "zero")

slvrs = {
   {elc1SlvrDir0, elc2SlvrDir0, emSlvrDir0},
   {elc1SlvrDir1, elc2SlvrDir1, emSlvrDir1},
}

slvrsLax = {
   {elc1SlvrDir0Lax, elc2SlvrDir0Lax, emSlvrDir0Lax},
   {elc1SlvrDir1Lax, elc2SlvrDir1Lax, emSlvrDir1Lax},
}

-- updater for source terms
sourceSlvr = Updater.ImplicitFiveMomentSrc2D {
   onGrid = grid,
   numFluids = 2,
   charge = {chargeElc, chargeElc},
   mass = {massElc, massElc},
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, t)
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
   for i,slvr in ipairs({elc1SlvrDir0, elc2SlvrDir0, emSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((eulerEqn:checkInvariantDomain(elcFluid1Temp) == false)
    or (eulerEqn:checkInvariantDomain(elcFluid2Temp) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qTemp, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elc1SlvrDir1, elc2SlvrDir1, emSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((eulerEqn:checkInvariantDomain(elcFluid1New) == false)
    or (eulerEqn:checkInvariantDomain(elcFluid2New) == false)) then
       useLaxSolver = true
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid1, elcFluid2, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(elcFluid1New, elcFluid2New, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({elc1SlvrDir0Lax, elc2SlvrDir0Lax, emSlvrDir0Lax}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   applyBc(qTemp, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elc1SlvrDir1Lax, elc2SlvrDir1Lax, emSlvrDir1Lax}) do
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
   updateSource(elcFluid1, elcFluid2, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(elcFluid1New, elcFluid2New, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------------------------------------------------
-- Diagnostics and Output --------------------------------------------
emEnergy = DataStruct.DynVector { numComponents = 1 }
emEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5*epsilon0*(ex*ex + ey*ey + ez*ez) +
		     0.5/mu0*(bx*bx + by*by + bz*bz)
	       end,
}
emEnergyCalc:setIn( {emField} )
emEnergyCalc:setOut( {emEnergy} )

ExEnergy = DataStruct.DynVector { numComponents = 1 }
ExEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5*epsilon0*ex*ex
	       end,
}
ExEnergyCalc:setIn( {emField} )
ExEnergyCalc:setOut( {ExEnergy} )

EyEnergy = DataStruct.DynVector { numComponents = 1 }
EyEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5*epsilon0*ey*ey
	       end,
}
EyEnergyCalc:setIn( {emField} )
EyEnergyCalc:setOut( {EyEnergy} )

BzEnergy = DataStruct.DynVector { numComponents = 1 }
BzEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5/mu0*bz*bz
	       end,
}
BzEnergyCalc:setIn( {emField} )
BzEnergyCalc:setOut( {BzEnergy} )

totalPtclEnergy = DataStruct.DynVector { numComponents = 1 }
totalPtclEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (rho1, mx1, my1, mz1, p1, rho2, mx2, my2, mz2, p2, ex, ey, ez, bx, by, bz, e1, e2)
		  return p1 + p2
	       end,
}
totalPtclEnergyCalc:setIn( {qNew} )
totalPtclEnergyCalc:setOut( {totalPtclEnergy} )

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({emEnergyCalc, ExEnergyCalc, EyEnergyCalc,
			 BzEnergyCalc, totalPtclEnergyCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )

   emEnergy:write( string.format("emEnergy_%d.h5", frame) )
   ExEnergy:write( string.format("ExEnergy_%d.h5", frame) )
   EyEnergy:write( string.format("EyEnergy_%d.h5", frame) )
   BzEnergy:write( string.format("BzEnergy_%d.h5", frame) )
   totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame) )
end

----------------------------
-- TIME-STEPPING FUNCTION --
----------------------------
function runSimulation(tStart, tEnd, numFrames, initDt)

   local frame = 1
   local tFrame = (tEnd-tStart)/numFrames
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
        q:copy(qDup)
      elseif (useLaxSolver == true) then
        -- negative density/pressure occured
        log (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
        q:copy(qDup)
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

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, numFrames, initDt)


