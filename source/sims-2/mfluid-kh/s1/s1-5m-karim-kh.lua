-- Turbulence Dissipation Challenge Kelvin-Helmholtz multi-fluid
-- setup. This is based on PIC setup in
--
-- Karimabadi, H., Roytershteyn, V., Wan, M., Matthaeus, W. H.,
-- Daughton, W., Wu, P., et al. (2013). Coherent structures,
-- intermittent turbulence, and dissipation in high-temperature
-- plasmas. Physics of Plasmas, 20(1), 012303. doi:10.1063/1.4773205
--

Pi = Lucee.Pi
log = Lucee.logInfo

-- physical parameters
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)
mgnErrorSpeedFactor = 1.0

-- problem parameters
n0 = 1.0 -- initial number density
theta = 2.86*Pi/180 -- inclination angle of field
plasmaBeta = 0.1 -- plasma beta
wpe_OmegaCe = 2.0
Lvfac = 4 -- Lv = Lvfac*di

-- secondary quantities
wpe = math.sqrt(n0*elcCharge^2/(epsilon0*elcMass))
wpi = math.sqrt(n0*ionCharge^2/(epsilon0*ionMass))
di = lightSpeed/wpi

B0 = wpe*elcMass/math.abs(elcCharge)/wpe_OmegaCe
Valf = B0*math.sin(theta)/math.sqrt(mu0*n0*ionMass)
U0 = 10*Valf
T0 = 0.5*plasmaBeta*(B0^2/(2*mu0))/n0 -- Te=Ti=T

Lv = Lvfac*di -- in physical units
Ly = 4*Pi*Lv/0.5 -- this ensures exactly 2 modes in domain
Lx = Ly/2.0

TiOverTe = 1.0

OmegaCe0 = elcCharge*B0/elcMass
OmegaCi0 = ionCharge*B0/ionMass

-- resolution and time-stepping
NX = 500
NY = 1000
cfl = 0.9
tStart = 0.0
tEnd = 500/OmegaCi0
nFrames = 100

log(string.format("di=%g", di))
log(string.format("wpe/OmegaCe=%g", -wpe/OmegaCe0))
log(string.format("OmegaCi=%g", OmegaCi0))
log(string.format("Lx=%g, Ly=%g", Lx, Ly))
log(string.format("plasmaBeta=%g", plasmaBeta))
log(string.format("Vthe/c=%g", math.sqrt(2*T0/elcMass)/lightSpeed))
log(string.format("Valf/c=%g", Valf/lightSpeed))
log(string.format("tEnd=%g,  nFrames=%d",tEnd,nFrames))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {-Lx/2, 0.0},
   upper = {Lx/2, Ly},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {1},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- aliases to various sub-systems
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidX = qX:alias(0, 5)
ionFluidX = qX:alias(5, 10)
emFieldX = qX:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   local tanh = math.tanh
   local cos = math.cos
   local sin = math.sin

   local me = elcMass
   local mi = ionMass
   local qe = elcCharge
   local qi = ionCharge
   local g1 = gasGamma-1.0
   local beta = plasmaBeta
   local Te = T0
   local Ti = T0

   local Bx = 0.0
   local By = B0*sin(theta)
   local Bz = B0*cos(theta)
   local Ux = 0.0 -- 0.15*U0*sin(0.5*y/Lv)*math.exp(-x^2/Lv^2)
   local Uy = U0*tanh(x/Lv) + 0.15*U0*sin(0.5*y/Lv)*math.exp(-x^2/Lv^2)

   local Ex = -B0*U0*tanh(x/Lv)
   local Ey = 0.0
   local Ez = 0.0

   local rhoe = me*n0
   local ee = n0*Te/g1 + 0.5*rhoe*(Ux^2 + Uy^2)
   
   local rhoi = mi*n0
   local ei = n0*Ti/g1 + 0.5*rhoi*(Ux^2 + Uy^2)

   return rhoe, rhoe*Ux, rhoe*Uy, 0.0, ee, rhoi, rhoi*Ux, rhoi*Uy, 0.0, ei, Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

bcElcCopy = BoundaryCondition.Copy { components = {0, 4} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcIonCopy = BoundaryCondition.Copy { components = {5, 9} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {6, 7, 8} }
bcElcFld = BoundaryCondition.ZeroTangent { components = {10, 11, 12} }
bcMgnFld = BoundaryCondition.ZeroNormal { components = {13, 14, 15} }
bcPot = BoundaryCondition.Copy { components = {16, 17} }
--FIXME: fact in bcPot

-- create boundary condition object
function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {
	 bcElcCopy, bcElcWall, 
	 bcIonCopy, bcIonWall,
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
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-- ds solvers for regular Euler equations along X
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
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

-- ds solvers for regular Euler equations along Y
elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
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

-- ds solvers for Lax Euler equations along X
elcLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
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

-- ds solvers for Lax Euler equations along Y
elcLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
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
sourceSlvr = Updater.ImplicitFiveMomentSrc2D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
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
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidX) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidX) == false)) then
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

   if ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
       useLaxSolver = true
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



