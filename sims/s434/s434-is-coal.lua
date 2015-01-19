-- Island coalescence problem

Pi = Lucee.Pi
log = Lucee.logInfo

-- physical parameters
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)
mgnErrorSpeedFactor = 1.0

n0 = 1.0
wpe = math.sqrt(n0*elcCharge^2/(epsilon0*elcMass))
wpi = math.sqrt(n0*ionCharge^2/(epsilon0*ionMass))
de = lightSpeed/wpe
di = lightSpeed/wpi

-- average wave number for pressure damping
elcCollAverageWaveNumber = 1/de
ionCollAverageWaveNumber = 1/de

lambda = 5.0*di
Lx = 4*Pi*lambda
Ly = Lx/2
Valf = 0.1
plasmaBeta = 1.0
TiOverTe = 1.0
nbOverN0 = 0.2
pert = 0.1
islandWidth = 0.4
vte = lightSpeed*0.35
T0 = vte^2*elcMass/2

B0 = Valf*math.sqrt(mu0*n0*ionMass)
OmegaCe0 = elcCharge*B0/elcMass
OmegaCi0 = ionCharge*B0/ionMass
psi0 = pert*B0

-- resolution and time-stepping
NX = 1280
NY = 640
cfl = 0.9
tStart = 0.0
tEnd = 2.5*Lx/Valf
nFrames = 50

log(string.format("di=%g", di))
log(string.format("wpe/OmegaCe=%g", -wpe/OmegaCe0))
log(string.format("Lx=%gdi", Lx/di))
log(string.format("plasmaBeta=%g", plasmaBeta))
log(string.format("Valf/c=%g", Valf/lightSpeed))
log(string.format("Vthe/c=%g", math.sqrt(2*T0/elcMass)/lightSpeed))
log(string.format("tEnd=%g,  nFrames=%d",tEnd,nFrames))

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
   periodicDirs = {0},
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
-- initial conditions
function init(x,y,z)
   local tanh = math.tanh
   local cosh = math.cosh
   local sinh = math.sinh
   local cos = math.cos
   local sin = math.sin

   local me = elcMass
   local mi = ionMass
   local qe = elcCharge
   local qi = ionCharge
   local ep = islandWidth
   local beta = plasmaBeta
   local l = lambda
   local TeFrac = 1.0 / (1.0 + TiOverTe)
   local TiFrac = 1.0 - TeFrac
   local sech2 = (1.0/cosh(y/l))^2
   local _2pi = 2.0*Pi
   local denom = ep*cos(x/l) + cosh(y/l)

   local Bxb = B0*sinh(y/l)/denom
   local Byb = B0*ep*sin(x/l)/denom
   local Bx = Bxb - psi0*(Lx/Ly/2)*cos(_2pi*(x-.5*Lx)/Lx)*sin(Pi*y/Ly) 
   local By = Byb + psi0*sin(_2pi*(x-.5*Lx)/Lx)*cos(Pi*y/Ly)
   local Bz = 0.0

   -- n = n0*sech^2(y/lambda)
   -- J = -(B0/lambda)*sech^2(y/lambda)
   -- assume Je/Ji = Te/Ti >>> Je = J*TeFrac, Ji = J*TiFrac
   -- emom = ne*me*ue = (me/qe)*Je
   -- n0*Ttotal = B0^2/2, Te = Ttotal*TeFrac, Ti = Ttotal*TiFrac
   local n = n0*nbOverN0 + n0*(1.0-ep*ep)/(ep*cos(x/l)+cosh(y/l))^2
   local Ttotal = beta*(B0*B0)/2.0/n0
   local Te = Ttotal*TeFrac
   local Ti = Ttotal*TiFrac
   local Jz  = (B0/l)*(ep*cos(x/l)/denom - cosh(y/l)/denom + ep*ep*sin(x/l)*sin(x/l)/denom/denom + sinh(y/l)*sinh(y/l)/denom/denom ) + psi0*Pi*(Lx/Ly^2/2.0 + 2/Lx)*cos(_2pi *(x - .5*Lx)/Lx)*cos(Pi*y/Ly)
   local Jze = Jz*TeFrac
   local Jzi = Jz*TiFrac

   local rhoe = me*n
   local momze = (me/qe)*Jze
   local pe = n*Te
   local pezz =  n*Te + momze*momze/rhoe
   
   local rhoi = mi*n
   local momzi = (mi/qi)*Jzi
   local pi = n*Ti
   local pizz =  n*Te + momzi*momzi/rhoi

   return rhoe, 0.0, 0.0, momze, pe, 0.0, 0.0, pe, 0.0, pezz, rhoi, 0.0, 0.0, momzi, pi, 0.0, 0.0, pi, 0.0, pizz, 0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields
bcElcCopy = BoundaryCondition.Copy { components = {0} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcElcPrCopyY = BoundaryCondition.Copy { components = {4, 6, 7, 9} }
bcElcPrFlipY = BoundaryCondition.Copy { components = {5, 8}, fact = {-1, -1} }

bcIonCopy = BoundaryCondition.Copy { components = {10} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {11, 12, 13} }
bcIonPrCopyY = BoundaryCondition.Copy { components = {14, 16, 17, 19} }
bcIonPrFlipY = BoundaryCondition.Copy { components = {15, 18}, fact = {-1, -1} }

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
bcBottom = createBc(1, "lower")
bcTop = createBc(1, "upper")

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   for i,bc in ipairs({bcBottom, bcTop}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
elcEqn = HyperEquation.TenMoment {
}
ionEqn = HyperEquation.TenMoment {
}
-- (Lax equations are used to fix negative pressure/density)
elcEqnLaxEqn = HyperEquation.TenMoment {
   numericalFlux = "lax",   
}
ionEqnLaxEqn = HyperEquation.TenMoment {
   numericalFlux = "lax",
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
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
   equation = elcEqnLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqnLaxEqn,
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
   equation = elcEqnLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEqnLaxEqn,
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
    or (qX:hasNan())
    or (ionEqn:checkInvariantDomain(ionFluidX) == false)) then
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
    or (qNew:hasNan())
    or (ionEqn:checkInvariantDomain(ionFluidNew) == false)) then
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

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------
-- dynvector to store integrated flux

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


