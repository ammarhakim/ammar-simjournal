-- Neutral flow over 2D cylinder (Mercury)

log = Lucee.logInfo

-- physical parameters
gasGamma = 5.0/3.0
ionMass = Lucee.ProtonMass -- ion Mass [kg]

-- Simulation parameters
swTemp = 15*Lucee.Ev2Kelvin -- SW temperature [K]
swNumDensity = 40e6 -- [#/m^3]
swDensity = swNumDensity*ionMass -- [kg/m^3]
swPressure = swNumDensity*Lucee.BoltzmannConstant*swTemp -- [Pa]
swSpeed = {-400e3, 50e3, 0.0 } -- [m/s]

-- Radius of Mercury (used as lenght scale)
Rm = 2439.7e3 -- [m]
-- Domain shape
--LX = {-64.0*Rm, 24.0*Rm} -- [m]
--LY = {-32.0*Rm, 32.0*Rm} -- [m]

LX = {-10.0*Rm, 10.0*Rm} -- [m]
LY = {-10.0*Rm, 10.0*Rm} -- [m]

-- Transit time for solar wind across domain
tSwTransit = (LX[2]-LX[1])/math.abs(swSpeed[1])

-- Resolution, time-stepping
NX = 10*10
NY = 10*10
cfl = 0.9
tStart = 0.0
tEnd = 4*tSwTransit
nFrames = 20

-- print some diagnostics (seems Jia et. al. define it without the 0.5)
swDynPressure = 0.5*swDensity*(swSpeed[1]^2 + swSpeed[2]^2 + swSpeed[3]^2)
swSoundSpeed = math.sqrt(gasGamma*swPressure/swDensity)
cellSz = (LX[2]-LX[1])/NX/Rm

log(string.format("Solar-wind number density %g", swNumDensity))
log(string.format("Solar-wind temperature %g [K]", swTemp))
log(string.format("Solar-wind pressure %g [Pa]", swPressure))
log(string.format("Solar-wind dynamic pressure %g [Pa]", swDynPressure))
log(string.format("Solar-wind transit across domain %g [s]", tSwTransit))
log(string.format("Solar-wind sound speed %g [m/s]", swSoundSpeed))
log(string.format("Solar-wind Mach number %g [m/s]", math.abs(swSpeed[1]/swSoundSpeed)))
log(string.format("Cell size %g [Rm]", cellSz))


------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = { LX[1], LY[1] },
   upper = { LX[2], LY[2] },
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

-- Function to describe circle
function circle(x, y, x0, y0, rad)
   return (x-x0)^2+(y-y0)^2<rad^2
end

-- in/out field representing embedded object
inOut = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}
inOut:set(
   function (x,y,z)
      local xc, yc = 0.0, 0.0
      local rad = Rm
      return  circle(x, y, 0.0, 0.0, rad) and -1.0 or 1.0
   end
)
inOut:sync()
-- write field
inOut:write("inOut.h5")

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   local rho, u, v, w, pr =  swDensity, swSpeed[1], 0.0, 0.0, swPressure
   local ke = 0.5*rho*(u^2+v^2+w^2)
   return rho, rho*u, rho*v, rho*w, pr/(gasGamma-1)+ke
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- wall BC
bcFluidCopy = BoundaryCondition.Copy { components = {0, 4} }
bcFluidWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }

-- set lower and upper boundaries to walls
bcLowerWallUpdater = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcFluidWall, bcFluidCopy},
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "lower",
}
bcUpperWallUpdater = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcFluidWall, bcFluidCopy},
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "upper",
}

-- updater for embedded BC (solid wall)
embeddedBcUpdater = Updater.StairSteppedBc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcFluidCopy, bcFluidWall},
   -- in/out field
   inOutField = inOut,
}

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt, dir)
   local bcList = {bcLowerWallUpdater, bcUpperWallUpdater}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end

   fld:applyCopyBc(0, "upper")
   fld:applyCopyBc(0, "lower")

   -- apply BCs on embedded boundary
   embeddedBcUpdater:setDir(dir)
   embeddedBcUpdater:setOut( {fld} )
   embeddedBcUpdater:advance(tCurr+myDt)

   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
   useIntermediateWave = true,
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
   updateDirections = {0},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}
-- ds solvers for regular Euler equations along Y
fluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}

-- ds solvers for Lax Euler equations along X
fluidLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}

-- ds solvers for Lax Euler equations along Y
fluidLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false

   -- apply BCs in X before doing X-sweep
   applyBc(q, tCurr, t-tCurr, 0)

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

   -- apply BCs in Y before doing Y-sweep
   applyBc(qX, tCurr, t-tCurr, 1)

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

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)

   -- apply BCs in X before doing X-sweep
   applyBc(q, tCurr, t-tCurr, 0)

   for i,slvr in ipairs({fluidLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
 
   -- apply BCs in Y before doing Y-sweep
   applyBc(qX, tCurr, t-tCurr, 1)

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

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
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
applyBc(q, 0.0, 0.0, 0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)
