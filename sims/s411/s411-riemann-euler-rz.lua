
log = Lucee.logInfo

-- physical parameters
gasGamma = 1.4

Lx = 1.5
Ly = 100.0 --- this is arbitrary (theta direction)
Lz = 1.0

-- resolution and time-stepping
NX = 150
NY = 1
NZ = 100
cfl = 0.9
tStart = 0.0
tEnd = 0.7
nFrames = 2

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc3D.CartGeneral {}
-- computational domain
grid = Grid.RectCart3D {
   lower = {0.0, 0.0, 0.0},
   upper = {Lx, Ly, Lz},
   cells = {NX, NY, NX},
   decomposition = decomp,
   periodicDirs = {},
}

-- solution
q = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 5,
   ghost = {2, 2},
}

-- aliases to various sub-systems
fluid = q:alias(0, 5)
fluidX = qX:alias(0, 5)
fluidNew = qNew:alias(0, 5)

-----------------------
-- INITIAL CONDITION --
-----------------------

-- initial conditions
function init(r,theta,z)
   -- See Langseth and LeVeque, section 3.2

   local rhoi, pri = 1.0, 5.0
   local rho0, pr0 = 1.0, 1.0
   local rho, pr

   local rloc = math.sqrt(r^2+(z-0.4)^2)
   if (rloc<0.2) then
      rho, pr = rhoi, pri
   else
      rho, pr = rho0, pr0
   end

   return rho, 0.0, 0.0, 0.0, pr/(gasGamma-1)
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

bcFluidCopy = BoundaryCondition.Copy { components = {0, 4} }
bcFluidWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }

-- create boundary condition object to apply wall BCs
function createWallBc(myDir, myEdge)
   local bc = Updater.Bc3D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {
	 bcFluidCopy, bcFluidWall,
      },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- walls
bc_Z0 = createWallBc(2, "lower")
bc_Z1 = createWallBc(2, "upper")

-- axis BCs
bcAxis = BoundaryCondition.Copy { 
   components = {0, 1, 2, 3, 4},
   fact = {1, -1, -1, 1, 1},
}

bcAxisCalc = Updater.Bc3D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = { bcAxis },
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "lower",
}

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   local bcList = {bc_Z0, bc_Z1, bcAxisCalc}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- open BCs on right
   fld:applyCopyBc(0, "upper")
   -- (no need for any BCs in Y direction)

   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
eulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}

-- ds solvers for regular Euler equations along X
fluidSlvrDir0 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
-- ds solvers for regular Euler equations along Z
fluidSlvrDir2 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {2}
}

-- ds solvers for Lax Euler equations along X
fluidLaxSlvrDir0 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
-- ds solvers for Lax Euler equations along Z
fluidLaxSlvrDir2 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {2}
}

-- gravitational source
axisSrc = PointSource.EulerAxisymmetric {
   -- takes and returns fluid variables
   inpComponents = {0, 1, 2, 3, 4},
   outComponents = {0, 1, 2, 3, 4},
   gasGamma = gasGamma,
}
-- updater to add gravitational force to fluid
axisSrcSlvr = Updater.GridOdePointIntegrator3D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {axisSrc},
}

-- function to update source terms
function updateSource(qIn, tCurr, t)
   -- gravity source
   axisSrcSlvr:setOut( {qIn} )
   axisSrcSlvr:setCurrTime(tCurr)
   axisSrcSlvr:advance(t)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false

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

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr)

   -- for axisymmetric problems, there is no Y update

   -- Z-direction updates
   for i,slvr in ipairs({fluidSlvrDir2}) do
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

   -- update source terms
   updateSource(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   -- X-direction updates
   for i,slvr in ipairs({fluidLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   applyBc(qX, tCurr, t-tCurr)

   -- for axisymmetric problems, there is no Y update

   -- Z-direction updates
   for i,slvr in ipairs({fluidLaxSlvrDir2}) do
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
   updateSource(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- dynvector to store fluid energy
fluidEnergy = DataStruct.DynVector { numComponents = 1 }
fluidEnergyCalc = Updater.IntegrateField3D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (rho, rhou, rhov, rhow, er)
		  return er
	       end,
}
fluidEnergyCalc:setIn( {fluid} )
fluidEnergyCalc:setOut( {fluidEnergy} )

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({fluidEnergyCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )
   fluidEnergy:write( string.format("fluidEnergy_%d.h5", frame) )
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
        qNew:copy(qNewDup)
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

-- set input/output arrays for various solvers

-- Regular Euler solvers
fluidSlvrDir0:setIn( {fluid} )
fluidSlvrDir0:setOut( {fluidX} )

fluidSlvrDir2:setIn( {fluidX} )
fluidSlvrDir2:setOut( {fluidNew} )

-- Lax Euler solvers
fluidLaxSlvrDir0:setIn( {fluid} )
fluidLaxSlvrDir0:setOut( {fluidX} )

fluidLaxSlvrDir2:setIn( {fluidX} )
fluidLaxSlvrDir2:setOut( {fluidNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 1.0
runSimulation(tStart, tEnd, nFrames, initDt)


