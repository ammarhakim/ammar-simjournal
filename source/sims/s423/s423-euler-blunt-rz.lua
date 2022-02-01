
log = Lucee.logInfo

-- physical parameters
gasGamma = 1.4

Lx = 1.0
Ly = 1.0 --- this is arbitrary (theta direction)
Lz = 2.0

-- inflow conditions
rhoIn = 1.0
prIn = 1.0
csIn = math.sqrt(gasGamma*prIn/rhoIn)
uIn = 2.0*csIn -- supersonic inflow
erIn = 0.5*rhoIn*uIn^2 + prIn/(gasGamma-1)

-- resolution and time-stepping
NX = 100
NY = 1
NZ = 200
cfl = 0.9
tStart = 0.0
tEnd = 5*2.0/uIn
nFrames = 10

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc3D.CartGeneral {}
-- computational domain
grid = Grid.RectCart3D {
   lower = {0.0, 0.0, 0.0},
   upper = {Lx, Ly, Lz},
   cells = {NX, NY, NZ},
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

-- shape library

-- ellipse
Ellipse = { a = 1.0, b = 1.0 }
function Ellipse:isInside(x,y)
   return x^2/self.a^2+y^2/self.b^2 < 1
end
function Ellipse:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

-- Box
Box = { extents = {1.0, 1.0} }
function Box:isInside(x,y)
   return (math.abs(x) < 0.5*self.extents[1]) and (math.abs(y) < 0.5*self.extents[2])
end
function Box:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

-- union adds two shapes together
Union = { a = nil, b = nil }
function Union:isInside(x,y)
   return self.a:isInside(x,y) or self.b:isInside(x,y)
end
function Union:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

-- intersect 
Intersect = { a = nil, b = nil }
function Intersect:isInside(x,y)
   return self.a:isInside(x,y) and self.b:isInside(x,y)
end
function Intersect:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

-- shift object to new origin
Shift = { a = nil, x0 = 0.0, y0 = 0.0 }
function Shift:isInside(x,y)
   return self.a:isInside(x-self.x0, y-self.y0)
end
function Shift:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

shape = Shift:new {
   a = Ellipse:new { a = 0.25, b = 0.5 },
   x0 = 0.0, y0 = Lz/2.0,
}

-- Elipsoid
Ellipsoid = { a = 1.0, b = 1.0, c = 1.0 }
function Ellipsoid:isInside(x,y,z)
   return (x/self.a)^2+(y/self.b)^2+(z/self.c)^2 < 1
end
function Ellipsoid:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

-- in/out field representing embedded object
inOut = DataStruct.Field3D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}
inOut:set(
   function (x,y,z)
      return shape:isInside(x,z) and -1 or 1
   end
)
inOut:sync()
-- write field
inOut:write("inOut.h5")

-----------------------
-- INITIAL CONDITION --
-----------------------

-- initial conditions
function init(r,theta,z)
   return 1e-2*rhoIn, 0.0, 0.0, 0.0, 1e-2*prIn/(gasGamma-1)
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

bcInflow = BoundaryCondition.Const { 
   components = {0, 1, 2, 3, 4},
   values = {rhoIn, 0.0, 0.0, -rhoIn*uIn, erIn},
}
bcInflowUpdater = Updater.Bc3D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcInflow},
   -- direction to apply
   dir = 2,
   -- edge to apply on
   edge = "upper",
}

bcFluidCopy = BoundaryCondition.Copy { components = {0, 4} }
bcFluidWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }

-- updater for embedded BC (solid wall)
embeddedBcUpdater = Updater.StairSteppedBc3D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcFluidCopy, bcFluidWall},
   -- in/out field
   inOutField = inOut,
}

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
function applyBc(fld, tCurr, myDt, dir)
   local bcList = {bcInflowUpdater, bcAxisCalc}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- open BCs on right
   fld:applyCopyBc(0, "upper")
   -- (no need for any BCs in Y direction)
   fld:applyCopyBc(2, "lower")

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
fluidSlvrDir0 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}, -- directions to update
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}
-- ds solvers for regular Euler equations along Z
fluidSlvrDir2 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerEqn,
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {2},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}

-- ds solvers for Lax Euler equations along X
fluidLaxSlvrDir0 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
}
-- ds solvers for Lax Euler equations along Z
fluidLaxSlvrDir2 = Updater.WavePropagation3D {
   onGrid = grid,
   equation = eulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {2},
   hasStairSteppedBoundary = true, -- we are solving with embedded boundary
   inOutField = inOut,
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

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr, 2)

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
   applyBc(q, tCurr, t-tCurr, 0)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, 0)

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


