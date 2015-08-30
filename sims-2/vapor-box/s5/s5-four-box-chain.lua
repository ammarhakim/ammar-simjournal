-- Gkeyll input file for vapor box simulation. There are no baffles
-- and flows are created by assuming a liquid Lithium wall on all
-- sides of the box.

log = Lucee.logInfo
-- Celcius -> Kelvin
function C2K(Tc) return 273.15+Tc end

-- physical parameters
gasGamma = 5.0/3.0
amu = 1.66053892e-27 -- [kg]
kb = Lucee.BoltzmannConstant -- [J/K]
-- Atomic mass of lithum
mLi = 6.941*amu -- [kg]
-- Atomospheric pressure
atmPress = 101325.0 -- [Pa]

-- initial vapor density and pressure (some arbitrary value). These
-- values are used to initialize a uniform vapor inside the box.
nInit = 4.0e19 -- #/m^3
Tinit = C2K(800) -- [K]
pInit = nInit*kb*Tinit -- [Pa]

-- number of sub-boxes
nBox = 5
boxSize = 0.4 -- [m]
baffleSize = 0.1 -- [m]
holeSize = 0.1 -- [m]

-- vapor box size
Lx = nBox*boxSize+(nBox-1)*baffleSize -- [m]
Ly = 0.4 -- [m]

-- Temperature of each box in [K]. Each box is held to this constant
-- temperature
Tbox = {C2K(950), C2K(787.5), C2K(625), C2K(462.5), C2K(300)}

-- Vapor pressure given local wall temperature [K]. It is valid for
-- liquid Lithium and taken from P. Browning and P.E. Potter, "AN
-- ASSESSMENT OF THE EXPERIMENTALLY DETERMINED VAPOUR PRESSURES OF THE
-- LIQUID ALKALI METALS", in the "Handbook of Thermodynamic and
-- Transport Properties of Alkali Metals", R.W. Ohse, Ed., 1985.
function vaporPressure(Twall)
   return math.exp(26.89-18880/Twall-0.4942*math.log(Twall))
end

function vaporDensity(Twall)
   return vaporPressure(Twall)/(kb*Twall)
end

-- vapor thermal speed
cs0 = math.sqrt(kb*Tinit/mLi) -- [m/2]

-- construct geometry
function inBox(x, y, xlo, xup)
   return (x>xlo[1] and x<xup[1] and y>xlo[2] and y<xup[2]) and 1.0 or -1.0
end
function union2(a, b)
   return math.max(a,b)
end
function union4(a, b, c, d)
   return math.max(a,b,c,d)
end

-- make the nth baffle
function makeBaffle(x,y,n)
   local xl = boxSize*n+baffleSize*(n-1)
   local xu = xl+0.1
   return union2(
      inBox(x,y, {xl,0.00}, {xu,0.15}),
      inBox(x,y, {xl,0.25}, {xu,0.40})
   )
end

-- This function defines the geometry of the baffle. This is used to
-- construct the geometry of domain interior.
function baffleBoxGeometry(x,y)
   return union4(
      makeBaffle(x,y,1), makeBaffle(x,y,2), makeBaffle(x,y,3), makeBaffle(x,y,4)
   )
end

bc = function(x,y,z,t)
   local bs = boxSize+0.5*baffleSize   
   local bd = boxSize+baffleSize
   local Twall = 0
   if x>(bs+3*bd) then
      Twall = Tbox[5]
   elseif x>(bs+2*bd) then
	 Twall = Tbox[4]
   elseif x>(bs+bd) then
      Twall = Tbox[3]
   elseif x>bs then
      Twall = Tbox[2]
   else
      Twall = Tbox[1]
   end
   return Twall
end

-- resolution and time-stepping
NY = 80
NX = NY*nBox
cfl = 0.75
tStart = 0.0
tEnd = 5*Lx/cs0 -- simulation time measure in thermal transits
nFrames = 10

-- print out some initial diagnostics
log (string.format("Initial mass density [kg/m^3]: %g", nInit*mLi))
log (string.format("Initial pressure [Pa]: %g", pInit))
log (string.format("Initial thermal speed of vapor [m/s]: %g", cs0))
log (string.format("Simulation run to [s]: %g", tEnd))
log (string.format("Left box temperature [K]: %g", Tbox[1]))
log (string.format("Right box temperature [K]: %g", Tbox[2]))
log (string.format("Left wall density [#/m^3]: %g", vaporDensity(Tbox[1])))
log (string.format("Right wall density [#/m^3]: %g", vaporDensity(Tbox[2])))
log (string.format("Left wall vapor temperature [K]: %g", vaporDensity(Tbox[1])))
log (string.format("Right wall vapor temperature [K]: %g", vaporDensity(Tbox[2])))
log (string.format("Simulation run to [s]: %g", tEnd))
log ("\n")

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {Lx, Ly},
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
   --writeGhost = {1, 1}
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
-- in/out field representing embedded object
inOut = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}
Twall = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}
Twall:set(bc)
Twall:write("Twall.h5")

-- aliases to various sub-systems
fluid = q:alias(0, 5)
fluidX = qX:alias(0, 5)
fluidNew = qNew:alias(0, 5)

-- construct object defining the geometry
inOut:set(
   function (x,y,z)
      -- negative sign removes the baffle from the outer rectangle to
      -- create the domain geometry
      return -baffleBoxGeometry(x,y)
   end
)
inOut:sync()
-- write field
inOut:write("inOut.h5")

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions (should return rho, rho*u, rho*v, rho*w, E)
function init(x,y,z)
   return nInit*mLi, 0.0, 0.0, 0.0, pInit/(gasGamma-1)
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- bcVaporFunc = BoundaryCondition.FieldFunction {
--    inpComponents = {0, 1, 2},
--    components = {0, 1, 2, 3, 4},
--    bc = function(x,y,z,t, rho, rhou, rhov)
--       local Twall = TwallLeft + x/Lx*(TwallRight-TwallLeft)
--       local pEq = vaporPressure(Twall)
--       local nEq = pEq/(kb*Twall)
--       local rhoEq = nEq*mLi
--       local uin = rhou/rho
--       local vin = rhov/rho
--       return rhoEq, rhoEq*uin, rhoEq*vin, 0.0, pEq/(gasGamma-1)+0.5*rhoEq*(uin^2+vin^2)
--    end,
-- }

bcVaporFunc = BoundaryCondition.Function {
   components = {0, 1, 2, 3, 4},
   bc = function(x,y,z,t)
      local bs = boxSize+0.5*baffleSize
      local bd = boxSize+baffleSize
      local Twall = 0
      if x>(bs+3*bd) then
	 Twall = Tbox[5]
      elseif x>(bs+2*bd) then
	 Twall = Tbox[4]
      elseif x>(bs+bd) then
	 Twall = Tbox[3]
      elseif x>bs then
	 Twall = Tbox[2]
      else
	 Twall = Tbox[1]
      end
      local pEq = vaporPressure(Twall)
      local nEq = pEq/(kb*Twall)
      local rhoEq = nEq*mLi
      local uin = 0.0
      local vin = 0.0
      return rhoEq, rhoEq*uin, rhoEq*vin, 0.0, pEq/(gasGamma-1)+0.5*rhoEq*(uin^2+vin^2)
   end,
}

-- create boundary condition object to apply exact BCs on right/top edges
function createVaporBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = { bcVaporFunc },
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- updater for embedded BC (vapor BCs)
embeddedBcUpdater = Updater.StairSteppedBc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcVaporFunc},
   -- in/out field
   inOutField = inOut,
}

-- create updaters to apply boundary conditions
bcLeft = createVaporBc(0, "lower")
bcRight = createVaporBc(0, "upper")
bcBottom = createVaporBc(1, "lower")
bcTop = createVaporBc(1, "upper")

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   local bcList = {bcLeft, bcRight, bcTop, bcBottom}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end

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
   numericalFlux = "lax", -- on of "lax" or "roe"
   useIntermediateWave = true, -- this only matters is numericalFlux="lax"
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
   updateDirections = {0}, -- directions to update
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
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({fluidLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   applyBc(qX, tCurr, t-tCurr)

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
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- dynvector to store fluid energy
fluidEnergy = DataStruct.DynVector { numComponents = 1 }
fluidEnergyCalc = Updater.IntegrateField2D {
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
   --fluidEnergy:write( string.format("fluidEnergy_%d.h5", frame) )
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
applyBc(q, 0.0, 0.0)
qNew:copy(q)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)


