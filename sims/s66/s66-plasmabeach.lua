-- Program to solve Two-Fluid equations

-- simulation parameters
cfl = 0.9

-- global settings
xlower = 0.0
xupper = 1.0
nx = 100

-- compute coordinate of interior last edge
dx = (xupper-xlower)/nx
xLastEdge = xupper-dx

-- compute drive frequency
deltaT = dx/Lucee.SpeedOfLight
driveOmega = Lucee.Pi/10/deltaT

-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {1000},
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 16,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- create aliases to various sub-system
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 16)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 16)

-- initial conditions to apply
function initElc(x,y,z)
   local sloc = 0.5 -- location of shock
   if (x<sloc) then
      return 1.0*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
   else
      return 0.125*elcMass/ionMass, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
   end
end
elcFluid:set(initElc)

function initIon(x,y,z)
   local sloc = 0.5 -- location of shock
   if (x<sloc) then
      return 1.0, 0.0, 0.0, 0.0, 5.0e-5/(gasGamma-1)
   else
      return 0.125, 0.0, 0.0, 0.0, 5.0e-6/(gasGamma-1)
   end
end
ionFluid:set(initIon)

function initEm(x,y,z)
   local sloc = 0.5 -- location of shock
   if (x<sloc) then
      return 0.0, 0.0, 0.0, 0.75e-2, 0.0, 1.0e-2
   else
      return 0.0, 0.0, 0.0, 0.75e-2, 0.0, -1.0e-2
   end
end
emField:set(initEm)

-- copy initial conditions over
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- define various equations to solve
elcEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
ionEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
maxwellEqn = HyperEquation.Maxwell {
   -- speed of light
   lightSpeed = lightSpeed
}

-- updater for electron equations
elcEulerSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- initialize updater
elcEulerSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
elcEulerSlvr:setIn( {elcFluid} )
elcEulerSlvr:setOut( {elcFluidNew} )

-- updater for Maxwell equations
maxSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- initialize updater
maxSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
maxSlvr:setIn( {emField} )
maxSlvr:setOut( {emFieldNew} )

-- Lorentz force on electrons
elcLorentzForce = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {0, 1, 2, 3, 10, 11, 12, 13, 14, 15},
   -- sets electron momentum and energy source
   outComponents = {1, 2, 3, 4},

   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
}

-- electron current contribution to fields
elcCurrent = PointSource.Current {
   -- takes electron momentum
   inpComponents = {1, 2, 3},
   -- sets current contribution to dE/dt equations
   outComponents = {10, 11, 12},

   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
   -- premittivity of free space
   epsilon0 = epsilon0,
}

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.GridOdePointIntegrator1D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {elcLorentzForce, elcCurrent},
}
-- initialize updater
sourceSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
sourceSlvr:setOut( {qNew} )

-- function to take one time-step
function solveTwoFluidSystem(tCurr, t)
   -- advance electrons
   elcEulerSlvr:setCurrTime(tCurr)
   elcStatus, elcDtSuggested = elcEulerSlvr:advance(t)
   -- advance fields
   maxSlvr:setCurrTime(tCurr)
   maxStatus, maxDtSuggested = maxSlvr:advance(t)

   -- check if any updater failed
   local status, dtSuggested = true, t-tCurr
   if ( (elcStatus == false) or (ionStatus == false) or (maxStatus == false) ) then
      status = false
      dtSuggested = math.min(elcDtSuggested, ionDtSuggested, maxDtSuggested)
   else
      status = true
      dtSuggested = math.min(elcDtSuggested, ionDtSuggested, maxDtSuggested)
   end

   -- update source terms
   if (status) then
      sourceSlvr:setCurrTime(tCurr)
      sourceSlvr:advance(t)
   end

   return status, dtSuggested
end

myDt = 100.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 5e-9

tCurr = tStart
step = 1

-- main loop
while true do
   -- copy qNew in case we need to take this step again
   qNewDup:copy(qNew)

   -- if needed adjust dt to hit tEnd exactly
   if (tCurr+myDt > tEnd) then
      myDt = tEnd-tCurr
   end

   print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
   -- advance fluids and fields
   status, dtSuggested = solveTwoFluidSystem(tCurr, tCurr+myDt)

   if (dtSuggested < myDt) then
      -- time-step too large
      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      myDt = dtSuggested
      qNew:copy(qNewDup)
   else
      -- apply copy BCs on lower and upper edges
      qNew:applyCopyBc(0, "lower")
      qNew:applyCopyBc(0, "upper")

      -- copy updated solution back
      q:copy(qNew)

      tCurr = tCurr + myDt
      step = step + 1
      -- check if done
      if (tCurr >= tEnd) then
	 break
      end
   end
end

-- write final solution
q:write("q_1.h5")