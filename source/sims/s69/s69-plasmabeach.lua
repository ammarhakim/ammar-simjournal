-- Program to solve two-fluid equations for plasma wave beach
-- problem. The ion fluid is not evolved.

-- global parameters
cfl = 0.4
gasGamma = 5.0/3.0
xlower = 0.0
xupper = 1.0
nx = 400

-- compute coordinate of interior last edge
dx = (xupper-xlower)/nx
xLastEdge = xupper-dx

dx100 = (xupper-xlower)/100
-- compute drive frequency
deltaT = dx100/Lucee.SpeedOfLight
driveOmega = Lucee.Pi/10/deltaT

-- computational domain
grid = Grid.RectCart1D {
   lower = {xlower},
   upper = {xupper},
   cells = {nx},
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- create aliases to various sub-system
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-- function to apply initial conditions
function init(x,y,z)
   local wpdt = 25*(1-x)^5 -- plasma frequency
   local factor = deltaT^2*Lucee.ElementaryCharge^2/(Lucee.ElectronMass*Lucee.Epsilon0)
   local ne = wpdt^2/factor
   local te = 1.0*Lucee.Ev2Kelvin -- electron temperature [K]
   local pre = ne*Lucee.BoltzmannConstant*te
   return Lucee.ElectronMass*ne, 0, 0, 0, pre/(gasGamma-1),
   Lucee.ProtonMass*ne, 0, 0, 0, pre/(gasGamma-1), 
   0, 0, 0, 0, 0, 0, 0, 0
end
-- set initial conditions
q:set(init)

-- copy initial conditions over
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- define various equations to solve
elcEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}

maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = Lucee.SpeedOfLight,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 0.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 0.0,
}

-- updater for electron equations
elcEulerSlvr = Updater.WavePropagation1D {
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

-- updater for Maxwell equations
maxSlvr = Updater.WavePropagation1D {
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

-- Lorentz force on electrons
elcLorentzForce = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {0, 1, 2, 3, 10, 11, 12, 13, 14, 15},
   -- sets electron momentum and energy source
   outComponents = {1, 2, 3, 4},

   -- species charge and mass
   charge = -Lucee.ElementaryCharge,
   mass = Lucee.ElectronMass,
}

-- electron current contribution to fields
elcCurrent = PointSource.Current {
   -- takes electron momentum
   inpComponents = {1, 2, 3},
   -- sets current contribution to dE/dt equations
   outComponents = {10, 11, 12},

   -- species charge and mass
   charge = -Lucee.ElementaryCharge,
   mass = Lucee.ElectronMass,
   -- premittivity of free space
   epsilon0 = Lucee.Epsilon0,
}

-- Current source from an "antenna"
antennaSrc = PointSource.Function {
   -- contributes to E_y component
   outComponents = {11},
   -- source term to apply
   source = function (x,y,z,t)
	       local J0 = 1.0e-12 -- Amps/m^3
	       if (x>xLastEdge) then
		  return -J0*math.sin(driveOmega*t)/Lucee.Epsilon0
	       else
		  return 0.0
	       end
	    end,
}

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.GridOdePointIntegrator1D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {elcLorentzForce, elcCurrent, antennaSrc},
}
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
   if ( (elcStatus == false) or (maxStatus == false) ) then
      status = false
      dtSuggested = math.min(elcDtSuggested, maxDtSuggested)
   else
      status = true
      dtSuggested = math.min(elcDtSuggested, maxDtSuggested)
   end

   -- update source terms
   if (status) then
      sourceSlvr:setCurrTime(tCurr)
      sourceSlvr:advance(t)
   end

   return status, dtSuggested
end

-- advance solution from tStart to tEnd, using optimal time-steps.
function advanceFrame(tStart, tEnd, initDt)

   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local nanOccured = false
   while true do
      -- copy qNew in case we need to take this step again
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format(" Taking step %d at time %g with dt %g", step, tCurr, myDt))
      -- advance fluids and fields
      status, dtSuggested = solveTwoFluidSystem(tCurr, tCurr+myDt)

      if (dtSuggested < myDt) then
	 -- time-step too large
	 print (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
      else
	 -- apply outflow BCs
	 qNew:applyCopyBc(0, "lower")
	 qNew:applyCopyBc(0, "upper")

	 -- check if a nan occured
	 if (qNew:hasNan()) then
	    print (string.format(" ** Nan occured at %g! Stopping simulation", tCurr))
	    nanOccured = true
	    break
	 end

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
   
   return dtSuggested, nanOccured
end

dtSuggested = 100.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 5.0e-9

nFrames = 100
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
   -- main loop
for frame = 1, nFrames do
   print (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested, nanOccured = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   qNew:write( string.format("q_%d.h5", frame) )
   if (nanOccured) then
      -- no need to continue if nan has occured
      break
   end
   tCurr = tCurr+tFrame
   print ("")
end
