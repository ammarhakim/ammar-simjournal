-- Program to solve two-fluid equations for ICW problem.

-- global parameters
cfl = 0.48
gasGamma = 5.0/3.0

-- tokamak parameters
Rmaj0 = 3.0 -- [m]
aWall = 1.0 -- [m]
Bmag0 = 4.0 -- [T]

Rmajbgn = 2.83 -- [m]
Zbgn = 0.0 -- [m]

Rmajend = 3.17 -- [m]
Zend = 0.0 -- [m]

cosRmaj =  1.0
cosZ = 0.0

-- number of species and details about them (electron, Deuterium and
-- Helium-3)
nSpecies = 3
chargeNumber = {-1, 1, 2}
massNumber = {0.0005446, 2.0, 3.0}
density0 = {0.400e20, 0.360e20, 0.0}
-- determine Helium 3 density from quasi-neutrality
density0[3] = -(chargeNumber[1]*density0[1] + chargeNumber[2]*density0[2])/chargeNumber[3]
densityWall = {0, 0, 0}

-- some useful quantities
elcCharge = chargeNumber[1]*Lucee.ElementaryCharge
elcMass = massNumber[1]*Lucee.ProtonMass

ionDCharge = chargeNumber[2]*Lucee.ElementaryCharge
ionDMass = massNumber[2]*Lucee.ProtonMass

ionH3Charge = chargeNumber[3]*Lucee.ElementaryCharge
ionH3Mass = massNumber[3]*Lucee.ProtonMass

length = math.sqrt((Rmajend-Rmajbgn)^2+(Zend-Zbgn)^2)
if (cosZ==0.0) then
   xRangeBgn = Rmajbgn
   xRangeEnd = Rmajend
elseif (cosRmaj == 0.0) then
   xRangeBgn = Zbgn
   xRangeEnd = Zend
else
   xRangeBgn = 0.0
   xRangeEnd = length
end

xlower = xRangeBgn
xupper = xRangeEnd
nx = 200

-- compute coordinate of interior last edge
dx = (xupper-xlower)/nx
xLastEdge = xupper-dx
xFirstEdge = xlower+dx

driveF = 40.5e6 -- [Hz]
driveOmega = 2*Lucee.Pi*driveF -- [r/s]

-- computational domain
grid = Grid.RectCart1D {
   lower = {xlower},
   upper = {xupper},
   cells = {nx},
}

-- solution: this has 6 additional components to store the static EM
-- field
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*nSpecies+8+6,
   ghost = {2, 2},
}
q:clear(0.0)

-- updated solution: this has 6 additional components to store the
-- static EM field
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 5*nSpecies+8+6,
   ghost = {2, 2},
}
qNew:clear(0.0)
-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- indices for each subsystem
es,ee, ds,de, h3s,h3e, ems,eme, sebs,sebe = 0,5, 5,10, 10,15, 15,23, 23,29

-- create aliases to various sub-system
elcFluid = q:alias(es, ee) -- electron
ionDFluid = q:alias(ds, de) -- Deuterium
ionH3Fluid = q:alias(h3s, h3e) -- Helium-3
emField = q:alias(ems, eme) -- Em field

elcFluidNew = qNew:alias(es, ee) -- electron
ionDFluidNew = qNew:alias(ds, de) -- Deuterium
ionH3FluidNew = qNew:alias(h3s, h3e) -- Helium-3
emFieldNew = qNew:alias(ems, eme) -- Em field

-- initialize various species
function initSpeciesElc(x,y,z)
   local spIdx = 1 -- electrons
   -- compute location corresponding to x in tokamak
   local Rmajor = Rmajbgn + (x-xRangeBgn)*cosRmaj
   local Z = Zbgn + (x-xRangeBgn)*cosZ
   local Rminor = math.sqrt((Rmajor-Rmaj0)^2 + Z^2)
   local densExp0 = 2.0
   local densExpWall = 1.0
   local fDensProfile = (1-(Rminor/aWall)^densExp0)^densExpWall

   -- now compute number density
   local ns = fDensProfile*(density0[spIdx]-densityWall[spIdx]) + densityWall[spIdx]
   local ts = 0.01*Lucee.Ev2Kelvin -- temperature [K]
   local ps = ns*Lucee.BoltzmannConstant*ts

   return massNumber[spIdx]*Lucee.ProtonMass*ns, 0, 0, 0, ps/(gasGamma-1)
end
elcFluid:set(initSpeciesElc)

function initSpeciesIonD(x,y,z)
   local spIdx = 2 -- Deuterium
   -- compute location corresponding to x in tokamak
   local Rmajor = Rmajbgn + (x-xRangeBgn)*cosRmaj
   local Z = Zbgn + (x-xRangeBgn)*cosZ
   local Rminor = math.sqrt((Rmajor-Rmaj0)^2 + Z^2)
   local densExp0 = 2.0
   local densExpWall = 1.0
   local fDensProfile = (1-(Rminor/aWall)^densExp0)^densExpWall

   -- now compute number density
   local ns = fDensProfile*(density0[spIdx]-densityWall[spIdx]) + densityWall[spIdx]
   local ts = 0.01*Lucee.Ev2Kelvin -- temperature [K]
   local ps = ns*Lucee.BoltzmannConstant*ts

   return massNumber[spIdx]*Lucee.ProtonMass*ns, 0, 0, 0, ps/(gasGamma-1)
end
ionDFluid:set(initSpeciesIonD)

function initSpeciesIonH3(x,y,z)
   local spIdx = 3 -- Helium-3
   -- compute location corresponding to x in tokamak
   local Rmajor = Rmajbgn + (x-xRangeBgn)*cosRmaj
   local Z = Zbgn + (x-xRangeBgn)*cosZ
   local Rminor = math.sqrt((Rmajor-Rmaj0)^2 + Z^2)
   local densExp0 = 2.0
   local densExpWall = 1.0
   local fDensProfile = (1-(Rminor/aWall)^densExp0)^densExpWall

   -- now compute number density
   local ns = fDensProfile*(density0[spIdx]-densityWall[spIdx]) + densityWall[spIdx]
   local ts = 0.01*Lucee.Ev2Kelvin -- temperature [K]
   local ps = ns*Lucee.BoltzmannConstant*ts

   return massNumber[spIdx]*Lucee.ProtonMass*ns, 0, 0, 0, ps/(gasGamma-1)
end
ionH3Fluid:set(initSpeciesIonH3)

-- alias to static magnetic field
staticEB = q:alias(sebs, sebe)
function initStaticEB(x,y,z)
   -- compute location corresponding to x in tokamak
   local Rmajor = Rmajbgn + (x-xRangeBgn)*cosRmaj
   local Z = Zbgn + (x-xRangeBgn)*cosZ
   local Rminor = math.sqrt((Rmajor-Rmaj0)^2 + Z^2)
   local densExp0 = 2.0
   local densExpWall = 1.0
   local fBfieldProfile = Rmaj0/Rmajor

   local Bz = Bmag0*fBfieldProfile

   return 0.0, 0.0, 0.0, 0.1*Bz, 0.0, Bz
end
staticEB:set(initStaticEB)
staticEB:write("staticEB.h5")

-- copy initial conditions over
qNew:copy(q)

-- function to write different fields
function writeFields(dumpNo)
   elcFluid:write( string.format("elc_%d.h5", dumpNo) )
   ionDFluid:write( string.format("D_%d.h5", dumpNo) )
   ionH3Fluid:write( string.format("H3_%d.h5", dumpNo) )
   emField:write( string.format("EM_%d.h5", dumpNo) )
end

-- write initial conditions
writeFields(0)

-- define various equations to solve
elcEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
ionDEulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}
ionH3EulerEqn = HyperEquation.Euler {
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
-- initialize updater
elcEulerSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
elcEulerSlvr:setIn( {elcFluid} )
elcEulerSlvr:setOut( {elcFluidNew} )

-- updater for Deuterium equations
ionDEulerSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionDEulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- initialize updater
ionDEulerSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
ionDEulerSlvr:setIn( {ionDFluid} )
ionDEulerSlvr:setOut( {ionDFluidNew} )

-- updater for Helium-3 equations
ionH3EulerSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionH3EulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- initialize updater
ionH3EulerSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
ionH3EulerSlvr:setIn( {ionH3Fluid} )
ionH3EulerSlvr:setOut( {ionH3FluidNew} )

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
   -- takes density, momentum and EM fields
   inpComponents = {es, es+1, es+2, es+3, ems, ems+1, ems+2, ems+3, ems+4, ems+5},
   -- sets electron momentum and energy source
   outComponents = {es+1, es+2, es+3, es+4},

   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
}
-- Lorentz force on electrons from static magnetic field
elcLorentzForceStaticB = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {es, es+1, es+2, es+3, sebs, sebs+1, sebs+2, sebs+3, sebs+4, sebs+5},
   -- sets momentum and energy source
   outComponents = {es+1, es+2, es+3, es+4},

   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
}
-- electron current contribution to fields
elcCurrent = PointSource.Current {
   -- takes electron momentum
   inpComponents = {es+1, es+2, es+3},
   -- sets current contribution to dE/dt equations
   outComponents = {ems, ems+1, ems+2},

   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
   -- premittivity of free space
   epsilon0 = Lucee.Epsilon0,
}

-- Lorentz force on Deuterium
ionDLorentzForce = PointSource.LorentzForce {
   -- takes density, momentum and EM fields
   inpComponents = {ds, ds+1, ds+2, ds+3, ems, ems+1, ems+2, ems+3, ems+4, ems+5},
   -- sets electron momentum and energy source
   outComponents = {ds+1, ds+2, ds+3, ds+4},

   -- species charge and mass
   charge = ionDCharge,
   mass = ionDMass,
}
-- Lorentz force on Deuterium from static magnetic field
ionDLorentzForceStaticB = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {ds, ds+1, ds+2, ds+3, sebs, sebs+1, sebs+2, sebs+3, sebs+4, sebs+5},
   -- sets momentum and energy source
   outComponents = {ds+1, ds+2, ds+3, ds+4},

   -- species charge and mass
   charge = ionDCharge,
   mass = ionDMass,
}
-- Deuterium current contribution to fields
ionDCurrent = PointSource.Current {
   -- takes electron momentum
   inpComponents = {ds+1, ds+2, ds+3},
   -- sets current contribution to dE/dt equations
   outComponents = {ems, ems+1, ems+2},

   -- species charge and mass
   charge = ionDCharge,
   mass = ionDMass,
   -- premittivity of free space
   epsilon0 = Lucee.Epsilon0,
}

-- Lorentz force on Helium-3
ionH3LorentzForce = PointSource.LorentzForce {
   -- takes density, momentum and EM fields
   inpComponents = {h3s, h3s+1, h3s+2, h3s+3, ems, ems+1, ems+2, ems+3, ems+4, ems+5},
   -- sets electron momentum and energy source
   outComponents = {h3s+1, h3s+2, h3s+3, h3s+4},

   -- species charge and mass
   charge = ionH3Charge,
   mass = ionH3Mass,
}
-- Lorentz force on Helium-3 from static magnetic field
ionH3LorentzForceStaticB = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {h3s, h3s+1, h3s+2, h3s+3, sebs, sebs+1, sebs+2, sebs+3, sebs+4, sebs+5},
   -- sets momentum and energy source
   outComponents = {h3s+1, h3s+2, h3s+3, h3s+4},

   -- species charge and mass
   charge = ionH3Charge,
   mass = ionH3Mass,
}
-- Helium-3 current contribution to fields
ionH3Current = PointSource.Current {
   -- takes electron momentum
   inpComponents = {h3s+1, h3s+2, h3s+3},
   -- sets current contribution to dE/dt equations
   outComponents = {ems, ems+1, ems+2},

   -- species charge and mass
   charge = ionH3Charge,
   mass = ionH3Mass,
   -- premittivity of free space
   epsilon0 = Lucee.Epsilon0,
}

-- Current source from an "antenna"
antennaSrc = PointSource.Function {
   -- contributes to E_y component
   outComponents = {ems+1},
   -- source term to apply
   source = function (x,y,z,t)
	       local J0 = 1.0 -- Amps/m^3
	       local tTurnOn = 10.0/driveF -- source is full magnitude at this time
	       if (x<xFirstEdge) then
		  local ramp = math.sin(0.5*Lucee.Pi*math.min(1, t/tTurnOn))
		  return -J0*ramp^2*math.sin(driveOmega*t)/Lucee.Epsilon0
	       else
		  return 0.0
	       end
	    end,
}

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.GridOdePointIntegrator1D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {elcLorentzForce, elcCurrent, elcLorentzForceStaticB,
	    ionDLorentzForce, ionDCurrent, ionDLorentzForceStaticB,
	    ionH3LorentzForce, ionH3Current, ionH3LorentzForceStaticB,
	    antennaSrc},
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
   -- advance Deuterium
   ionDEulerSlvr:setCurrTime(tCurr)
   ionDStatus, ionDDtSuggested = ionDEulerSlvr:advance(t)
   -- advance Helium-3
   ionH3EulerSlvr:setCurrTime(tCurr)
   ionH3Status, ionH3DtSuggested = ionH3EulerSlvr:advance(t)
   -- advance fields
   maxSlvr:setCurrTime(tCurr)
   maxStatus, maxDtSuggested = maxSlvr:advance(t)

   -- check if any updater failed
   local status, dtSuggested = true, t-tCurr
   if ( (elcStatus == false) or (ionDStatus == false) or (ionH3Status == false) or (maxStatus == false) ) then
      status = false
      dtSuggested = math.min(elcDtSuggested, ionDDtSuggested, ionH3DtSuggested, maxDtSuggested)
   else
      status = true
      dtSuggested = math.min(elcDtSuggested, ionDDtSuggested, ionH3DtSuggested, maxDtSuggested)
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
tEnd = 100/driveF

nFrames = 100
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
   -- main loop
for frame = 1, nFrames do
   print (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested, nanOccured = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   writeFields(frame)
   if (nanOccured) then
      -- no need to continue if nan has occured
      break
   end
   tCurr = tCurr+tFrame
   print ("")
end
