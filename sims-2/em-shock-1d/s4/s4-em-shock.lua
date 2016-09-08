-- Vlasov-Maxwell solver 

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- permiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

Te_Ti = 9.0 -- ratio of electron to ion temperature
machNum = 1.5 -- Mach number computed from ion thermal speed
n0 = 1.0 -- initial number density
elcTemp = 1.0e-2 -- electron temperature
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge
beta = 2.0 -- plasma beta

ionTemp = elcTemp/Te_Ti -- ion temperature
ionMass = 1836.2 -- ion mass
ionCharge = 1.0 -- ion charge

-- magnetic field
B0 = math.sqrt(2*mu0*n0*(elcTemp+ionTemp)/beta)
BzL = B0 -- magnetic field on left
BzR = B0 -- magnetic field on right

-- thermal speeds
cs = math.sqrt(elcTemp/ionMass)
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
wpi = math.sqrt(ionCharge^2*n0/(epsilon0*ionMass))
lambdaD = vtElc/wpe

-- electron and ion drift speeds
elcDrift = machNum*cs
ionDrift = elcDrift -- no net current

-- domain size and simulation time
LX = 20*lightSpeed/wpe
tStart = 0.0 -- start time 
tEnd = 500.0/wpe
nFrames = 50
nFramesDistf = 50 --write out distribution function with less cadence due to size

-- Resolution, time-stepping etc.
NX = 64
NVX = 16
NVY = 16

cfl = 0.5/(2*polyOrder+1)

-- compute max thermal speed to set velocity space extents
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc
VL_ION, VU_ION = -6.0*vtIon-ionDrift, 6.0*vtIon+ionDrift

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Sound speed=%g", cs))
log(string.format("Mach number=%g", machNum))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Plasma frequency=%g", wpe))
log(string.format("Debye length=%g", lambdaD))
log(string.format("Ion inertial lenght=%g", lightSpeed/wpe*math.sqrt(ionMass/elcMass)))
log(string.format("Electron inertial lenght=%g", lightSpeed/wpe))
log(string.format("Magnetic field=%g", B0))
log(string.format("Domain size=%g", LX))
log(string.format("Cell size=%g", LX/NX))
log(string.format("Number of Debye Lengths per grid point=%g", (LX/lambdaD)/NX))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Electron/Ion drift speed=%g", elcDrift))
log(string.format("Electron domain extents = [%g,%g]", VL_ELC, VU_ELC))
log(string.format("Ion domain extents = [%g,%g]", VL_ION, VU_ION))
log(string.format("Edge of velocity space/speed of light = %g", VU_ELC))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc3D.CartProd { cuts = {4,2,2} }
confDecomp = DecompRegionCalc1D.SubCartProd3D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid for electrons
phaseGridElc = Grid.RectCart3D {
   lower = {0.0, VL_ELC, VL_ELC},
   upper = {LX, VU_ELC, VU_ELC},
   cells = {NX, NVX, NVY},
   decomposition = phaseDecomp,   
}
-- phase space grid for ions
phaseGridIon = Grid.RectCart3D {
   lower = {0.0, VL_ION, VL_ION},
   upper = {LX, VU_ION, VU_ION},
   cells = {NX, NVX, NVY},
   decomposition = phaseDecomp,   
}

-- configuration space grid (same for electrons and ions)
confGrid = Grid.RectCart1D {
   lower = {0.0},
   upper = {LX},
   cells = {NX},
   decomposition = confDecomp,
}

-- phase-space basis functions for electrons and ions
phaseBasisElc = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridElc,
   polyOrder = polyOrder,
}
phaseBasisIon = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridIon,
   polyOrder = polyOrder,
}
-- configuration-space basis functions (shared by both species)
confBasis = NodalFiniteElement1D.SerendipityElement {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- extra fields for performing RK update
distfNewElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}

-- Electron number density
numDensityElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion number density
numDensityIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron momentum
momentumElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
--Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
--Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
--Electron pressure tensor
pressureTensorElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 3*confBasis:numNodes(),
   ghost = {1, 1},
}

--Ion pressure tensor
pressureTensorIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 3*confBasis:numNodes(),
   ghost = {1, 1},
}
-- net current
current = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for adding to EM fields (this perhaps is not the best way to do it)
emSource = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}

-- EM field
em = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
em1 = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
emNew = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}

--------------------------------
-- INITIAL CONDITION UPDATERS --
--------------------------------

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdriftx, vdrifty, vt, vx, vy)
   local v2 = (vx - vdriftx)^2 + (vy - vdrifty)^2
   return n0/(2*Lucee.Pi*vt^2)*math.exp(-v2/(2*vt^2))
end

-- updater to initialize electron distribution function
initDistfElc = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,vx,vy,t)
      local sloc = 0.5*LX
      local fxv = maxwellian(n0, elcDrift, 0.0, vtElc, vx, vy)
      if x>sloc then
	 fxv = maxwellian(n0, -elcDrift, 0.0, vtElc, vx, vy)
      end
      return fxv
   end
}

-- updater to initialize ion distribution function
initDistfIon = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,vx,vy,t)
      local sloc = 0.5*LX
      local fxv = maxwellian(n0, ionDrift, 0.0, vtIon, vx, vy)
      if x>sloc then
	 fxv = maxwellian(n0, -ionDrift, 0.0, vtIon, vx, vy)
      end
      return fxv
   end
}

-- field initial condition to apply
function initEMFields(x,y,z)
   local sloc = 0.5*LX
   local Bz = BzL
   if x>sloc then
      Bz = BzR
   end
   return 0.0, 0.0, 0.0, 0.0, 0.0, Bz, 0.0, 0.0
end

-- updater to apply initial conditions for fields
initField = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      return initEMFields(x,y,z)
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   cfl = cfl,
   charge = elcCharge,
   mass = elcMass,
   polyOrder = polyOrder,
}
vlasovPositivityElc = Updater.NodalDgScalingLimiter1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
}
vlasovSolverIon = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   cfl = cfl,
   charge = ionCharge,
   mass = ionMass,
   polyOrder = polyOrder,
}
vlasovPositivityIon = Updater.NodalDgScalingLimiter1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
}

-- Maxwell equation object
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = lightSpeed,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 0.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 1.0,
   -- numerical flux to use: one of "upwind" or "central"
   numericalFlux = "upwind",
}

-- updater to solve Maxwell equations
maxwellSlvr = Updater.NodalDgHyper1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- equation system to solver
   equation = maxwellEqn,
   -- CFL number
   cfl = cfl,
}

-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 0, -- moment to compute
}
numDensityCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 0, --moment to compute
}

-- Updater to compute momentum
momentumCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 1,
}

-- Updater to compute scalar energy
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 2,
   scalarPtclEnergy = true,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 2,
   scalarPtclEnergy = true,
}

-- Updater to compute independent components of pressure tensor
pressureTensorCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 2,
}
pressureTensorCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 2,
}

-- This strange looking updater copies the currents into the EM source
-- field. Perhaps this is not the best way to do things, and one can
-- imagine a source updater which adds current sources to the dE/dt
-- Maxwell equation
copyToEmSource = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {0,1,2},
   targetComponents = {0,1,2},
}

-------------------------
-- Boundary Conditions --
-------------------------

-- boundary applicator objects for fluids and fields
distfElcBc = BoundaryCondition.NodalDgCopy3D {
   components = {0},
   fact = {1},
   basis = phaseBasisElc,
}

-- create boundary condition object
function createBcElc(myDir, myEdge)
   local bc = Updater.Bc3D {
      onGrid = phaseGridElc,
      -- boundary conditions to apply
      boundaryConditions = {distfElcBc},
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   bc:setOut( {fldElc} )
   return bc
end

-- boundary applicator objects for fluids and fields
distfIonBc = BoundaryCondition.NodalDgCopy3D {
   components = {0},
   fact = {1},
   basis = phaseBasisIon,
}

-- create boundary condition object
function createBcIon(myDir, myEdge)
   local bc = Updater.Bc3D {
      onGrid = phaseGridIon,
      -- boundary conditions to apply
      boundaryConditions = {distfIonBc},
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   bc:setOut( {fldIon} )
   return bc
end

BcXLeftElc = createBcElc(0, "lower")
BcXRightElc = createBcElc(0, "upper")

BcXLeftIon = createBcIon(0, "lower")
BcXRightIon = createBcIon(0, "upper")

function applyDistFuncBc(tCurr, myDt, fldElc, fldIon)

   local bcListElc = {BcXLeftElc, BcXRightElc}
   for i,bc in ipairs(bcListElc) do
      bc:setOut( {fldElc} )
      bc:advance(tCurr+myDt)
   end
   local bcListIon = {BcXLeftIon, BcXRightIon}
   for i,bc in ipairs(bcListIon) do
      bc:setOut( {fldIon} )
      bc:advance(tCurr+myDt)
   end

   -- sync the distribution function across processors
   fldElc:sync()
   fldIon:sync()

end

-- boundary applicator objects for fluids and fields
emBc = BoundaryCondition.NodalDgCopy1D {
   components = {0,1,2,3,4,5,6,7},
   fact = {1,1,1,1,1,1,1,1},
   basis = confBasis,
}

-- create boundary condition object
function createBcEm(myDir, myEdge)
   local bc = Updater.Bc1D {
      onGrid = confGrid,
      -- boundary conditions to apply
      boundaryConditions = {emBc},
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

BcXLeftEM = createBcEm(0, "lower")
BcXRightEM = createBcEm(0, "upper")

-- apply BCs to EM fields
function applyEmBc(tCurr, myDt, EM)
   for i,bc in ipairs({BcXLeftEM, BcXRightEM}) do
      --runUpdater(bc, tCurr, myDt, {}, {EM})
   end
   EM:applyCopyBc(0, "lower")
   EM:applyCopyBc(0, "upper")
   -- sync EM fields across processors
   EM:sync()
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

----------------------
-- SOLVER UTILITIES --
----------------------


-- generic function to run an updater
function runUpdater(updater, tCurr, myDt, inpFlds, outFlds)
   updater:setCurrTime(tCurr)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(tCurr+myDt)
end


function calcNumDensity(calculator, tCurr, myDt, distfIn, numDensOut)
   return runUpdater(calculator, tCurr, myDt, {distfIn}, {numDensOut})
end
-- function to calculate momentum density
function calcMomentum(calculator, tCurr, myDt, distfIn, momentumOut)
   return runUpdater(calculator, tCurr, myDt, {distfIn}, {momentumOut})
end
-- function to calculate energy density
function calcPtclEnergy(calculator, tCurr, myDt, distfIn, energyOut)
  return runUpdater(calculator, tCurr, myDt, {distfIn}, {energyOut})
end
-- function to calculate pressure tensor
function calcPressureTensor(calculator, tCurr, myDt, distfIn, pressureTensorOut)
  return runUpdater(calculator, tCurr, myDt, {distfIn}, {pressureTensorOut})
end

-- function to compute moments from distribution function
function calcMoments(tCurr, myDt, distfElcIn, distfIonIn)
   -- number density
   calcNumDensity(numDensityCalcElc, tCurr, myDt, distfElcIn, numDensityElc)
   calcNumDensity(numDensityCalcIon, tCurr, myDt, distfIonIn, numDensityIon)
   -- momemtum
   calcMomentum(momentumCalcElc, tCurr, myDt, distfElcIn, momentumElc)
   calcMomentum(momentumCalcIon, tCurr, myDt, distfIonIn, momentumIon)
   -- scalar energy
   calcPtclEnergy(ptclEnergyCalcElc, tCurr, myDt, distfElcIn, ptclEnergyElc)
   calcPtclEnergy(ptclEnergyCalcIon, tCurr, myDt, distfIonIn, ptclEnergyIon)
   -- pressure tensor
   calcPressureTensor(pressureTensorCalcElc, tCurr, myDt, distfElcIn, pressureTensorElc)
   calcPressureTensor(pressureTensorCalcIon, tCurr, myDt, distfIonIn, pressureTensorIon)
end

-- function to update Vlasov equation
function updateVlasovEqn(vlasovSlvr, tCurr, myDt, distfIn, emIn, distfOut)
   return runUpdater(vlasovSlvr, tCurr, myDt, {distfIn, emIn}, {distfOut})
end

function enforcePositivity(positivityFix, tCurr, myDt, distfOut)
   return 0
   --return runUpdater(positivityFix, tCurr, myDt, {}, {distfOut})
end

-- solve maxwell equation
function updateMaxwellEqn(tCurr, myDt, emIn, emOut)
   maxwellSlvr:setCurrTime(tCurr)
   maxwellSlvr:setIn( {emIn} )
   maxwellSlvr:setOut( {emOut} )
   return maxwellSlvr:advance(tCurr+myDt)
end

-- function to compute diagnostics
function calcDiagnostics(tCurr, myDt)
   -- Nothing to do
end

----------------------------
-- Time-stepping routines --
----------------------------

-- take single RK step
function rkStage(tCurr, myDt, elcIn, ionIn, emIn, elcOut, ionOut, emOut)
   -- update distribution functions and homogenous Maxwell equations
   local stElc, dtElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, elcIn, emIn, elcOut)
   enforcePositivity(vlasovPositivityElc, tCurr, myDt, elcOut)
   local stIon, dtIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, ionIn, emIn, ionOut)
   enforcePositivity(vlasovPositivityIon, tCurr, myDt, ionOut)
   local stEm, dtEm = updateMaxwellEqn(tCurr, myDt, emIn, emOut)
   
   -- compute currents from old values
   calcMomentum(momentumCalcElc, tCurr, myDt, elcIn, momentumElc)
   calcMomentum(momentumCalcIon, tCurr, myDt, ionIn, momentumIon)
   current:combine(elcCharge, momentumElc, ionCharge, momentumIon)
   -- copy into EM sources
   runUpdater(copyToEmSource, tCurr, myDt, {current}, {emSource})
   -- add in current source to Maxwell equation output
   emOut:accumulate(-myDt/epsilon0, emSource)
   
   if (stElc == false) or (stIon == false) or (stEm == false)  then
      return false, math.min(dtElc, dtIon, dtEm)
   end
   return true, math.min(dtElc, dtIon, dtEm)
end

function rk3(tCurr, myDt)
   local myStatus, myDtSuggested
   -- RK stage 1
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distfElc, distfIon, em, distf1Elc, distf1Ion, em1)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   -- apply BC
   applyDistFuncBc(tCurr, myDt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, myDt, em1)

   -- RK stage 2
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Elc, distf1Ion, em1, distfNewElc, distfNewIon, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   em1:combine(3.0/4.0, em, 1.0/4.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, myDt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, myDt, em1)   

   -- RK stage 3
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Elc, distf1Ion, em1, distfNewElc, distfNewIon, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   em1:combine(1.0/3.0, em, 2.0/3.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, myDt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, myDt, em1)   

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)
   em:copy(em1)

   return true, myDtSuggested
end

-- make duplicates in case we need them
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()
emDup = em:duplicate()

-- function to advance solution from tStart to tEnd
function runSimulation(tStart, tEnd, nFrames, initDt)
   local frame = 1
   local frameDistf = 1
   local tFrame = (tEnd-tStart)/nFrames
   local tFrameDistf = (tEnd-tStart)/nFramesDistf
   local nextIOt = tFrame
   local nextIOtDistf = tFrameDistf
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- the grand loop 
   while true do
      distfDupElc:copy(distfElc)
      distfDupIon:copy(distfIon)
      emDup:copy(em)
      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      -- advance particles and fields
      log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)
      if (status == false) then
	 -- time-step too large
	 log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 distfElc:copy(distfDupElc)
	 distfIon:copy(distfDupIon)
	 em:copy(emDup)
      else
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)

         -- write out data
         if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then

            log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
            writeFields(frame, tCurr+myDt)
            frame = frame + 1
            nextIOt = nextIOt + tFrame
            step = 0

            if(tCurr + myDt > nextIOtDistf or tCurr+myDt >= tEnd) then

               log (string.format(" Writing distribution function at time %g (frame %d) ...\n", tCurr+myDt, frameDistf))
               writeDistributionFunction(frameDistf, tCurr+myDt)
               frameDistf = frameDistf + 1
               nextIOtDistf = nextIOtDistf + tFrameDistf

            end
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

-- Write out data frame 'frameNum' with at specified time 'tCurr'
function writeFields(frameNum, tCurr)

   -- EM field
   em:write(string.format("em_%d.h5", frameNum), tCurr)

   -- compute moments and write them out
   calcMoments(tCurr, 0.0, distfElc, distfIon)
   numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)
   momentumElc:write(string.format("momentumElc_%d.h5", frameNum), tCurr)
   momentumIon:write(string.format("momentumIon_%d.h5", frameNum), tCurr)
   ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)
   ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr)
   pressureTensorElc:write(string.format("pressureTensorElc_%d.h5", frameNum), tCurr)
   pressureTensorIon:write(string.format("pressureTensorIon_%d.h5", frameNum), tCurr)

end
-- Write out distribution function frame 'frameNum' with at specified time 'tCurr'
function writeDistributionFunction(frameNum, tCurr)
   -- distribution functions
   distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
   distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- apply initial conditions for electrons and ion
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
-- apply initial conditions for electrons and ion
runUpdater(initField, 0.0, 0.0, {}, {em})

-- apply BCs
applyDistFuncBc(0.0, 0.0, distfElc, distfIon)
applyEmBc(0.0, 0.0, em)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial fields
writeFields(0, 0.0)
writeDistributionFunction(0, 0.0)

-- run the whole thing
initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)
