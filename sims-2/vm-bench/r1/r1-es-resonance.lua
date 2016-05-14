-- Vlasov-Poisson solver in Poisson-Bracket formulation

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
epsilon0 = 1.0 -- permittivity of free space

Te_Ti = 1.0 -- ratio of electron to ion temperaute
n0 = 1.0 -- initial number density
elcTemp = 1.0 -- electron temperature
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

ionTemp = elcTemp/Te_Ti -- ion temperature
ionMass = 1836.2 -- ion mass
ionCharge = 1.0 -- ion charge

-- thermal speeds
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- wave-number
knumber = 0.5/lambdaD
-- field amplitude
phi0 = 0.01
-- velocity at which we want resonance
vRes = 1.0*vtElc
-- field frequency
omega0 = knumber*vRes
-- domain size and simulation time
LX = 2*Lucee.Pi/knumber
tStart = 0.0 -- start time 
tEnd = 20.0/wpe
nFrames = 10

-- Resolution, time-stepping etc.
NX = 64
NV = 32
polyOrder = 2

cfl = 0.5/(2*polyOrder+1)

-- compute max thermal speed to set velocity space extents
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc
VL_ION, VU_ION = -10.0*vtIon, 10.0*vtIon

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Plasma frequency=%g", wpe))
log(string.format("Drive frequency=%g", omega0))
log(string.format("Resonance velocity=%g", vRes))
log(string.format("Debye length=%g", lambdaD))
log(string.format("Cell size=%g", LX/NX))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Electron domain extents = [%g,%g]", VL_ELC, VU_ELC))
log(string.format("Ion domain extents = [%g,%g]", VL_ION, VU_ION))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc2D.CartProd { cuts = {2,1} }
confDecomp = DecompRegionCalc1D.SubCartProd2D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid for electrons
phaseGridElc = Grid.RectCart2D {
   lower = {-LX/2, VL_ELC},
   upper = {LX/2, VU_ELC},
   cells = {NX, NV},
   decomposition = phaseDecomp,
   periodicDirs = {0},
}
-- phase space grid for ions
phaseGridIon = Grid.RectCart2D {
   lower = {-LX/2, VL_ION},
   upper = {LX/2, VU_ION},
   cells = {NX, NV},
   decomposition = phaseDecomp,
   periodicDirs = {0},
}

-- configuration space grid (same for electrons and ions)
confGrid = Grid.RectCart1D {
   lower = {-LX/2},
   upper = {LX/2},
   cells = {NX},
   decomposition = confDecomp,
   periodicDirs = {0},      
}

-- phase-space basis functions for electrons and ions
phaseBasisElc = NodalFiniteElement2D.Serendipity {
   onGrid = phaseGridElc,
   polyOrder = polyOrder,
}
phaseBasisIon = NodalFiniteElement2D.Serendipity {
   onGrid = phaseGridIon,
   polyOrder = polyOrder,
}
-- configuration-space basis functions (shared by both species)
confBasis = NodalFiniteElement1D.Lobatto {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field2D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- extra fields for performing RK update
distfNewElc = DataStruct.Field2D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field2D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
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
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- charge density
chargeDensity = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- field to store potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- drive potential
phiDrive = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field2D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- create field to store kinetic energy term in Hamiltonian
hamilKeElc = DataStruct.Field2D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
hamilKeIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

--------------------------------
-- INITIAL CONDITION UPDATERS --
--------------------------------

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- updater to initialize distribution function
initDistfElc = Updater.ProjectOnNodalBasis2D {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
      return maxwellian(n0, 0.0, vtElc, v)
   end
}

-- updater to initialize distribution function
initDistfIon = Updater.ProjectOnNodalBasis2D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
      return maxwellian(n0, 0.0, vtIon, v)
   end
}

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilKeElc = Updater.EvalOnNodes2D {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilKeIon = Updater.EvalOnNodes2D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}

-- updater to compute drive potential
calcDrivePhi = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
      return  phi0*math.cos(knumber*x-omega0*t)
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.PoissonBracket {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   hamilNodesShared = false, -- Hamiltonian is not continuous
   zeroFluxDirections = {1},
}
vlasovSolverIon = Updater.PoissonBracket {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   hamilNodesShared = false, -- Hamiltonian is not continuous
   zeroFluxDirections = {1},
}

-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridElc,
   basis2d = phaseBasisElc,
   basis1d = confBasis,
   moment = 0, -- moment to compute
}
numDensityCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridIon,
   basis2d = phaseBasisIon,
   basis1d = confBasis,
   moment = 0, --moment to compute
}

-- updater to compute phi from charge density
phiFromChargeDensityCalc = Updater.FemPoisson1D {
   onGrid = confGrid,
   basis = confBasis,
   sourceNodesShared = false, -- charge density is discontinous
   solutionNodesShared = false, -- solution is  discontinous
   -- left boundary is wall, so fix potential to ground
   bcLeft = { T = "D", V = 0.0 },
   -- right boundary is wall, so fix potential to ground
   bcRight = { T = "D", V = 0.0 },
}

-- Updater to compute electron momentum
momentumCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridElc,
   basis2d = phaseBasisElc,
   basis1d = confBasis,
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridIon,
   basis2d = phaseBasisIon,
   basis1d = confBasis,
   moment = 1,
}

-- Updater to compute electron momentum
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridElc,
   basis2d = phaseBasisElc,
   basis1d = confBasis,
   moment = 2,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGridIon,
   basis2d = phaseBasisIon,
   basis1d = confBasis,
   moment = 2,
}

-- updater to copy potential (1D field) to Hamiltonian (2D) field
-- (used in constructing full Hamiltonian, which also includes the KE
-- part)
copyTo2DElc = Updater.CopyNodalFields1D_2D {
   onGrid = phaseGridElc,
   sourceBasis = confBasis,
   targetBasis = phaseBasisElc
}
copyTo2DIon = Updater.CopyNodalFields1D_2D {
   onGrid = phaseGridIon,
   sourceBasis = confBasis,
   targetBasis = phaseBasisIon
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fldElc})
   end
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fldIon})
   end

   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyCopyBc(0, "lower")
      fld:applyCopyBc(0, "upper")
   end
   -- sync the distribution function across processors
   fldElc:sync()
   fldIon:sync()  
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

totalPtclElc = DataStruct.DynVector { numComponents = 1, }
totalPtclIon = DataStruct.DynVector { numComponents = 1, }

-- updater compute total number of electrons in domain
totalPtclCalcElc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- updater compute total number of ions in domain
totalPtclCalcIon = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- for DG fields common nodes not shared
}

----------------------
-- SOLVER UTILITIES --
----------------------

-- generic function to run an updater
function runUpdater(updater, currTime, timeStep, inpFlds, outFlds)
   updater:setCurrTime(currTime)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(currTime+timeStep)
end

-- function to calculate number density
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end
-- function to calculate momentum density
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end
-- function to calculate energy density
function calcPtclEnergy(calculator, curr, dt, distfIn, energyOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {energyOut})
end

-- function to compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
   -- number density
   calcNumDensity(numDensityCalcElc, curr, dt, distfElcIn, numDensityElc)
   calcNumDensity(numDensityCalcIon, curr, dt, distfIonIn, numDensityIon)
   -- momentum
   calcMomentum(momentumCalcElc, curr, dt, distfElcIn, momentumElc)
   calcMomentum(momentumCalcIon, curr, dt, distfIonIn, momentumIon)
   -- energy
   calcPtclEnergy(ptclEnergyCalcElc, curr, dt, distfElcIn, ptclEnergyElc)
   calcPtclEnergy(ptclEnergyCalcIon, curr, dt, distfIonIn, ptclEnergyIon)
end

-- function to copy 1D field to 2D field
function copyPhi(copier, curr, dt, phi1, phi2)
   runUpdater(copier, curr, dt, {phi1}, {phi2})
   phi2:sync()
end

-- function to update Vlasov equation
function updateVlasovEqn(pbSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- function to compute phi from number density
function calcPhiFromChargeDensity(curr, dt, distElcIn, distIonIn, phiOut)
   calcMoments(curr, dt, distElcIn, distIonIn)
   -- charge density: -rhoc/epsilon0
   chargeDensity:combine(-ionCharge/epsilon0, numDensityIon, -elcCharge/epsilon0, numDensityElc)
   local myS, myDt = runUpdater(phiFromChargeDensityCalc, curr, dt, {chargeDensity}, {phiOut})
   -- add in drive field
   runUpdater(calcDrivePhi, curr, dt, {}, {phiDrive})
   --phiOut:accumulate(1.0, phiDrive)
   phiOut:copy(phiDrive)
   return myS, myDt
end

-- compute hamiltonian for electrons
function calcHamiltonianElc(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DElc, curr, dt, phiIn, hamilOut)
   hamilOut:scale(elcCharge/elcMass)
   hamilOut:accumulate(1.0, hamilKeElc)
end
-- compute hamiltonian for ions
function calcHamiltonianIon(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DIon, curr, dt, phiIn, hamilOut)
   hamilOut:scale(ionCharge/ionMass)
   hamilOut:accumulate(1.0, hamilKeIon)
end

-- function to compute diagnostics
function calcDiagnostics(curr, dt)
   runUpdater(totalPtclCalcElc, curr, dt, {numDensityElc}, {totalPtclElc})
   runUpdater(totalPtclCalcIon, curr, dt, {numDensityIon}, {totalPtclIon})
end

----------------------------
-- Time-stepping routines --
----------------------------

function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   -- RK stage 1
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distfElc, hamilElc, distf1Elc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distfIon, hamilIon, distf1Ion)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 2
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 3
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   calcPhiFromChargeDensity(tCurr, myDt, distfElc, distfIon, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   return true, math.min(dtSuggestedElc, dtSuggestedIon)
end

-- make a duplicate in case we need it
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()

-- function to advance solution from tStart to tEnd
function runSimulation(tStart, tEnd, nFrames, initDt)
   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- the grand loop 
   while true do
      distfDupElc:copy(distfElc)
      distfDupIon:copy(distfIon)
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
   -- distribution functions
   distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
   distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)
   -- potential
   phi1d:write(string.format("phi_%d.h5", frameNum), tCurr)   
   -- moments
   numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)
   momentumElc:write(string.format("momentumElc_%d.h5", frameNum), tCurr)
   momentumIon:write(string.format("momentumIon_%d.h5", frameNum), tCurr)
   -- diagnostics
   totalPtclElc:write(string.format("totalPtclElc_%d.h5", frameNum), tCurr)
   totalPtclIon:write(string.format("totalPtclIon_%d.h5", frameNum), tCurr)
   ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)
   ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr)
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- -- apply initial conditions for electrons and ion
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
-- initialize KE parts of Hamiltonian
runUpdater(initHamilKeElc, 0.0, 0.0, {}, {hamilKeElc})
runUpdater(initHamilKeIon, 0.0, 0.0, {}, {hamilKeIon})

-- -- compute initial moments
calcMoments(0.0, 0.0, distfElc, distfIon)

-- -- calculate initial potential and Hamiltonians
applyBc(0.0, 0.0, distfElc, distfIon)
calcPhiFromChargeDensity(0.0, 0.0, distfElc, distfIon, phi1d)
calcHamiltonianElc(0.0, 0.0, phi1d, hamilElc)
calcHamiltonianIon(0.0, 0.0, phi1d, hamilIon)

--compute initial diagnostics
calcDiagnostics(0.0, 0.0)

-- write out initial fields
writeFields(0, 0.0)

-- run the whole thing
initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)

-- print some timing information
Lucee.logInfo(string.format("Total time in poisson solver = %g", phiFromChargeDensityCalc:totalAdvanceTime()) )
Lucee.logInfo(string.format(
		 "Total time in poisson bracket = %g", vlasovSolverElc:totalAdvanceTime()+vlasovSolverIon:totalAdvanceTime()))
Lucee.logInfo(string.format(
		 "Total time in number density calc = %g", numDensityCalcElc:totalAdvanceTime()+numDensityCalcIon:totalAdvanceTime()))
     
