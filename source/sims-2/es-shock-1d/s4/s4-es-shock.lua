-- Electrostatic shock problem

log = Lucee.logInfo

-- problem parameters
Te_Ti = 1.0 -- ratio of electron to ion temperaute
machNum = 3.0 -- Mach number computed from ion thermal speed
n0 = 1.0 -- initial number density

-- physical parameters
elcTemp = 1.0
elcMass = 1.0
elcCharge = -1.0

ionTemp = elcTemp/Te_Ti
ionMass = 1836.2
ionCharge = 1.0

-- permittivity of free space
epsilon0 = 1.0

-- thermal speeds
cs = math.sqrt(elcTemp/ionMass)
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- electron and ion drift speeds
elcDrift = machNum*cs
ionDrift = elcDrift -- no net current

-- domain size and simulation time
LX = 200*lambdaD
tEnd = 500.0/wpe
nFrames = 100

-- Resolution, time-stepping etc.
NX = 200
NV = 64
polyOrder = 2

cfl = 0.5/(2*polyOrder+1)

-- compute max thermal speed to set velocity space extents
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc
VL_ION, VU_ION = -6.0*vtIon-ionDrift, 6.0*vtIon+ionDrift  -- -6.0*vtIon, 6.0*vtIon

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Sound speed=%g", cs))
log(string.format("Mach number=%g", machNum))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Plasma frequency=%g", wpe))
log(string.format("Debye length=%g", lambdaD))
log(string.format("Cell size=%g", LX/NX))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Electron/Ion drift speed=%g", elcDrift))
log(string.format("Electron domain extents = [%g,%g]", VL_ELC, VU_ELC))
log(string.format("Ion domain extents = [%g,%g]", VL_ION, VU_ION))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- phase space grid for electrons
gridElc = Grid.RectCart2D {
   lower = {0.0, VL_ELC},
   upper = {LX, VU_ELC},
   cells = {NX, NV},
}
-- phase space grid for ions
gridIon = Grid.RectCart2D {
   lower = {0.0, VL_ION},
   upper = {LX, VU_ION},
   cells = {NX, NV},
}

-- create FEM nodal basis
basisElc = NodalFiniteElement2D.Serendipity {
   onGrid = gridElc,
   polyOrder = polyOrder,
}
basisIon = NodalFiniteElement2D.Serendipity {
   onGrid = gridIon,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- A generic function to run an updater.
--
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

-- updater to initialize distribution function
initDistfElc = Updater.ProjectOnNodalBasis2D {
   onGrid = gridElc,
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
		 local sloc = 0.5*LX
		 local fxv = maxwellian(n0, elcDrift, vtElc, v)
		 if x>sloc then
		    fxv = maxwellian(n0, -elcDrift, vtElc, v)
		 end
		 return fxv
	      end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})

-- updater to initialize distribution function
initDistfIon = Updater.ProjectOnNodalBasis2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
		 local sloc = 0.5*LX
		 local fxv = maxwellian(n0, ionDrift, vtIon, v)
		 if x>sloc then
		    fxv = maxwellian(n0, -ionDrift, vtIon, v)
		 end
		 return fxv
	      end
}
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})

-- extra fields for performing RK update
distfNewElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create field to store kinetic energy term in Hamiltonian
hamilKeElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
hamilKeIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilKeElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeElc, 0.0, 0.0, {}, {hamilKeElc})

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilKeIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeIon, 0.0, 0.0, {}, {hamilKeIon})

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.PoissonBracket {
   onGrid = gridElc,
   basis = basisElc,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}
vlasovSolverIon = Updater.PoissonBracket {
   onGrid = gridIon,
   basis = basisIon,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}

-- spatial grid
gridCS = Grid.RectCart1D {
   lower = {0.0},
   upper = {LX},
   cells = {NX},
}

-- spatial FEM nodal basis
basisCS = NodalFiniteElement1D.Lobatto {
   onGrid = gridCS,
   polyOrder = polyOrder,
}

-- Electron number density
numDensityElc = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}
-- Ion number density
numDensityIon = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}
-- Electron momentum
momentumElc = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}

-- charge density
chargeDensity = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}

-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 0,
}
numDensityCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

-- dynvector for total particle count
totalPtclElc = DataStruct.DynVector { numComponents = 1, }
totalPtclIon = DataStruct.DynVector { numComponents = 1, }

-- to compute total number of particles in domain
totalPtclCalcElc = Updater.IntegrateNodalField1D {
   --grid for updater
   onGrid = gridCS,
   --basis functions to use
   basis = basisCS,
   --are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}

totalPtclCalcIon = Updater.IntegrateNodalField1D {
   --grid for updater
   onGrid = gridCS,
   --basis functions to use
   basis = basisCS,
   --are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}

--set input fields
totalPtclCalcElc:setIn( {numDensityElc} )
--set output dynvector
totalPtclCalcElc:setOut( {totalPtclElc} )

--set input fields
totalPtclCalcIon:setIn( {numDensityIon} )
--set output dynvector
totalPtclCalcIon:setOut( {totalPtclIon} )

-- Updater to compute electron momentum
momentumCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 1,
}

-- Updater to compute electron momentum
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 2,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basisCS,
   -- desired moment (0, 1 or 2)
   moment = 2,
}

-- calculate number density at current time
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end

-- calculate momentum at current time
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end

-- calculate energy
function calcPtclEnergy(calculator, curr, dt, distfIn, energyOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {energyOut})
end

-- compute moments from distribution function
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

-- compute initial moments
calcMoments(0.0, 0.0, distfElc, distfIon)

-- field to store continous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}
phi1dDg = DataStruct.Field1D {
   onGrid = gridCS,
   numComponents = basisCS:numNodes(),
   ghost = {1, 1},
}

-- updater to copy 1D field to 2D field
copyTo2DElc = Updater.Copy1DTo2DNodalField {
   onGrid = gridElc,
}
copyTo2DIon = Updater.Copy1DTo2DNodalField {
   onGrid = gridIon,
}

-- function to copy 1D field to 2D field
function copyPhi(copier, curr, dt, phi1, phi2)
   return runUpdater(copier, curr, dt, {phi1}, {phi2})
end

-- updater to compute phi from charge density
phiFromChargeDensityCalc = Updater.FemPoisson1D {
   onGrid = gridCS,
   basis = basisCS,
   -- flag to indicate if nodes in src field are shared
   sourceNodesShared = false, -- charge density is discontinous
   -- left boundary is wall, so fix potential to ground
   bcLeft = { T = "D", V = 0.0 },
   -- right boundary is wall, so fix potential to ground
   bcRight = { T = "D", V = 0.0 },
}

-- update Vlasov equation, given appropriate updater 
function updateVlasovEqn(pbSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- compute phi from number density
function calcPhiFromChargeDensity(curr, dt, distElcIn, distIonIn, phiOut)
   calcMoments(curr, dt, distElcIn, distIonIn)
   -- charge density: -rhoc/epsilon0
   chargeDensity:combine(-ionCharge/epsilon0, numDensityIon, -elcCharge/epsilon0, numDensityElc)
   return runUpdater(phiFromChargeDensityCalc, curr, dt, {chargeDensity}, {phiOut})
end

-- create updater to initialize chi
copyCToD = Updater.CopyContToDisCont1D {
   onGrid = gridCS,
   basis = basisCS,
}
function copyPotential(tCurr, dt, cgIn, dgOut)
   runUpdater(copyCToD, tCurr, dt, {cgIn}, {dgOut})
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

-- A HACK
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val}
   end
end
function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   end
end

bcConst = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}
bcFunctionElc = BoundaryCondition.Function {
   components = getCountTbl(polyOrder),
   bc = function(x,v,z,t)
	   local sloc = 0.5*LX
	   local elcThermal = math.sqrt(elcTemp/elcMass)
	   local fv = maxwellian(n0, elcDrift, vtElc, v)
	   if x>sloc then
	      fv = maxwellian(n0, -elcDrift, vtElc, v)
	   end
	   if polyOrder == 1 then
	      return fv,fv,fv,fv
	   else
	      return fv,fv,fv,fv,fv,fv,fv,fv
	   end
	end,
}
bcFunctionIon = BoundaryCondition.Function {
   components = getCountTbl(polyOrder),
   bc = function(x,v,z,t)
	   local sloc = 0.5*LX
	   local ionThermal = math.sqrt(ionTemp/ionMass)
	   local fv = maxwellian(n0, ionDrift, vtIon, v)
	   if x>sloc then
	      fv = maxwellian(n0, -ionDrift, vtIon, v)
	   end
	   if polyOrder == 1 then
	      return fv,fv,fv,fv
	   else
	      return fv,fv,fv,fv,fv,fv,fv,fv
	   end
	end,
}

-- function to make make BC updaters
function makeBcObjElc()
   local bcLower = Updater.Bc2D {
      onGrid = gridElc,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc2D {
      onGrid = gridElc,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

function makeBcObjIon()
   local bcLower = Updater.Bc2D {
      onGrid = gridIon,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc2D {
      onGrid = gridIon,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

-- make objects to apply BCs
bcLowerElc, bcUpperElc = makeBcObjElc()
bcLowerIon, bcUpperIon = makeBcObjIon()

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon)
   for i,bc in ipairs({bcLowerElc, bcUpperElc}) do
      --runUpdater(bc, curr, dt, {}, {fldElc})
   end
   for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
      --runUpdater(bc, curr, dt, {}, {fldIon})
   end

   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyCopyBc(0, "lower")
      fld:applyCopyBc(0, "upper")
   end

   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyCopyBc(1, "lower")
      fld:applyCopyBc(1, "upper")
   end
end

-- calculate initial potential and Hamiltonians
applyBc(0.0, 0.0, distfElc, distfIon)
calcPhiFromChargeDensity(0.0, 0.0, distfElc, distfIon, phi1d)
calcHamiltonianElc(0.0, 0.0, phi1d, hamilElc)
calcHamiltonianIon(0.0, 0.0, phi1d, hamilIon)

function calcDiagnostics(curr, dt)
   totalPtclCalcElc:setCurrTime(curr)
   totalPtclCalcElc:advance(curr+dt)

   totalPtclCalcIon:setCurrTime(curr)
   totalPtclCalcIon:advance(curr+dt)
end

--compute initial diagnostics
calcDiagnostics(0.0, 0.0)

-- function to take a time-step using SSP-RK3 time-stepping scheme
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
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDupElc:copy(distfElc)
      distfDupIon:copy(distfIon)

      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 distfElc:copy(distfDupElc)
	 distfIon:copy(distfDupIon)
	 myDt = dtSuggested
      else
	 calcDiagnostics(tCurr, myDt)

	 tCurr = tCurr + myDt
	 myDt = dtSuggested
	 step = step + 1
	 if (tCurr >= tEnd) then
	    break
	 end
      end

   end
   return dtSuggested
end

-- Write out data frame 'frameNum' with at specified time 'tCurr'
function writeFields(frameNum, tCurr)
   distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
   distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)

   numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)

   momentumElc:write(string.format("momentumElc_%d.h5", frameNum), tCurr)
   momentumIon:write(string.format("momentumIon_%d.h5", frameNum), tCurr)

   totalPtclElc:write(string.format("totalPtclElc_%d.h5", frameNum), tCurr )
   totalPtclIon:write(string.format("totalPtclIon_%d.h5", frameNum), tCurr )

   ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)
   ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr)

   copyPotential(0.0, 0.0, phi1d, phi1dDg)
   phi1dDg:write(string.format("phi_%d.h5", frameNum), tCurr)
end

-- write out initial fields
writeFields(0, 0.0)

-- parameters to control time-stepping
tStart = 0.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
tFrame = (tEnd-tStart)/nFrames
tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end


