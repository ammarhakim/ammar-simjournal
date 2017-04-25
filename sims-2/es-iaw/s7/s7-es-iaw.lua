-- Vlasov-Maxwell solver: Units are SI

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
Tratio = 2.0
knumber = 0.5

n0 = 1.0
ionCharge = 1.0
ionMass = 1.0
elcCharge = -1.0 -- electron charge
elcMass = 1.0/1836.2 -- electron mass
epsilon0 = 1.0
gasGamma = 3

-- thermal speeds
ionTemp = 1.0
elcTemp = ionTemp/Tratio
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)

-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))

-- perturbation
pert = 1e-4
-- domain size and simulation time
LX = 2*Lucee.Pi/knumber

-- Resolution, time-stepping etc.
NX = 32
NV = 32
tStart = 0.0 -- start time 
tEnd = 12/(vtIon*knumber)
nFrames = 20

cfl = 0.5/(2*polyOrder+1)

-- set velocity space extents
VL_ION, VU_ION = -6.0*vtIon, 6.0*vtIon

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Cell size=%g", LX/NX))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Plasma frequency=%g", wpe))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Ion velocity domain extents = [%g,%g]", VL_ION, VU_ION))
log(string.format("Ion config domain extents = [%g,%g]", -LX/2, LX/2))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc2D.CartProd { cuts = {1,1} }
confDecomp = DecompRegionCalc1D.SubCartProd2D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid for ions
phaseGridIon = Grid.RectCart2D {
   lower = {-LX/2, VL_ION},
   upper = {LX/2, VU_ION},
   cells = {NX, NV},
   periodicDirs = {0},   
   decomposition = phaseDecomp,   
}

-- configuration space grid (same for electrons and ions)
confGrid = Grid.RectCart1D {
   lower = {-LX/2},
   upper = {LX/2},
   cells = {NX},
   periodicDirs = {0},
   decomposition = confDecomp,
}

phaseBasisIon = NodalFiniteElement2D.SerendipityElement {
   onGrid = phaseGridIon,
   polyOrder = polyOrder,
}

-- configuration-space basis functions (shared by both species)
confBasis = NodalFiniteElement1D.SerendipityElement {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

confBasisEs = NodalFiniteElement1D.Lobatto {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

distfNewIon = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- Ion number density
numDensityIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- Electron number density (used for Poisson solve)
numDensityElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- RHS of Poission equations
poissonRHS = DataStruct.Field1D {
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

-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- electron fluid 5 moment field
elcFluid = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 5*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
elcFluid1 = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 5*confBasis:numNodes(),
   ghost = {1, 1},
}
elcFluidNew = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 5*confBasis:numNodes(),
   ghost = {1, 1},
}

-- electrostatic potential calculated from number density
potential = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- Electric field calculated from potential
Ex = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- EM field
em = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
em:clear(0.0)
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
-- placeholder fields for source solve
elcFluidSource = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 5*confBasis:numNodes(),
   ghost = {1, 1},
}
emSourceFluid = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for adding to EM fields (this perhaps is not the best way to do it)
emSourceCurrent = DataStruct.Field1D {
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
function maxwellian(n0, vdrift_x, vt, vx)
   local v2 = (vx-vdrift_x)^2
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-v2/(2*vt^2))
end

-- updater to initialize ion distribution function
initDistfIon = Updater.ProjectOnNodalBasis2D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,vx,t)
      local alpha = pert -- perturbation
      local k = knumber 
      local nIon = n0*(1+alpha*math.cos(k*x))
      return maxwellian(nIon, 0.0, vtIon, vx)
   end
}

-- updater to initialize electron fluid 5 moment fields
initElcFluid = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local alpha = pert -- perturbation
      local k = knumber
      local nElc = n0
      local elcEnergy = nElc*elcTemp/(gasGamma-1)
      return elcMass*nElc, 0.0, 0.0, 0.0, elcEnergy
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- Updater for ion Vlasov equation

vlasovSolverIon = Updater.EigenNodalVlasov1X1V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   cfl = cfl,
   charge = ionCharge,
   mass = ionMass,
   polyOrder = polyOrder,
}

-- electron fluid equation object
elcFluidEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}

-- updater to solve Maxwell equations
elcFluidSlvr = Updater.NodalDgHyper1D {
   onGrid = confGrid,
   basis = confBasis,
   equation = elcFluidEqn,
   cfl = cfl,
}

SrcSlvr = Updater.DGExplicitIncrementFiveMomentSrc1D {
   onGrid = confGrid,
   basis = confBasis,
   numFluids = 1,
   charge = {elcCharge},
   mass = {elcMass},
   epsilon0 = epsilon0,
}

-- updater to compute phi from number density
phiFromNumDensityCalc = Updater.FemPoisson1D {
   onGrid = confGrid,
   basis = confBasisEs,
   sourceNodesShared = false, -- charge density is discontinous
   solutionNodesShared = false, -- solution is  discontinous
   periodicDirs = {0},
}

ExFromPhiCalc = Updater.NodalGradient1D {
   onGrid = confGrid,
   basis = confBasis,
}

-- Updater to compute ion number density
numDensityCalcIon = Updater.DistFuncMomentCalc1X1V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   moment = 0,
}

-- Updater to compute ion momentum
momentumCalcIon = Updater.DistFuncMomentCalc1X1V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 1,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1X1V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 2,
   scalarPtclEnergy = true,
}

-- This strange looking updater copies the currents into the EM source
-- field. Perhaps this is not the best way to do things, and one can
-- imagine a source updater which adds current sources to the dE/dt
-- Maxwell equation
copyToEmSource = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {0},
   targetComponents = {0},
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields


-- apply boundary conditions to distribution functions
function applyDistFuncBc(curr, dt, fldIon)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fldIon})
   end
   
   -- sync the distribution function across processors
   fldIon:sync()
end

-- apply BCs to EM and electron fluid fields
function applyConfigBc(curr, dt, elcFluid, EM)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {elcFluid, EM})
   end
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {elcFluid, EM})
   end

   -- sync electron fluid fields across processors
   elcFluid:sync()

   -- sync EM fields across processors
   EM:sync()
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

emEnergy = DataStruct.DynVector { numComponents = 1 }
emEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*(ex^2+ey^2+ez^2)
   end,
}
emEnergyCalc:setIn( {em} )
emEnergyCalc:setOut( {emEnergy} )

-- dynvector for total ptcl energy for ions
totalPtclEnergyIon = DataStruct.DynVector { numComponents = 1, }
-- to compute total particle energy for ions
totalPtclEnergyIonCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (u)
      return u
   end,
}
-- set input field
totalPtclEnergyIonCalc:setIn( {ptclEnergyIon} )
-- set output dynvector
totalPtclEnergyIonCalc:setOut( {totalPtclEnergyIon} )

-- dynvector for total ptcl energy for ions
totalPtclEnergyElc = DataStruct.DynVector { numComponents = 1, }
-- to compute total particle energy for ions
totalPtclEnergyElcCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (rho, rhou, rhov, rhow, Er)
      return Er
   end,
}
-- set input field
totalPtclEnergyElcCalc:setIn( {elcFluid} )
-- set output dynvector
totalPtclEnergyElcCalc:setOut( {totalPtclEnergyElc} )

ionNumInCell = DataStruct.DynVector { numComponents = confBasis:numNodes() }
ionNumInCellCalc = Updater.RecordFieldInCell1D {
   onGrid = confGrid,
   cellIndex = {NX/3},
}
ionNumInCellCalc:setIn( {numDensityIon} )
ionNumInCellCalc:setOut( {ionNumInCell} )

phiInCell = DataStruct.DynVector { numComponents = confBasis:numNodes() }
phiInCellCalc = Updater.RecordFieldInCell1D {
   onGrid = confGrid,
   cellIndex = {NX/3},
}
phiInCellCalc:setIn( {potential} )
phiInCellCalc:setOut( {phiInCell} )

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
-- function to calculate momentum density
function calcPtclEnergy(calculator, curr, dt, distfIn, ptclEnergyOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {ptclEnergyOut})
end

-- function to compute moments from distribution function
function calcMoments(curr, dt, distfIonIn)
   calcNumDensity(numDensityCalcIon, curr, dt, distfIonIn, numDensityIon)
   calcMomentum(momentumCalcIon, curr, dt, distfIonIn, momentumIon)
   calcPtclEnergy(ptclEnergyCalcIon, curr, dt, distfIonIn, ptclEnergyIon)
end

-- function to update Vlasov equation
function updateVlasovEqn(vlasovSlvr, curr, dt, distfIn, emIn, distfOut)
   return runUpdater(vlasovSlvr, curr, dt, {distfIn, emIn}, {distfOut})
end

-- compute phi from number density
function calcPhiFromNumDensity(curr, dt, numDensIonIn, numDensElcIn, phiOut)
   phiOut:clear(0.0)
   -- RHS = -rho_c/epsilon0 (fluids evolve electron mass density)
   poissonRHS:combine(-ionCharge/epsilon0, numDensIonIn, -elcCharge/(elcMass*epsilon0), numDensElcIn)
   return runUpdater(phiFromNumDensityCalc, curr, dt, {poissonRHS}, {phiOut})
end

-- function to calculate electric field from potential
function calcExFromPhi(calculator, curr, dt, potentialIn, ExOut)
   potentialIn:scale(-1.0)
   return runUpdater(calculator, curr, dt, {potentialIn}, {ExOut})
end

-- solve electron fluid equation
function updateElcFluidEqn(curr, dt, elcFluidIn, elcFluidOut)
   return runUpdater(elcFluidSlvr, curr, dt, {elcFluidIn}, {elcFluidOut})
end

-- solve electron fluid sources
function updateSource(curr, dt, elcFluidIn, emIn, elcFluidOut, emOut)
   return runUpdater(SrcSlvr, curr, dt, {elcFluidIn, emIn}, {elcFluidOut, emOut})
end

-- function to compute diagnostics
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({emEnergyCalc, totalPtclEnergyIonCalc, totalPtclEnergyElcCalc, phiInCellCalc, ionNumInCellCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

----------------------------
-- Time-stepping routines --
----------------------------

-- take single RK step
function rkStage(tCurr, dt, ionIn, elcFluidIn, emIn, ionOut, elcFluidOut, emOut)
   -- update distribution functions and homogenous Maxwell equations
   local stElc, dtElc = updateElcFluidEqn(tCurr, dt, elcFluidIn, elcFluidOut)
   local stIon, dtIon = updateVlasovEqn(vlasovSolverIon, tCurr, dt, ionIn, emIn, ionOut)

   -- calculate source contribution
   updateSource(tCurr, dt, elcFluidIn, emIn, elcFluidSource, emSourceFluid)

   -- accumulate source contribution, but only for the electron fluid
   -- this is because we're using Poisson's equation for the electric field
   elcFluidOut:accumulate(dt, elcFluidSource)

   -- copying electron mass density so can be used for Poisson solve
   runUpdater(copyToEmSource, tCurr, dt, {elcFluidOut}, {numDensityElc})

   -- calculating updated ion number density
   calcNumDensity(numDensityCalcIon, tCurr, dt, ionOut, numDensityIon)

   -- calculating phi from Poisson's equation 
   -- note that the mass factor in the electron number density is scaled out within the function
   calcPhiFromNumDensity(tCurr, dt, numDensityIon, numDensityElc, potential)

   -- calculating -grad(phi) to obtain the electric field
   calcExFromPhi(ExFromPhiCalc, tCurr, dt, potential, Ex)

   -- copying the updated electric field into the appropriate data structure
   runUpdater(copyToEmSource, tCurr, dt, {Ex}, {emOut})
    if (stElc == false) or (stIon == false) then
      return false, math.min(dtElc, dtIon)
   end
   return true, math.min(dtElc, dtIon)
end

function rk3(tCurr, myDt)
   local myStatus, myDtSuggested

   -- RK stage 1
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distfIon, elcFluid, em, distf1Ion, elcFluid1, em1)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Ion)
   applyConfigBc(tCurr, dt, elcFluid1, em1)

   -- RK stage 2
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Ion, elcFluid1, em1, distfNewIon, elcFluidNew, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   elcFluid1:combine(3.0/4.0, elcFluid, 1.0/4.0, elcFluidNew)
   em1:combine(3.0/4.0, em, 1.0/4.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Ion)
   applyConfigBc(tCurr, dt, elcFluid1, em1)   

   -- RK stage 3
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Ion, elcFluid1, em1, distfNewIon, elcFluidNew, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   elcFluid1:combine(1.0/3.0, elcFluid, 2.0/3.0, elcFluidNew)
   em1:combine(1.0/3.0, em, 2.0/3.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Ion)
   applyConfigBc(tCurr, dt, elcFluid1, em1)   

   distfIon:copy(distf1Ion)
   elcFluid:copy(elcFluid1)
   em:copy(em1)

   return true, myDtSuggested
end

-- make duplicates in case we need them
distfDupIon = distfIon:duplicate()
elcFluidDup = elcFluid:duplicate()
emDup = em:duplicate()

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
      distfDupIon:copy(distfIon)
      elcFluidDup = elcFluid:duplicate()
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

	 distfIon:copy(distfDupIon)
         elcFluidDup = elcFluid:duplicate()
	 em:copy(emDup)
      else
	 -- compute diagnostics
         calcPtclEnergy(ptclEnergyCalcIon, tCurr, myDt, distfIon, ptclEnergyIon)
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
   distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)

   -- EM field
   em:write(string.format("em_%d.h5", frameNum), tCurr)   
   elcFluid:write(string.format("elcFluid_%d.h5", frameNum), tCurr)

   -- compute moments and write them out
   calcMoments(tCurr, 0.0, distfIon)
   numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)
   momentumIon:write(string.format("momentumIon_%d.h5", frameNum), tCurr)
   ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr)
   -- diagnostics
   emEnergy:write( string.format("emEnergy_%d.h5", frameNum) )
   totalPtclEnergyIon:write( string.format("totalPtclEnergyIon_%d.h5", frameNum) )
   totalPtclEnergyElc:write( string.format("totalPtclEnergyElc_%d.h5", frameNum) )
   ionNumInCell:write( string.format("ionNumInCell_%d.h5", frameNum) )
   phiInCell:write( string.format("phiInCell_%d.h5", frameNum) )
end


----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- apply initial conditions for electrons and ion
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
runUpdater(initElcFluid, 0.0, 0.0, {}, {elcFluid})

applyConfigBc(0.0, 0.0, elcFluid, em)
applyDistFuncBc(0.0, 0.0, distfIon)

runUpdater(copyToEmSource, 0.0, 0.0, {elcFluid}, {numDensityElc})
calcNumDensity(numDensityCalcIon, 0.0, 0.0, distfIon, numDensityIon)

calcPhiFromNumDensity(0.0, 0.0, numDensityIon, numDensityElc, potential)
calcExFromPhi(ExFromPhiCalc, 0.0, 0.0, potential, Ex)
runUpdater(copyToEmSource, 0.0, 0.0, {Ex}, {em})

applyConfigBc(0.0, 0.0, elcFluid, em)

-- apply BCs
applyDistFuncBc(0.0, 0.0, distfIon)
applyConfigBc(0.0, 0.0, elcFluid, em)
-- compute initial diagnostics
calcMoments(0.0, 0.0, distfIon)
calcDiagnostics(0.0, 0.0)
-- write out initial fields
writeFields(0, 0.0)

-- run the whole thing
initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)

