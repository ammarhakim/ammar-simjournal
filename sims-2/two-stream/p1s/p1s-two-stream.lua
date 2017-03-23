-- Vlasov-Poisson solver in Poisson-Bracket formulation

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
knumber = 0.5 -- wave-number
elcCharge = -1.0 -- signed electron charge
elcMass = 1.0 -- electron mass
elVTerm = 0.1 -- electron thermal velocity
epsilon0 = 1.0 -- permittivity of free space
vDrift = 1.0 -- dfirft velocity
perturbation = 1.0e-6 -- distribution function perturbation

-- resolution and time-stepping
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber -- configuration space extents
VL, VU = -6.0, 6.0 -- velocity space extents (this is in units of vthermal for electrons)
NX, NV = 64, 16 -- mesh size

cfl = 0.2 -- CFL number
tStart = 0.0 -- start time 
tEnd = 50 -- end time
nFrames = 10 -- number of output frames to write

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition objects for phase-space and configuration space
phaseDecomp = DecompRegionCalc2D.CartProd { cuts = {2,2} }
confDecomp = DecompRegionCalc1D.SubCartProd2D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase-space domain
phaseGrid = Grid.RectCart2D {
   lower = {XL, VL},
   upper = {XU, VU},
   cells = {NX, NV},
   decomposition = phaseDecomp,
   periodicDirs = {0},
}
-- configuration space grid
confGrid = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
   decomposition = confDecomp,
   periodicDirs = {0},
}

-- phase-space basis functions
phaseBasis = NodalFiniteElement2D.SerendipityElement {
   onGrid = phaseGrid,
   polyOrder = polyOrder,
}
-- configuration space basis functions
confBasis = NodalFiniteElement1D.Lobatto {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = 1*phaseBasis:numNodes(),
   ghost = {1, 1},
}
-- distribution function for electrons
distf = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = 1*phaseBasis:numNodes(),
   ghost = {1, 1},
}
-- backup copy of distribution function in case we need to repeat time-step
distfDup = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = 1*phaseBasis:numNodes(),
   ghost = {1, 1},
}

-- extra fields for performing RK update
distfNew = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
   ghost = {1, 1},
}
distf1 = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
   ghost = {1, 1},
}

-- number density of ions and electrons
numDensityIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
numDensity = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- ptcl energy
ptclEnergy = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 1*confBasis:numNodes(),
   ghost = {1, 1},
}

-- field to store electrostatic potential
phi1d = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- Hamiltonian
hamil = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
   ghost = {1, 1},
}
-- Kinetic-energy term in Hamiltonian
hamilKE = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
   ghost = {1, 1},
}

--------------------------------
-- INITIAL CONDITION UPDATERS --
--------------------------------

-- Maxwellian with specified thermal and drift speeds
function maxwellianDistf(vt, vDrift, x, v)
   return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vDrift )^2/(2*vt^2))
end
-- Two-stream Maxwellian with specified thermal and drift speeds
function twoStreamDistf(vt, vDrift, x, v)
   return 1/math.sqrt(8*Lucee.Pi*vt^2)*(math.exp(-(v-vDrift )^2/(2*vt^2))+math.exp(-(v+vDrift)^2/(2*vt^2)))
end

-- updater to initialize ion distribution function
initDistfIon = Updater.EvalOnNodes2D {
   onGrid = phaseGrid,
   -- basis functions to use
   basis = phaseBasis,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization   
   evaluate = function(x,y,z,t)
      return maxwellianDistf(elVTerm, 0.0, x, y)
   end
}
-- updater to initialize electron distribution function
initDistf = Updater.EvalOnNodes2D {
   onGrid = phaseGrid,
   -- basis functions to use
   basis = phaseBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
      local alpha = perturbation
      local k = knumber
      return (1+alpha*math.cos(k*x))*twoStreamDistf(elVTerm, vDrift, x, y)
   end
}

-- updater to initialize hamiltonian
initHamilKE = Updater.EvalOnNodes2D {
   onGrid = phaseGrid,
   -- basis functions to use
   basis = phaseBasis,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- updater for Vlasov equation for electrons
pbSlvr = Updater.PoissonBracket {
   onGrid = phaseGrid,
   basis = phaseBasis,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   hamilNodesShared = false, -- Hamiltonian is not continuous
   zeroFluxDirections = {1},
}

-- updater to compute phi from number density
phiFromNumDensityCalc = Updater.FemPoisson1D {
   onGrid = confGrid,
   basis = confBasis,
   sourceNodesShared = false, -- source is DG field
   solutionNodesShared = false, -- solution is DG field
   -- periodic directions
   periodicDirs = {0},
}

-- updater to compute number density
numDensityCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid
   onGrid = phaseGrid,
   -- 2D phase-space basis functions
   basis2d = phaseBasis,
   -- 1D spatial basis functions
   basis1d = confBasis,
   -- desired moment (0, 1 or 2)
   moment = 0,
}
-- updater to compute ptcl energy
ptclEnergyCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = phaseGrid,
   -- 2D phase-space basis functions
   basis2d = phaseBasis,
   -- 1D spatial basis functions
   basis1d = confBasis,
   -- desired moment (0, 1 or 2)
   moment = 2,
}

-- updater to copy potential (1D field) to Hamiltonian (2D) field
-- (used in constructing full Hamiltonian, which also includes the KE
-- part)
copyTo2D = Updater.CopyNodalFields1D_2D {
   onGrid = phaseGrid,
   sourceBasis = confBasis,
   targetBasis = phaseBasis
}

-- scaling positivity limiter
scalingPositivityLimiter = Updater.ScalingPositivityLimiter1X1V {
   onGrid = phaseGrid,
   basis  = phaseBasis,
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields

-- function to apply boundary conditions to specified field
function applyPhaseBc(fld, tCurr, myDt)
   local bcList = {}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- sync ghost cells
   fld:sync()
end
function applyConfBc(fld, tCurr, myDt)
   local bcList = {}
   for i,bc in ipairs(bcList) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- sync ghost cells
   fld:sync()
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

-- dynvectors for various diagnostics
totalNumDensity = DataStruct.DynVector { numComponents = 1, }
totalPtclEnergy = DataStruct.DynVector { numComponents = 1, }
fieldEnergy = DataStruct.DynVector { numComponents = 1, }

-- updater to compute integral of various fields
totalFieldCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- for DG fields common nodes not shared
   integrand = function (er)
      return er
   end,
}
-- updater to compute field energy
fieldEnergyCalc = Updater.NormGrad1D {
   onGrid = confGrid,
   basis = confBasis,
}

----------------------
-- SOLVER UTILITIES --
----------------------
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

-- update Poisson bracket operator
function poissonBracket(curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- compute phi from number density
function calcPhiFromNumDensity(curr, dt, numDensIn, phiOut)
   phiOut:clear(0.0)
   -- compute total charge (ion charge is assumed to be 1.0)
   numDensIn:accumulate(-1.0, numDensityIon)
   numDensIn:scale(1/epsilon0)
   applyConfBc(numDensIn)
   
   return runUpdater(phiFromNumDensityCalc, curr, dt, {numDensIn}, {phiOut})
end

-- function to copy 1D field to 2D field
function copyPhi(curr, dt, phi1, phi2)
   runUpdater(copyTo2D, curr, dt, {phi1}, {phi2})
   phi2:sync()
end

-- calculate number density
function calcNumDensity(curr, dt, distfIn, numDensOut)
   runUpdater(numDensityCalc, curr, dt, {distfIn}, {numDensOut})
   applyConfBc(numDensOut)
end

-- calculate ptcl energy
function calcPtclEnergy(curr, dt, distfIn, energyOut)
   runUpdater(ptclEnergyCalc, curr, dt, {distfIn}, {energyOut})
   applyConfBc(energyOut)
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   calcNumDensity(curr, dt, distfIn, numDensity)
   calcPtclEnergy(curr, dt, distfIn, ptclEnergy)
end

-- compute hamiltonian
function calcHamiltonian(curr, dt, distIn, hamilOut)
   calcMoments(curr, dt, distIn)
   calcPhiFromNumDensity(curr, dt, numDensity, phi1d)
   hamilOut:clear(0.0)
   copyPhi(curr, dt, phi1d, hamilOut) -- potential energy contribution
   hamilOut:scale(elcCharge/elcMass)
   hamilOut:accumulate(1.0, hamilKE)
   hamilOut:sync()
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   runUpdater(totalFieldCalc, curr, dt, {numDensity}, {totalNumDensity})
   runUpdater(totalFieldCalc, curr, dt, {ptclEnergy}, {totalPtclEnergy})
   runUpdater(fieldEnergyCalc,  curr, dt, {phi1d}, {fieldEnergy})
end

-- write data to H5 files
function writeFields(frame, t)
   distf:write( string.format("distf_%d.h5", frame), t )
   numDensity:write( string.format("numDensity_%d.h5", frame), t )
   phi1d:write( string.format("phi1d_%d.h5", frame), t )
   hamil:write( string.format("hamil_%d.h5", frame), t )
   totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame), t )
   totalNumDensity:write( string.format("totalNumDensity_%d.h5", frame), t ) -- actually rho_c/epsilon_0
   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame), t )
end

----------------------------
-- Time-stepping routines --
----------------------------

function rk3(tCurr, myDt)
   local status, dtSuggested
   -- RK stage 1
   status, dtSuggested = poissonBracket(tCurr, myDt, distf, hamil, distf1)
   if (status == false) then
      return status, dtSuggested
   end
   applyPhaseBc(distf1)
   runUpdater(scalingPositivityLimiter, tCurr, myDt, {}, {distf1})
   calcHamiltonian(tCurr, myDt, distf1, hamil)
   
   -- RK stage 2
   status, dtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   if (status == false) then
      return status, dtSuggested
   end
   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyPhaseBc(distf1)
   runUpdater(scalingPositivityLimiter, tCurr, myDt, {}, {distf1})
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 3
   status, dtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   if (status == false) then
      return status, dtSuggested
   end
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyPhaseBc(distf1)
   runUpdater(scalingPositivityLimiter, tCurr, myDt, {}, {distf1})
   distf:copy(distf1)
   calcHamiltonian(tCurr, myDt, distf, hamil)

   return status, dtSuggested
end

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
      distfDup:copy(distf)
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
	 distf:copy(distfDup)
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

----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- apply initial conditions for electrons and ion
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
runUpdater(initDistf, 0.0, 0.0, {}, {distf})
applyPhaseBc(distfIon)
applyPhaseBc(distf)

-- initialize KE part of Hamiltonian
runUpdater(initHamilKE, 0.0, 0.0, {}, {hamilKE})
applyPhaseBc(hamilKE)

-- compute number density of ions
calcNumDensity(0.0, 0.0, distfIon, numDensityIon)

-- calculate initial Hamiltonian
calcHamiltonian(0.0, 0.0, distf, hamil)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)

-- write initial conditions
writeFields(0, 0.0)
distfIon:write("distfIon.h5") -- ions don't evolve, so write only once
hamilKE:write("hamilKE.h5") -- KE part of Hamiltonian doesn't evolve, so write it out once

initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)
