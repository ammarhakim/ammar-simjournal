-- Neutral Boltzmann-BGK mode

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order

-- left/right state for shock
nl, ul, pl = 1.0, 0.75, 1.0
nr, ur, pr = 1.0, -0.75, 0.75

Lx = 1.0 -- domain size
mfp = 0.01*Lx -- mean-free path

-- thermal velocity to give same energy as in fluid internal energy
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- use left state as reference
nu = vThermal/mfp -- collision frequency

tStart = 0.0
tEnd = 0.15
nFrames = 5

VL, VU = -6.0*vThermal, 6.0*vThermal -- velocity space extents

-- Resolution, time-stepping etc.
NX = 128
NV = 32

cfl = 0.5/(2*polyOrder+1)

-- print some diagnostics
log(string.format("tEnd = %g,  nFrames=%d", tEnd, nFrames))
log(string.format("Mean-free path = %g", mfp))
log(string.format("Collision frequency = %g", nu))
log(string.format("Time between collisions = %g", 1/nu))
log(string.format("vThermal_l = %g, vThermal_r = %g", vThermal_l, vThermal_r))
log(string.format("\n"))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc2D.CartProd { cuts = {1,1} }
confDecomp = DecompRegionCalc1D.SubCartProd2D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid
phaseGrid = Grid.RectCart2D {
   lower = {0.0, VL},
   upper = {Lx, VU},
   cells = {NX, NV},
   decomposition = phaseDecomp,
}

-- configuration space grid
confGrid = Grid.RectCart1D {
   lower = {0.0},
   upper = {Lx},
   cells = {NX},
   decomposition = confDecomp,
}

-- phase-space basis functions
phaseBasis = NodalFiniteElement2D.SerendipityElement {
   onGrid = phaseGrid,
   polyOrder = polyOrder,
}
-- configuration-space basis functions
confBasis = NodalFiniteElement1D.LagrangeTensor {
   onGrid = confGrid,
   polyOrder = polyOrder,
   nodeLocation = "uniform",
}

-- distribution function for 
distf = DataStruct.Field2D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
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

-- Maxwellian distribution corresponding
distfMaxwell = DataStruct.Field2D {
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

-- number density
numDensity = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- momentum
momentum = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergy = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
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

-- updater to initialize electron distribution function
initDistf = Updater.ProjectOnNodalBasis2D {
   onGrid = phaseGrid,
   basis = phaseBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,v,z,t)
      local n, u, vt = nl, ul, vThermal_l
      if (x>0.5) then
	 n, u, vt = nr, ur, vThermal_r
      end
      return maxwellian(n, u, vt, v)
   end
}

-- updater to initialize hamiltonian
initHamilKE = Updater.EvalOnNodes2D {
   onGrid = phaseGrid,
   basis = phaseBasis,
   shareCommonNodes = false,
   evaluate = function (x,v,z,t)
      return v^2/2
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- Updater for electron Vlasov equation
vlasovSolver = Updater.PoissonBracket {
   onGrid = phaseGrid,
   basis = phaseBasis,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   hamilNodesShared = false, -- Hamiltonian is not continuous
   zeroFluxDirections = {1},
}

-- updater to initialize Maxwellian
calcMaxwellian = Updater.MaxwellDistInit1X1V {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis
}

-- Updater to compute number density
numDensityCalc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGrid,
   basis2d = phaseBasis,
   basis1d = confBasis,
   moment = 0,
}
-- Updater to compute momentum
momentumCalc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGrid,
   basis2d = phaseBasis,
   basis1d = confBasis,
   moment = 1,
}
-- Updater to compute particle energy
ptclEnergyCalc = Updater.DistFuncMomentCalc1D {
   onGrid = phaseGrid,
   basis2d = phaseBasis,
   basis1d = confBasis,
   moment = 2,
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields

-- apply boundary conditions to distribution functions
function applyDistFuncBc(curr, dt, fld)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fld})
   end
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   fld:sync()
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

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

-- calculate number density
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end
-- calculate momentum density
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end
-- calculate ptcl energy
function calcPtclEnergy(calculator, curr, dt, distfIn, energyOut)
   runUpdater(calculator, curr, dt, {distfIn}, {energyOut})
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   calcNumDensity(numDensityCalc, curr, dt, distfIn, numDensity)
   calcMomentum(momentumCalc, curr, dt, distfIn, momentum)
   calcPtclEnergy(ptclEnergyCalc, curr, dt, distfIn, ptclEnergy)
end

-- update Vlasov equation
function updateVlasovEqn(vlasovSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(vlasovSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- compute diagnostics
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

----------------------------
-- Time-stepping routines --
----------------------------

-- take single RK step
function rkStage(tCurr, dt, distfIn, hamilIn, distfOut)
   local status, suggestedDt = updateVlasovEqn(vlasovSolver, tCurr, dt, distfIn, hamilIn, distfOut)
   -- add contribution from BGK collision operator
   calcMoments(tCurr, dt, distfIn)
   runUpdater(calcMaxwellian, tCurr, dt, {numDensity, momentum, ptclEnergy}, {distfMaxwell})
   distfOut:accumulate(nu*dt, distfMaxwell, -nu*dt, distfIn)
   
   return status, suggestedDt
end

function rk3(tCurr, myDt)
   local myStatus, myDtSuggested
   -- RK stage 1
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf, hamilKE, distf1)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   applyDistFuncBc(tCurr, dt, distf1)

   -- RK stage 2
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1, hamilKE, distfNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyDistFuncBc(tCurr, dt, distf1)

   -- RK stage 3
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1, hamilKE, distfNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyDistFuncBc(tCurr, dt, distf1)

   distf:copy(distf1)

   return true, myDtSuggested
end

-- make duplicates in case we need them
distfDup = distf:duplicate()

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

-- Write out data frame 'frameNum' with at specified time 'tCurr'
function writeFields(frameNum, tCurr)
   distf:write(string.format("distf_%d.h5", frameNum), tCurr)
   numDensity:write(string.format("numDensity_%d.h5", frameNum), tCurr)
   momentum:write(string.format("momentum_%d.h5", frameNum), tCurr)
   ptclEnergy:write(string.format("ptclEnergy_%d.h5", frameNum), tCurr)
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------

-- setup simulation
runUpdater(initDistf, 0.0, 0.0, {}, {distf})
runUpdater(initHamilKE, 0.0, 0.0, {}, {hamilKE})
applyDistFuncBc(0.0, 0.0, distf)
calcMoments(0.0, 0.0, distf)
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

-- run the whole thing
initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)

-- print some timing information
log(string.format("Total time in vlasov solver for electrons = %g", vlasovSolver:totalAdvanceTime()))
log(string.format("Total time moment computations = %g",
		  numDensityCalc:totalAdvanceTime()+momentumCalc:totalAdvanceTime()+ptclEnergyCalc:totalAdvanceTime()))
