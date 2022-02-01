-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.2/(2*polyOrder+1)

-- wave-number
knumber = 0.5
-- perpendicular wave number
kperp = math.sqrt(0.01)
-- normalized beta (electron plasma is 2*beta*me/mi)
beta = 0.1

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 16, 32

-- run specified updater
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

-- phase space grid 
grid = Grid.RectCart2D {
   lower = {XL, VL},
   upper = {XU, VU},
   cells = {NX, NV},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   onGrid = grid,
   polyOrder = polyOrder,
}

-- equilibrium distribution function
distfEquil = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

function maxwellian(x,v)
   return 1/math.sqrt(2*Lucee.Pi)*math.exp(-v^2/2)
end

-- updater to initialize distribution function
initDistfEquil = Updater.EvalOnNodes2D {
   onGrid = grid,
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 local v = y
		 return maxwellian(x,y)
	      end
}
initDistfEquil:setOut( {distfEquil} )
initDistfEquil:advance(0.0) -- time is irrelevant

-- write this out
distfEquil:write("distfEquil.h5")

-- perturbed distribution function
distf = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- extra fields for performing RK update
distfNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- to store increments from Poisson Bracket updaters
incrDistf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
incrDistf2 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- Hamiltonian
hamil = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- Kinetic energy term in Hamiltonian
hamilKE = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- updater to initialize hamiltonian
initHamilKE = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
initHamilKE:setOut( {hamilKE} )
-- initialize potential
initHamilKE:advance(0.0) -- time is irrelevant
hamilKE:applyPeriodicBc(0)

-- write this out a diagnostic
hamilKE:write("hamilKE.h5")

-- updater to initialize distribution function
initDistf = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 return 1e-6*math.cos(knumber*x)*maxwellian(x,y)
	      end
}
initDistf:setOut( {distf} )
initDistf:advance(0.0) -- time is irrelevant

-- updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   -- compute only increments
   onlyIncrement = true,
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.Lobatto {
   onGrid = grid_1d,
   polyOrder = polyOrder,
}

-- moments of perturbed distribution function
numDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
momentum = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}
momentumCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 1,
}

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   -- compute number density and momentum
   runUpdater(numDensityCalc, curr, dt, {distfIn}, {numDensity})
   runUpdater(momentumCalc, curr, dt, {distfIn}, {momentum})
   -- apply BCs to these fields
   numDensity:applyPeriodicBc(0)
   momentum:applyPeriodicBc(0)
end

-- field to store continous potentials in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}
Apar1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}

-- field to store continous potentials in 2D
phi2d = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}
Apar2d = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}

-- updater to copy 1D field to 2D field
copyTo2D = Updater.Copy1DTo2DNodalField {
   onGrid = grid,
}

-- function to copy 1D fields to 2D fields
function copyFieldTo2D(curr, dt, f1d, f2d)
   runUpdater(copyTo2D, curr, dt, {f1d}, {f2d})
end

-- updater to compute continous field from discontinious fields
contFromDisContCalc = Updater.ContFromDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}

-- function to calculate EM fields from moments
function calcEMField(curr, dt, numIn, momIn)
   -- calculate electrostatic potential
   runUpdater(contFromDisContCalc, curr, dt, {numIn}, {phi1d})
   phi1d:scale(-1/kperp^2)
   -- calculate parallel magnetic potential
   runUpdater(contFromDisContCalc, curr, dt, {momIn}, {Apar1d})
   Apar1d:scale(-beta/(kperp^2+beta))
end

-- updater to compute perturbed Hamiltonian
pertHamilCalc = Updater.LinEmGke1DPertHamil {
   onGrid = grid,
   basis = basis,
   charge = -1.0, -- species charge
   mass = 1.0, -- species mass
}

-- function to compute total linearized Hamiltonian
function calcHamiltonian(curr, dt, phi1dIn, Apar1dIn, hamilOut)
   -- calculate 2D fields from 1D fields
   copyFieldTo2D(curr, dt, phi1dIn, phi2d)
   copyFieldTo2D(curr, dt, Apar1dIn, Apar2d)

   hamilOut:clear(0.0)
   -- compute perturbed Hamiltonian
   runUpdater(pertHamilCalc, curr, dt, {phi2d, Apar2d}, {hamilOut})
   -- accumulate free-streaming contribution
   hamilOut:accumulate(1.0, hamilKE)
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }

-- to compute field energy
fieldEnergyCalc = Updater.NormGrad1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
}

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   runUpdater(fieldEnergyCalc, curr, dt, {phi1d}, {fieldEnergy})
end

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

function updateVlasovEqn(curr, dt, H0, Ht, distIn, distOut)
   -- compute increment from {f1, H0} term
   local s1, dt1 = runUpdater(pbSlvr, curr, dt, {distIn, H0}, {incrDistf1})
   -- compute increment from {f0, H0+H1} term
   local s2, dt2 = runUpdater(pbSlvr, curr, dt, {distfEquil, Ht}, {incrDistf2})

   local dtMin = math.min(dt1, dt2)
   -- check if either of the steps failed
   if ((s1 == false) or (s2 == false)) then
      return false, dtMin
   end

   -- combine terms
   distOut:combine(1.0, distIn, dt, incrDistf1, dt, incrDistf2)
   return true, dtMin
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1
   status, dtSuggested = updateVlasovEqn(tCurr, myDt, hamilKE, hamil, distf, distf1)
   if (status == false) then
      return status, dtSuggested
   end
   applyBc(distf1)
   -- compute Hamiltonian for use in next stage
   calcMoments(tCurr, myDt, distf1)
   calcEMField(tCurr, myDt, numDensity, momentum)
   calcHamiltonian(tCurr, myDt, phi1d, Apar1d, hamil)

   -- RK stage 2
   status, dtSuggested = updateVlasovEqn(tCurr, myDt, hamilKE, hamil, distf1, distfNew)
   if (status == false) then
      return status, dtSuggested
   end
   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(distf1)
   -- compute Hamiltonian for use in next stage
   calcMoments(tCurr, myDt, distf1)
   calcEMField(tCurr, myDt, numDensity, momentum)
   calcHamiltonian(tCurr, myDt, phi1d, Apar1d, hamil)

   -- RK stage 3
   status, dtSuggested = updateVlasovEqn(tCurr, myDt, hamilKE, hamil, distf1, distfNew)
   if (status == false) then
      return status, dtSuggested
   end
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyBc(distf1)
   distf:copy(distf1)
   -- compute Hamiltonian for use in next time-step
   calcMoments(tCurr, myDt, distf)
   calcEMField(tCurr, myDt, numDensity, momentum)
   calcHamiltonian(tCurr, myDt, phi1d, Apar1d, hamil)

   return status, dtSuggested
end

-- make a duplicate in case we need it
distfDup = distf:duplicate()

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDup:copy(distf)

      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 distf:copy(distfDup)
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

-- these fields are only for I/O and are not used anywhere else
phiDG = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
AparDG = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- create updater to initialize chi
copyCToD = Updater.CopyContToDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}

-- write data to H5 files
function writeFields(frame, tm)
   distf:write( string.format("distf_%d.h5", frame), tm )
   numDensity:write( string.format("numDensity_%d.h5", frame), tm )
   momentum:write( string.format("momentum_%d.h5", frame), tm )

   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame), tm )

   -- copy CG to DG fields before writing them out (makes it easier for plotting)
   runUpdater(copyCToD, 0.0, 0.0, {phi1d}, {phiDG})
   runUpdater(copyCToD, 0.0, 0.0, {Apar1d}, {AparDG})

   phiDG:write( string.format("phi_%d.h5", frame), tm )
   AparDG:write( string.format("Apar_%d.h5", frame), tm )
end

-- Compute initial set of fields
calcMoments(0.0, 0.0, distf)
calcEMField(0.0, 0.0, numDensity, momentum)
calcHamiltonian(0.0, 0.0, phi1d, Apar1d, hamil)
calcDiagnostics(0.0, 0.0)

-- write initial conditions
writeFields(0, 0.0)

-- parameters to control time-stepping
tStart = 0.0
tEnd = 10.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 4
tFrame = (tEnd-tStart)/nFrames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end