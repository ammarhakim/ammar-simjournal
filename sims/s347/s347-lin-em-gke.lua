-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.2

-- wave-number
knumber = 0.5
-- perpendicular wave number
kperp = 0.1
-- normalized beta (electron plasma is 2*beta*me/me)
beta = 1.0

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 16, 16

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
-- to store increment from equilibrium distribution 
incrEquil = DataStruct.Field2D {
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
   Apar1d:scale(-beta/(kperp^2-beta))
end

-- updater to compute perturbed Hamiltonian
pertHamilCalc = Updater.LinEmGke1DPertHamil {
   onGrid = grid,
   basis = basis,
   charge = -1.0, -- species charge
   mass = 1.0, -- species mass
}

-- function to compute total linearized Hamiltonian
function calcHamil(curr, dt, phi1dIn, Apar1dIn, hamilOut)
   -- calculate 2D fields from 1D fields
   copyFieldTo2D(curr, dt, phi1dIn, phi2d)
   copyFieldTo2D(curr, dt, Apar1dIn, Apar2d)

   hamilOut:clear(0.0)
   -- compute perturbed Hamiltonian
   runUpdater(pertHamilCalc, curr, dt, {phi2d, Apar2d}, {hamilOut})
   -- accumulate free-streaming contribution
   hamilOut:accumulate(1.0, hamilKE)
end

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

-- write data to H5 files
function writeFields(frame, tm)
   distf:write( string.format("distf_%d.h5", frame), tm )
   numDensity:write( string.format("numDensity_%d.h5", frame), tm )
   momentum:write( string.format("momentum_%d.h5", frame), tm )
end

-- Compute initial set of fields
calcMoments(0.0, 0.0, distf)
calcEMField(0.0, 0.0, numDensity, momentum)
calcHamil(0.0, 0.0, phi1d, Apar1d, hamil)

-- write initial conditions
writeFields(0, 0.0)