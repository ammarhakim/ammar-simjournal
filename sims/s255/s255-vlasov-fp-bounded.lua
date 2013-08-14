-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.2

-- domain extents
XL, XU = 0.0, 10.0
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 32, 16

-- phase space grid 
grid = Grid.RectCart2D {
   lower = {XL, VL},
   upper = {XU, VU},
   cells = {NX, NV},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 4 and 8.
   polyOrder = polyOrder,
}

-- distribution function
distf = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
distf:clear(0.0)

-- extra fields for performing RK update
distfNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*basis:numNodes(),
   ghost = {1, 1},
}
distf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*basis:numNodes(),
   ghost = {1, 1},
}

-- Hamiltonian
hamil = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*basis:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- updater to initialize hamiltonian
initHamil = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local phi = math.tanh((10-x))  --x^10/5*math.exp(x^2)/math.exp(1.0)
		 return y^2/2 + phi/1836.2
	      end
}
initHamil:setOut( {hamil} )
-- initialize potential
initHamil:advance(0.0) -- time is irrelevant

-- updater to initialize distribution function
initDistf = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local v, vt = y, 1.0
		 local vdrift = 0.0
		 return 1/math.sqrt(2*Lucee.Pi*vt)*math.exp(-(v-vdrift)^2/(2*vt^2))
	      end
}
initDistf:setOut( {distf} )
-- initialize potential
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

-- number density
numDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = 1*basis_1d:numNodes(),
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
-- output is number density
numDensityCalc:setOut( {numDensity} )

-- dynvector for total particle count
totalPtcl = DataStruct.DynVector { numComponents = 1, }

-- to compute total number of particles in domain
totalPtclCalc = Updater.IntegrateNodalField1D {
   onGrid = grid_1d,
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclCalc:setIn( {numDensity} )
-- set output dynvector
totalPtclCalc:setOut( {totalPtcl} )

-- dynvector for number density in a cell
numDensInCell = DataStruct.DynVector { 
   numComponents = basis_1d:numNodes(),
}

-- to compute number density in a cell
numDensInCellCalc = Updater.RecordFieldInCell1D {
   onGrid = grid_1d,
   cellIndex = {2},
}
-- set input field
numDensInCellCalc:setIn( {numDensity} )
-- set output dynvector
numDensInCellCalc:setOut( {numDensInCell} )

-- ptcl energy
ptclEnergy = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = 1*basis_1d:numNodes(),
   ghost = {1, 1},
}

-- updater to compute ptcl energy
ptclEnergyCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 2,
}
-- output is ptcl energy
ptclEnergyCalc:setOut( {ptclEnergy} )

-- dynvector for total ptcl energy
totalPtclEnergy = DataStruct.DynVector { numComponents = 1, }

-- to compute total particle energy
totalPtclEnergyCalc = Updater.IntegrateNodalField1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclEnergyCalc:setIn( {ptclEnergy} )
-- set output dynvector
totalPtclEnergyCalc:setOut( {totalPtclEnergy} )

function applyFld1D(fld)
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
end

bcLeftWall = BoundaryCondition.NodalDgFunction2D {
   basis = basis,
   components = {0},
   bc = function(x,y,z,t)
	   local v, vt = y, 1.0
	   local vdrift = 0.0
	   return 1/math.sqrt(2*Lucee.Pi*vt)*math.exp(-(v-vdrift)^2/(2*vt^2))
	end,
}
bcRightWall = BoundaryCondition.NodalDgFunction2D {
   basis = basis,
   components = {0},
   bc = function(x,y,z,t)
	   return 0.0
	end,
}
bcLower = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcLeftWall},
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "lower",
}
bcUpper = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcRightWall},
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "upper",
}

function applyFld2D(curr, dt, fld)
   for i,bc in ipairs({bcLower, bcUpper}) do
      bc:setCurrTime(curr)
      bc:setOut( {fld} )
      bc:advance(curr+dt)
   end
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:advance(curr+dt)

   ptclEnergyCalc:setCurrTime(curr)
   ptclEnergyCalc:setIn( {distfIn} )
   ptclEnergyCalc:advance(curr+dt)
end

-- compute initial moments
calcMoments(0.0, 0.0, distf)

-- function to apply boundary conditions
function applyBc(fld)
   applyFld2D(0.0, 0.0, fld)
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

-- apply BCs to initial conditions
applyBc(distf)

-- write initial conditions
hamil:write("hamil_0.h5")
distf:write("distf_0.h5")
numDensity:write("numDensity_0.h5")

-- update Poisson bracket operator
function poissonBracket(curr, dt, distfIn, hamilIn, distfOut)
   pbSlvr:setCurrTime(curr)
   pbSlvr:setIn( {distfIn, hamilIn} )
   pbSlvr:setOut( {distfOut} )
   return pbSlvr:advance(curr+dt)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   totalPtclCalc:setCurrTime(curr)
   totalPtclCalc:advance(curr+dt)

   numDensInCellCalc:setCurrTime(curr)
   numDensInCellCalc:advance(curr+dt)

   totalPtclEnergyCalc:setCurrTime(curr)
   totalPtclEnergyCalc:advance(curr+dt)
end

-- compute initial diagnostics
calcDiagnostics(0.0, 0.0) -- funky?

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1
   status, dtSuggested = poissonBracket(tCurr, myDt, distf, hamil, distf1)

   if (status == false) then
      return status, dtSuggested
   end

   applyBc(distf1)
   calcMoments(tCurr, myDt, distf1)

   -- RK stage 2
   status, dtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)

   if (status == false) then
      return status, dtSuggested
   end

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(distf1)
   calcMoments(tCurr, myDt, distf1)

   -- RK stage 3
   status, dtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)

   if (status == false) then
      return status, dtSuggested
   end

   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)

   applyBc(distf1)
   distf:copy(distf1)
   calcMoments(tCurr, myDt, distf)

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

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
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

-- write data to H5 files
function writeFields(frame, tm)
   distf:write( string.format("distf_%d.h5", frame), tm )
   numDensity:write( string.format("numDensity_%d.h5", frame), tm )
   totalPtcl:write(string.format("totalPtcl_%d.h5", frame), tm )
   totalPtclEnergy:write(string.format("totalPtclEnergy_%d.h5", frame), tm )
   numDensInCell:write(string.format("numDensInCell_%d.h5", frame), tm )
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 100.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 10
tFrame = (tEnd-tStart)/nFrames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
