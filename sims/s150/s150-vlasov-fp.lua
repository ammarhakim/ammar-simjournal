-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.2

-- domain extents
XL, XU = -1.0, 1.0
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 64, 128

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
   numCgNodesPerCell_1d = 1 -- for spatial basis
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
   numCgNodesPerCell_1d = 2 -- for spatial basis
end

-- Determine number of global nodes per cell for use in creating DG
-- fields.
if (polyOrder == 1) then
   numDgNodesPerCell = 4
   numDgNodesPerCell_1d = 2 -- for spatial basis
elseif (polyOrder == 2) then
   numDgNodesPerCell = 8
   numDgNodesPerCell_1d = 3 -- for spatial basis
end

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
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
-- clear out contents
distf:clear(0.0)

-- extra fields for performing RK update
distfNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
distf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- Hamiltonian
hamil = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
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
      local v = y
      return v^2/2 + x^2
   end
}
initHamil:setOut( {hamil} )
-- initialize potential
initHamil:advance(0.0) -- time is irrelevant
hamil:applyPeriodicBc(0)

function uniformMaxwellian(x,y,z,t)
   local v, vt = y, 1.0
   return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-v^2/(2*vt^2))
end

function circle(x,y,z,t)
   local v = y
   local xc, vc = Lucee.Pi, 1.0
   local rad = 0.25
   local r = math.sqrt((x-xc)^2 + (v-vc)^2)
   if r < rad then
      return 1.0
   else
      return 0.0
   end
end

-- updater to initialize distribution function
initDistf = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
      return uniformMaxwellian(x,y,z,t)
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
   -- grid on which elements should be constructured
   onGrid = grid_1d,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 2 and 3.
   polyOrder = polyOrder,
}

-- number density
numDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numDgNodesPerCell_1d,
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
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclCalc:setIn( {numDensity} )
-- set output dynvector
totalPtclCalc:setOut( {totalPtcl} )

-- ptcl energy
ptclEnergy = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numDgNodesPerCell_1d,
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

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:advance(curr+dt)

   numDensity:applyPeriodicBc(0)

   ptclEnergyCalc:setCurrTime(curr)
   ptclEnergyCalc:setIn( {distfIn} )
   ptclEnergyCalc:advance(curr+dt)

   ptclEnergy:applyPeriodicBc(0)
end

-- compute initial moments
calcMoments(0.0, 0.0, distf)

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

-- apply BCs to initial conditions
applyBc(distf)

-- write initial conditions
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
function writeFields(frame)
   distf:write( string.format("distf_%d.h5", frame) )
   numDensity:write( string.format("numDensity_%d.h5", frame) )
   totalPtcl:write(string.format("totalPtcl_%d.h5", frame) )
   totalPtclEnergy:write(string.format("totalPtclEnergy_%d.h5", frame) )
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 20.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 100
tFrame = (tEnd-tStart)/nFrames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
