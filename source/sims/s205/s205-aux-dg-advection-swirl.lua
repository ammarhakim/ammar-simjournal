-- Input file for advection with auxiliary variables

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.5/(2*polyOrder-1)

-- period for distortion
Tperiod = 1.5

-- grid on which equations are to be solved
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {1.0, 1.0},
   cells = {32, 32},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 4 and 8.
   polyOrder = polyOrder,
}

-- number of nodes per cell for CG field
numCgNodesPerCell = basis:numExclusiveNodes()
-- number of nodes per cell for DG field
numDgNodesPerCell = basis:numNodes()

-- vorticity
chi = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
-- clear out contents
chi:clear(0.0)

-- extra fields for performing RK update
chiNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}
chi1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1*numDgNodesPerCell,
   ghost = {1, 1},
}

-- velocity field
flowFieldStatic = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 2*numDgNodesPerCell, -- [ux, uy]
   ghost = {2, 2},
}
flowField = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 2*numDgNodesPerCell, -- [ux, uy]
   ghost = {2, 2},
}

-- create updater to initialize potential
initFlowField = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local mpi = Lucee.Pi
		 local ux = math.sin(mpi*x)^2*math.sin(2*mpi*y)
		 local uy = -math.sin(mpi*y)^2*math.sin(2*mpi*x)
		 return ux, uy
	      end
}
initFlowField:setOut( {flowFieldStatic} )
-- initialize potential
initFlowField:advance(0.0) -- time is irrelevant

-- write it out
flowFieldStatic:write("flowField_0.h5")

-- create updater to initialize vorticity
initChi = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local x0, y0, r0 = 0.25, 0.5, 0.15
		 local r = math.min(math.sqrt((x-x0)^2+(y-y0)^2), r0)/r0
		 local chump = 0.25*(1+math.cos(Lucee.Pi*r))

		 local xc, yc, rc = 0.75, 0.5, 0.15
		 local r2 = math.min(math.sqrt((x-xc)^2+(y-yc)^2), rc)/rc
		 local chump2 = 0.25*(1+math.cos(Lucee.Pi*r2))

		 return chump + chump2		 
	      end
}
initChi:setOut( {chi} )
-- initialize potential
initChi:advance(0.0) -- time is irrelevant

-- write initial value
chi:write("chi_0.h5")

-- define equation to solve
advectionEqn = HyperEquation.AuxAdvection2D {
}

-- updater to solve hyperbolic equations
advectSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = advectionEqn,
   -- CFL number
   cfl = cfl,
}

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

-- apply BCs to initial conditions
applyBc(chi)
chiNew:copy(chi)

function solveAdvection(curr, dt, chiIn, flowIn, chiOut)
   advectSlvr:setCurrTime(curr)
   advectSlvr:setIn( {chiIn, flowIn} ) -- flow-field is specified
   advectSlvr:setOut( {chiOut} )
   return advectSlvr:advance(curr+dt)
end

-- function to solve for potential
function poisson(curr, dt, chiIn, flowFieldOut)
   -- this is a time-dependent potential and not a true Poisson solve
   flowFieldOut:combine(math.cos(curr*Lucee.Pi/Tperiod), flowFieldStatic)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1 (chi1 <- chi + L(chi))
   status, dtSuggested = solveAdvection(tCurr, myDt, chi, flowField, chi1)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   -- apply BCs
   applyBc(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi1, flowField)

   -- RK stage 2
   status, dtSuggested = solveAdvection(tCurr, myDt, chi1, flowField, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)

   -- apply BCs
   applyBc(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi1, flowField)

   -- RK stage 3
   status, dtSuggested = solveAdvection(tCurr, myDt, chi1, flowField, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   -- apply BCs
   applyBc(chi1)

   -- copy over solution
   chi:copy(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi1, flowField)

   return status, dtSuggested
end

-- make a duplicate in case we need it
chiDup = chi:duplicate()

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   -- declare local variables
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- main loop
   while tCurr<=tEnd do
      -- copy chi over
      chiDup:copy(chi)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- take a time-step
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 -- copy in case current solutions were messed up
	 chi:copy(chiDup)
	 myDt = dtSuggested
      else
	 tCurr = tCurr + myDt
	 myDt = dtSuggested
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end

   end

   return dtSuggested
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 2*Tperiod
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 40
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   chi:write( string.format("chi_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
