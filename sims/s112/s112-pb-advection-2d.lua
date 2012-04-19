-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.025

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
end

-- Determine number of global nodes per cell for use in creating DG
-- fields.
if (polyOrder == 1) then
   numDgNodesPerCell = 4
elseif (polyOrder == 2) then
   numDgNodesPerCell = 8
end

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

-- potential
phi = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create updater to initialize potential
initPhi = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local ux, uy = 1.0, 1.0
		 return ux*y - uy*x
	      end
}
initPhi:setOut( {phi} )
-- initialize potential
initPhi:advance(0.0) -- time is irrelevant

-- create updater to initialize vorticity
initChi = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local xc, yc = 0.5, 0.5
		 local r2 = (x-xc)^2 + (y-yc)^2
		 return math.exp(-75*r2)
	      end
}
initChi:setOut( {chi} )
-- initialize potential
initChi:advance(0.0) -- time is irrelevant

-- apply BC to get ghost correct
chiNew:applyPeriodicBc(0)
chiNew:applyPeriodicBc(1)

-- create updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}

-- write initial value
chi:write("chi_0.h5")

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
end

function poissonBracket(curr, dt, chiIn, phiIn, chiOut)
   pbSlvr:setCurrTime(curr)
   pbSlvr:setIn( {chiIn, phiIn} )
   pbSlvr:setOut( {chiOut} )
   return pbSlvr:advance(curr+dt)
end

-- function to take a time-step using RK2 time-stepping scheme
function rk2(tCurr, myDt)

   local status, dtSuggested

   -- RK stage 1 (chi1 <- chi + L(chi))
   status, dtSuggested = poissonBracket(tCurr, myDt, chi, phi, chi1)   

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   -- apply periodic BC
   applyBc(chi1)

   -- RK stage 2 (chiNew <- chi1 + L(chi1))
   status, dtSuggested = poissonBracket(tCurr, myDt, chi1, phi, chiNew)   

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   -- do final update of chi1 <-- 0.5*(chi + chiNew)
   chi1:combine(0.5, chi, 0.5, chiNew)
   -- apply Bcs
   applyBc(chi1)

   -- copy over solution
   chi:copy(chi1)

   return status, dtSuggested
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1 (chi1 <- chi + L(chi))
   status, dtSuggested = poissonBracket(tCurr, myDt, chi, phi, chi1)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   -- apply BCs
   applyBc(chi1)

   -- RK stage 2
   status, dtSuggested = poissonBracket(tCurr, myDt, chi1, phi, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)

   -- apply BCs
   applyBc(chi1)

   -- RK stage 3
   status, dtSuggested = poissonBracket(tCurr, myDt, chi1, phi, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   -- apply BCs
   applyBc(chi1)

   -- copy over solution
   chi:copy(chi1)

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
      status, dtSuggested = rk2(tCurr, myDt)

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
tEnd = 1.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1
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
