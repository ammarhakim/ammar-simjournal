-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.25/(2*polyOrder+1)

-- grid on which equations are to be solved
grid = Grid.RectCart1D {
   lower = {0},
   upper = {2*Lucee.Pi},
   cells = {64},
}

-- create FEM nodal basis
basis = NodalFiniteElement1D.Lobatto {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}
-- for RK time-stepping
q1 = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}

-- duplicate for use in time-stepping
qNewDup = qNew:duplicate()

-- updater to apply initial conditions
initField = Updater.EvalOnNodes1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return math.sin(x)
	      end
}
initField:setOut( {q} )
-- initialize
initField:advance(0.0) -- time is irrelevant

-- updater to solve diffusion equation
diffSolver = Updater.HyperDiffusion1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- diffusion coefficent
   diffusionCoeff = 1.0,
   -- CFL number
   cfl = cfl,
}

-- total enstrophy diagnostic
totalEnstrophy = DataStruct.DynVector {
   -- number of components in diagnostic
   numComponents = 1,
}

-- updater to compute total energy
enstrophyCalc = Updater.TotalEnstrophy1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
}
-- set input/output (this never changes, so it once)
enstrophyCalc:setIn( {q} )
enstrophyCalc:setOut( {totalEnstrophy} )

-- compute initial enstrophy of system
enstrophyCalc:advance(0)

-- apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
end

applyBc(q)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- solve advection equation
function solveDiffusion(curr, dt, qIn, qOut)
   diffSolver:setCurrTime(curr)
   diffSolver:setIn( {qIn} )
   diffSolver:setOut( {qOut} )
   return diffSolver:advance(curr+dt)
end

function calcDiagnostics(curr, dt)
   enstrophyCalc:setCurrTime(curr)
   enstrophyCalc:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = solveDiffusion(tCurr, myDt, q, q1)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   applyBc(q1)

   -- RK stage 2
   local myStatus, myDtSuggested = solveDiffusion(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end
   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1)

   -- RK stage 3
   local myStatus, myDtSuggested = solveDiffusion(tCurr, myDt, q1, qNew)
   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   q1:combine(1.0/3.0, q, 2.0/3.0, qNew)
   applyBc(q1)
   q:copy(q1)

   return myStatus, myDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      qNewDup:copy(qNew)
      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 qNew:copy(qNewDup)
	 myDt = dtSuggested
      else
	 -- compute diagnostics
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

-- write data to H5 file
function writeFields(frame)
   q:write( string.format("q_%d.h5", frame) )
   totalEnstrophy:write( string.format("totalEnstrophy_%d.h5", frame) )
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
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

