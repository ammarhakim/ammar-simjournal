-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1/(2*polyOrder+1)

LX = 2*Lucee.Pi -- domain length
NX = 8 -- number of cells
dx = LX/NX -- cell size

-- grid on which equations are to be solved
grid = Grid.RectCart1D {
   lower = {0},
   upper = {LX},
   cells = {NX},
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = polyOrder+1,
   ghost = {2, 2},
}
-- for RK time-stepping
q1 = DataStruct.Field1D {
   onGrid = grid,
   numComponents = polyOrder+1,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = polyOrder+1,
   ghost = {2, 2},
}
-- heat source
src = DataStruct.Field1D {
   onGrid = grid,
   numComponents = polyOrder+1,
   ghost = {2, 2},
}

-- duplicate for use in time-stepping
qNewDup = qNew:duplicate()

-- updater to apply initial conditions
initField = Updater.ProjectOnBasis1D {
   onGrid = grid,
   numBasis = polyOrder+1,
   project = function (x,y,z,t)
		return 0.0
	     end
}
initField:setOut( {q} )
-- initialize
initField:advance(0.0) -- time is irrelevant

-- updater to apply initial conditions
initSrc = Updater.ProjectOnBasis1D {
   onGrid = grid,
   numBasis = polyOrder+1,
   project = function (x,y,z,t)
		return math.sin(x)
	      end
}
initSrc:setOut( {src} )
-- initialize
---initSrc:advance(0.0) -- time is irrelevant

-- function to initialize source with an exact projection
function initSrcExact(x,y,z)
   local xj = x
   local f0 = (math.cos((2*xj-dx)/2.0)-math.cos((2*xj+dx)/2.0))/dx
   local f1 = ((2*math.sin((2*xj+dx)/2)-dx*math.cos((2*xj+dx)/2)-2*math.sin((2*xj-dx)/2)-dx*math.cos((2*xj-dx)/2))/dx)*3.0/dx
   return f0, f1
end
src:set(initSrcExact)
src:write("src.h5")

-- updater to solve diffusion equation
diffSolver = Updater.DGDiffusion1D { 
   onGrid = grid,
   -- polynomial order
   polyOrder = 1,
   -- scheme
   scheme = "LDG-R",
   -- diffusion coefficent
   diffusionCoeff = 1.0,
   -- CFL number
   cfl = cfl,
}

-- L2 norm of solution (energy)
l2Norm = DataStruct.DynVector {
   -- number of components in diagnostic
   numComponents = 1,
}

l2NormCalc = Updater.ModalL2Norm {
   onGrid = grid,
   -- number of basis function
   numBasis = polyOrder+1,
}
l2NormCalc:setIn( {q} )
l2NormCalc:setOut( {l2Norm} )

-- compute L2 norm of initial conditions
l2NormCalc:advance(0.0)

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
   local myS, myDt = diffSolver:advance(curr+dt)
   -- accumulate source
   qOut:accumulate(dt, src)

   return myS, myDt
end

function calcDiagnostics(curr, dt)
   l2NormCalc:setCurrTime(curr)
   l2NormCalc:advance(curr+dt)
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
      if (tCurr+myDt > tEnd) then -- adjust to tEnd exactly
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 qNew:copy(qNewDup)
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

-- write data to H5 file
function writeFields(frame)
   q:write( string.format("q_%d.h5", frame) )
   l2Norm:write( string.format("l2Norm_%d.h5", frame) )
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 10.0
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

