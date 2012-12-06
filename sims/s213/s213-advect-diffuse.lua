-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.2/2

-- number of cells
NX, NY = 32, 1
-- extent of grid
LX, LY = 2*Lucee.Pi, 1.0
-- advection speeds
ux, uy = 1.0, 0.0
-- diffusion coefficient
alpha = 0.1

-- grid spacing
dx = LX/NX
dy = LY/NY

-- compute time-steps for hyperbolic and parabolic update
dtH = cfl*dx/math.max(ux,uy)
dtP = 0.00963829

Lucee.logInfo(string.format("Hyperbolic time-step is %g", dtH))
Lucee.logInfo(string.format("Parabolic  time-step is %g\n", dtP))

-- grid on which equations are to be solved
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {LX, LY},
   cells = {NX, NY},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, 2 or 3. Corresponding
   -- number of nodes are 2, 3, or 4.
   polyOrder = polyOrder,
}

-- number of nodes per cell for DG fields
numDgNodesPerCell = basis:numNodes()

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}

-- for RK time-stepping
q1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}
qDiv = DataStruct.Field2D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}

-- duplicate for use in time-stepping
qNewDup = qNew:duplicate()

-- to store gradients
gradQ = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 2*numDgNodesPerCell,
   ghost = {2, 2},
}

-- updater to apply initial conditions
initField = Updater.EvalOnNodes2D {
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

-- define equation to solve
advectionEqn = HyperEquation.Advection {
   -- advection velocity
   speeds = {1.0, 0.0, 0.0}
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

-- gradient equation
gradEqn = HyperEquation.GradAuxFlux2D {
   -- coefficient for gradient
   coefficient = alpha,
   -- flux to use (one of 3-point or 5-point)
   fluxType = "3-point",
}
-- divergence equation
divEqn = HyperEquation.DivAuxFlux2D {
   -- flux to use (one of 3-point or 5-point)
   fluxType = "3-point",
}

-- updater to compute gradients
gradSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- update-type
   onlyIncrement = true,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = gradEqn,
   -- CFL number
   cfl = cfl,
}

-- updater to compute divergence
divSlvr = Updater.NodalDgHyper2D {
   onGrid = grid,
   -- update-type
   onlyIncrement = true,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = divEqn,
   -- CFL number
   cfl = cfl,
}

-- apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
end

applyBc(q)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- solve advection equation
function solveAdvection(curr, dt, qIn, qOut)
   advectSlvr:setCurrTime(curr)
   advectSlvr:setIn( {qIn} )
   advectSlvr:setOut( {qOut} )
   return advectSlvr:advance(curr+dt)
end

-- compute gradients
function calcGradient(curr, dt, qIn, gradOut)
   gradSlvr:setCurrTime(curr)
   gradSlvr:setIn( {gradOut, qIn} )
   gradSlvr:setOut( {gradOut} )
   return gradSlvr:advance(curr+dt)
end

-- compute divergence
function calcDivergence(curr, dt, gradIn, divOut)
   divSlvr:setCurrTime(curr)
   divSlvr:setIn( {divOut, gradIn} )
   divSlvr:setOut( {divOut} )
   return divSlvr:advance(curr+dt)
end

-- compute diffusion term contribution
function calcDiffusion(curr, dt, qIn, qDiffOut)
   if (dt > dtP) then
      -- time-step exceeds CFL from diffusion
      return false, 0.99*dtP
   end

   -- compute gradient of field
   calcGradient(curr, dt, qIn, gradQ)
   applyBc(gradQ)
   -- compute divergence
   calcDivergence(curr, dt, gradQ, qDiffOut)
   applyBc(qDiffOut)

   return true, dtP
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q, q1)
   -- compute diffusion term
   local myDiffStatus, myDiffDtSuggested = calcDiffusion(tCurr, myDt, q, qDiv)
   -- accumulate into solution
   q1:accumulate(myDt, qDiv)

   if (myStatus == false or myDiffStatus == false) then
      return myStatus, math.min(myDtSuggested, myDiffDtSuggested)
   end

   applyBc(q1)

   -- RK stage 2
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q1, qNew)
   -- compute diffusion term
   local myDiffStatus, myDiffDtSuggested = calcDiffusion(tCurr, myDt, q1, qDiv)
   -- accumulate into solution
   qNew:accumulate(myDt, qDiv)

   if (myStatus == false or myDiffStatus == false) then
      return myStatus, math.min(myDtSuggested, myDiffDtSuggested)
   end

   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1)

   -- RK stage 3
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q1, qNew)
   -- compute diffusion term
   local myDiffStatus, myDiffDtSuggested = calcDiffusion(tCurr, myDt, q1, qDiv)
   -- accumulate into solution
   qNew:accumulate(myDt, qDiv)

   if (myStatus == false or myDiffStatus == false) then
      return myStatus, math.min(myDtSuggested, myDiffDtSuggested)
   end

   q1:combine(1.0/3.0, q, 2.0/3.0, qNew)
   applyBc(q1)
   q:copy(q1)

   return myStatus, math.min(myDtSuggested, myDiffDtSuggested)
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

      Lucee.logDebug (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 Lucee.logDebug (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 qNew:copy(qNewDup)
	 myDt = dtSuggested
      else
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
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 2*Lucee.Pi
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 4
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
