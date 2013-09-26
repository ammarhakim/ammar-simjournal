-- Program to solve Maxwell equations

-- CFL number
cfl = 0.75

-- computational domain
grid = Grid.RectCart2D {
   lower = {-1.0, -1.0},
   upper = {1.0, 1.0},
   cells = {151, 151},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   -- [Ex, Ey, Ez, Bx, By, Bz, phi_e, phi_m]
   numComponents = 8,
   ghost = {2, 2},
}
qX = DataStruct.Field2D { -- to store results from X sweep
   onGrid = grid,
   -- [Ex, Ey, Ez, Bx, By, Bz, phi_e, phi_m]
   numComponents = 8,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   -- [Ex, Ey, Ez, Bx, By, Bz, phi_e, phi_m]
   numComponents = 8,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- create an alias to in-plane electric field
eFldAlias = qNew:alias(0, 2) -- [Ex, Ey]

-- initial condition to apply
function init(x,y,z)
   local xc, yc = 0.0, 0.0
   local rad2 = (x-xc)^2 + (y-yc)^2
   return 0.0, 0.0, math.exp(-25*rad2), 0.0, 0.0, 0.0, 0.0, 0.0
end

-- apply initial conditions
q:set(init)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5", 0.0)

-- define equation to solve
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = 1.0,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 1.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 1.0,
}


-- updater for Euler equations
maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "no-limiter", 
   cfl = cfl,
   cflm = cfl*1.01,
   updateDirections = {0} -- directions to update
}
-- set input/output arrays (these do not change so set it once)
maxSlvrDir0:setIn( {q} )
maxSlvrDir0:setOut( {qX} )

maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "no-limiter", 
   cfl = cfl,
   cflm = cfl*1.01,
   updateDirections = {1} -- directions to update
}
maxSlvrDir1:setIn( {qX} )
maxSlvrDir1:setOut( {qNew} )

-- boundary condition to apply
bcElc = BoundaryCondition.ZeroTangent { components = {0, 1, 2} }
bcMgn = BoundaryCondition.ZeroNormal { components = {3, 4, 5} }
potBc = BoundaryCondition.Copy { components = {6, 7}, fact = {-1, -1} }

-- create boundary condition object
function createBc(myDir, myEdge)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {bcElc, bcMgn, potBc},
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   return bc
end

-- create updaters to apply boundary conditions
bcLeft = createBc(0, "lower")
bcRight = createBc(0, "upper")
bcBottom = createBc(1, "lower")
bcTop = createBc(1, "upper")

-- apply boundary conditions
function applyBcDir(dir, fld, tCurr, dt)
   if dir == 0 then
      for i,bc in ipairs({bcLeft, bcRight}) do
	 bc:setOut( {fld} )
	 bc:setCurrTime(tCurr)
	 bc:advance(tCurr+dt)
      end
   else
      for i,bc in ipairs({bcTop, bcBottom}) do
	 bc:setOut( {fld} )
	 bc:setCurrTime(tCurr)
	 bc:advance(tCurr+dt)
      end
   end
end

function solveMaxwell(tCurr, dt)
   -- X-direction update
   maxSlvrDir0:setCurrTime(tCurr)
   local s0, dt0 = maxSlvrDir0:advance(tCurr+dt)
   applyBcDir(0, qX, tCurr, dt)

   -- Y-direction update
   maxSlvrDir1:setCurrTime(tCurr)
   local s1, dt1 = maxSlvrDir1:advance(tCurr+dt)
   applyBcDir(1, qNew, tCurr, dt)

   return s1 and s2, math.min(dt0, dt1)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local status, dtSuggested
   local tCurr = tStart
   local myDt = initDt
   while true do
      qNewDup:copy(qNew)
      if (tCurr+myDt > tEnd) then -- hit tEnd exactly
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("  Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = solveMaxwell(tCurr, myDt)

      if (status == false) then
	 Lucee.logInfo (string.format("  ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
      else
	 q:copy(qNew)
	 tCurr = tCurr + myDt
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end
   end

   return dtSuggested
end

dtSuggested = 100.0 -- initial suggested time-step
-- parameters to control time-stepping
tStart = 0.0
tEnd = 3.0

nFrames = 60
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   q:write( string.format("q_%d.h5", frame), tCurr+tFrame )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

