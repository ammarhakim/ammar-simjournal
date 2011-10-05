-- Program to solve Maxwell equations

-- CFL number
cfl = 0.9

-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {100},
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   -- [Ex, Ey, Ez, Bx, By, Bz, phi_e, phi_m]
   numComponents = 8,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
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
   local sloc = 0.5
   if (x<sloc) then
      return 0.0, 1.0, 0.0, 1.0, -0.75, 0.0, 0.0, 0.0
   else
      return 0.0, -1.0, 0.0, 1.0, 0.75, 0.0, 0.0, 0.0
   end
end

-- apply initial conditions
q:set(init)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

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
maxSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered", 
   cfl = cfl,
   cflm = cfl*1.01
}

-- initialize updater
maxSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
maxSlvr:setIn( {q} )
maxSlvr:setOut( {qNew} )

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)

   local step = 1
   local tCurr = tStart
   local myDt = initDt
   while true do
      -- copy qNew in case we need to take this step again
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (
	 string.format("  Taking step %d at time %g with dt %g", 
		       step, tCurr, myDt))

      -- set current time
      maxSlvr:setCurrTime(tCurr)
      -- advance solution
      status, dtSuggested = maxSlvr:advance(tCurr+myDt)

      if (dtSuggested < myDt) then
	 -- time-step too large
	 Lucee.logInfo (
	    string.format("  ** Time step %g too large! Will retake with dt %g", 
			  myDt, dtSuggested))
	 myDt = dtSuggested
	 qNew:copy(qNewDup)
      else
	 -- apply copy BCs on lower and upper edges
	 qNew:applyCopyBc(0, "lower")
	 qNew:applyCopyBc(0, "upper")

	 -- copy updated solution back
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
tEnd = 0.25

nFrames = 2
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   q:write( string.format("q_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

