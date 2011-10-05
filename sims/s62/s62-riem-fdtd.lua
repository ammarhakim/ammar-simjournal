-- Input file to solve transverse magnetic mode Maxwell equations

cfl = 0.75
-- define some globals constants
lowerx = 0.0
upperx = 1.0
cellsx = 100
light_speed = 1.0

-- cell spacing
dx = (upperx-lowerx)/cellsx

grid = Grid.RectCart1D {
   lower = {lowerx},
   upper = {upperx},
   cells = {cellsx},
}

Elc = DataStruct.Field1D {
   -- (electric field is stored on edges)
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

ElcNew = DataStruct.Field1D {
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

Mgn = DataStruct.Field1D {
   -- (magnetic field is stored on faces)
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}

MgnNew = DataStruct.Field1D {
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}
MgnNew:clear(0.0)

-- create duplicates in case we need to take step again
ElcDup = ElcNew:duplicate()
MgnDup = MgnNew:duplicate()

function initElc(x, y, z)
   local sloc = 0.5
   local xez = x-0.5*dx -- node coordinate
   if (xez < sloc) then
      return 0.0, 1.0, 0.0
   else
      return 0.0, -1.0, 0.0
   end
end
-- initialize electric field
Elc:set(initElc)
ElcNew:copy(Elc)

function initMgn(x, y, z)
   local sloc = 0.5
   local xez = x-0.5*dx -- node coordinate
   if (xez < sloc) then
      return 1.0, -0.75, 0.0
   else
      return 1.0, 0.75, 0.0
   end
end
-- initialize electric field
Mgn:set(initMgn)
MgnNew:copy(Mgn)

-- write initial conditions
Elc:write("electricField_0.h5")

-- updater for electric field
elcUpdate = Updater.EdgeFaceCurl1D {
   onGrid = grid,
   -- ElcNew = Elc + c^2*dt*curl(Mgn)
   alpha = light_speed^2,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = cfl,
}
-- initialize updater
elcUpdate:initialize()
-- set input/out arrays (these do not change so set it once)
elcUpdate:setIn( {Elc, Mgn} )
elcUpdate:setOut( {ElcNew} )

-- updater for magnetic field
mgnUpdate = Updater.FaceEdgeCurl1D {
   onGrid = grid,
   -- MgnNew = Mgn - dt*curl(ElcNew)
   alpha = -1.0,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = cfl,
}
-- initialize updater
mgnUpdate:initialize()
-- set input/out arrays (these do not change so set it once)
mgnUpdate:setIn( {Mgn, ElcNew} )
mgnUpdate:setOut( {MgnNew} )

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)

   local step = 1
   local tCurr = tStart
   local myDt = initDt
   -- main loop
   while true do
      -- copy data in case we need to take step again
      ElcDup:copy(Elc)
      MgnDup:copy(Mgn)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- apply BCs on magnetic fields
      Mgn:applyCopyBc(0, "lower")
      Mgn:applyCopyBc(0, "upper")

      -- update electric field
      elcUpdate:setCurrTime(tCurr)
      elcStatus, elcDtSuggested = elcUpdate:advance(tCurr+myDt)

      -- apply BCs on electric field
      ElcNew:applyCopyBc(0, "lower")
      ElcNew:applyCopyBc(0, "upper")

      -- update magnetic field
      mgnUpdate:setCurrTime(tCurr)
      mgnStatus, mgnDtSuggested = mgnUpdate:advance(tCurr+myDt)

      if ( (elcStatus == false) or (mgnStatus == false) ) then
	 -- time-step too large
	 dtSuggested = math.min(elcDtSuggested, mgnDtSuggested)

	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 -- copy over data from duplicate
	 Elc:copy(ElcDup)
	 Mgn:copy(MgnDup)
      else
	 -- step succeeded, proceed to next step
	 dtSuggested = math.min(elcDtSuggested, mgnDtSuggested)

	 Elc:copy(ElcNew)
	 Mgn:copy(MgnNew)

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

-- parameters to control time-stepping
tStart = 0.0
tEnd = 0.25

-- compute time-step to move magnetic field to dt/2
dtSuggested = cfl*dx/light_speed
-- advance magnetic field by half time-step
mgnUpdate:setCurrTime(0.0)
mgnStatus, mgnDtSuggested = mgnUpdate:advance(0.5*dtSuggested)
Mgn:copy(MgnNew)

nFrames = 2
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   Elc:write( string.format("electricField_%d.h5", frame) )
   Mgn:write( string.format("magneticField_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
