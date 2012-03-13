-- Input file to solve transverse magnetic mode Maxwell
-- equations. This uses a dual of the standard Yee cell and stores the
-- magentic field on edges and electric field on faces.

-- define some globals constants
lowerx = -1.0
lowery = -1.0
upperx = 1.0
uppery = 1.0
cellsx = 100
cellsy = 100
light_speed = 1.0

-- cell spacing
dx = (upperx-lowerx)/cellsx
dy = (uppery-lowery)/cellsy

grid = Grid.RectCart2D {
   lower = {lowerx, lowery},
   upper = {upperx, uppery},
   cells = {cellsx, cellsy},
}

Elc = DataStruct.Field2D {
   -- (electric field is stored on faces)
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

ElcNew = DataStruct.Field2D {
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}

Mgn = DataStruct.Field2D {
   -- (magnetic field is stored on edges)
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}
Mgn:clear(0.0)

MgnNew = DataStruct.Field2D {
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}
MgnNew:clear(0.0)

-- create duplicates in case we need to take step again
ElcDup = ElcNew:duplicate()
MgnDup = MgnNew:duplicate()

function initEz(x, y, z)
   local xc, yc = 0.0, 0.0
   local rad2 = (x-xc)^2 + (y-yc)^2
   return 0.0, 0.0, math.exp(-25*rad2)
end
-- initialize electric field
Elc:set(initEz)
ElcNew:copy(Elc)

-- write initial conditions
Mgn:write("magneticField_0.h5")
Elc:write("electricField_0.h5")

-- updater for magnetic field
mgnUpdate = Updater.EdgeFaceCurl2D {
   onGrid = grid,
   -- MgnNew = Mgn - dt*curl(ElcNew)
   alpha = -1.0,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = 0.45,
   -- extra cells to update (needed for the dual Yee cell)
   ghostUpdate = {0, 1},
}
-- set input/out arrays (these do not change so set it once)
mgnUpdate:setIn( {Mgn, Elc} )
mgnUpdate:setOut( {MgnNew} )

-- updater for electric field
elcUpdate = Updater.FaceEdgeCurl2D {
   onGrid = grid,
   -- ElcNew = Elc + c^2*dt*curl(Mgn)
   alpha = light_speed^2,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = 0.45,
}
-- set input/out arrays (these do not change so set it once)
elcUpdate:setIn( {Elc, MgnNew} )
elcUpdate:setOut( {ElcNew} )

-- for use in electric field Bcs
bcCopyFlip = BoundaryCondition.Copy { components = {0, 1, 2}, fact = {-1, -1, -1} }
-- for use in magnetic field BCs
bxCopyFlip = BoundaryCondition.Copy { components = {0}, fact = {-1} }
byCopyFlip = BoundaryCondition.Copy { components = {1}, fact = {-1} }

-- function to create boundary condition object
function createBc(myDir, myEdge, outFld, appBc)
   local bc = Updater.Bc2D {
      onGrid = grid,
      -- boundary conditions to apply
      boundaryConditions = {appBc},
      -- direction to apply
      dir = myDir,
      -- edge to apply on
      edge = myEdge,
   }
   bc:setOut( {outFld} )
   return bc
end

-- create boundary conditions updaters
bcElcLeft = createBc(0, "lower", Elc, bcCopyFlip)
bcElcRight = createBc(0, "upper", Elc, bcCopyFlip)
bcElcTop = createBc(1, "upper", Elc, bcCopyFlip)
bcElcBottom = createBc(1, "lower", Elc, bcCopyFlip)

bcMgnLeft = createBc(0, "lower", MgnNew, bxCopyFlip)
bcMgnRight = createBc(0, "upper", MgnNew, bxCopyFlip)
bcMgnTop = createBc(1, "upper", MgnNew, byCopyFlip)
bcMgnBottom = createBc(1, "lower", MgnNew, byCopyFlip)

-- run an updater
function runUpdater(updater, tcurr, dt)
   updater:setCurrTime(tcurr)
   updater:advance(tcurr+dt)
end

-- take one time-step, without worrying about adapting it
function singleStep(tCurr, myDt)
   -- apply BCs on electric field
   runUpdater(bcElcLeft, tCurr, myDt)
   runUpdater(bcElcRight, tCurr, myDt)
   runUpdater(bcElcBottom, tCurr, myDt)
   runUpdater(bcElcTop, tCurr, myDt)

   -- update magnetic field (Mgn -> MgnNew)
   mgnUpdate:setCurrTime(tCurr)
   mgnStatus, mgnDtSuggested = mgnUpdate:advance(tCurr+myDt)

   -- apply BCs on magnetic fields
   runUpdater(bcMgnLeft, tCurr, myDt)
   runUpdater(bcMgnRight, tCurr, myDt)
   runUpdater(bcMgnTop, tCurr, myDt)
   runUpdater(bcMgnBottom, tCurr, myDt)

   -- update electric field (Elc -> ElcNew)
   elcUpdate:setCurrTime(tCurr)
   elcStatus, elcDtSuggested = elcUpdate:advance(tCurr+myDt)

   return elcStatus and mgnStatus, math.min(elcDtSuggested, mgnDtSuggested)
end

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

      -- take a time-step
      status, dtSuggested = singleStep(tCurr, myDt)

      if ( status == false ) then
	 -- time-step too large

	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 -- copy over data from duplicate
	 Elc:copy(ElcDup)
	 Mgn:copy(MgnDup)
      else
	 -- step succeeded, proceed to next step

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
tEnd = 3.0

-- compute time-step to move magnetic field to dt/2
dtSuggested = math.min(0.45*dx/light_speed, 0.45*dy/light_speed)
dt2 = 0.5*dtSuggested
-- move the magnetic field by -dt/2
runUpdater(bcElcLeft, tStart, -dt2)
runUpdater(bcElcRight, tStart, -dt2)
runUpdater(bcElcBottom, tStart, -dt2)
runUpdater(bcElcTop, tStart, -dt2)
mgnUpdate:setCurrTime(tStart)
mgnStatus, mgnDtSuggested = mgnUpdate:advance(tStart-dt2)
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
