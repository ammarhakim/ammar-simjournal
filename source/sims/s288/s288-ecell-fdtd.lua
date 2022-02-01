-- Dual Yee FDTD scheme, with electric field at faces as well as at
-- cell-center.

-- define some globals constants
lowerx = 0.0
lowery = 0.0
upperx = 80.0
uppery = 40.0
cellsx = 80
cellsy = 40
light_speed = Lucee.SpeedOfLight

-- cell spacing
dx = (upperx-lowerx)/cellsx
dy = (uppery-lowery)/cellsy

grid = Grid.RectCart2D {
   lower = {lowerx, lowery},
   upper = {upperx, uppery},
   cells = {cellsx, cellsy},
}

ElcFace = DataStruct.Field2D {
   -- (electric field stored on faces)
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}
ElcFaceNew = DataStruct.Field2D {
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}
ElcCell = DataStruct.Field2D {
   -- (electric field stored on cell centers)
   onGrid = grid,
   -- [Ex, Ey, Ez]
   numComponents = 3,
   ghost = {1, 1},
}
ElcCellNew = DataStruct.Field2D {
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
ElcFaceDup = ElcFaceNew:duplicate()
MgnDup = MgnNew:duplicate()

function initEz(x, y, z)
   local m, n = 8, 5
   local a = m*Lucee.Pi/(upperx-lowerx)
   local b = n*Lucee.Pi/(uppery-lowery)
   local Ez = math.sin(a*x)*math.sin(b*y)
   return 0.0, 0.0, Ez
end
-- initialize electric field
ElcFace:set(initEz)
ElcFaceNew:copy(ElcFace)

-- write initial conditions
Mgn:write("magneticField_0.h5")
ElcFace:write("electricField_0.h5")

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
mgnUpdate:setIn( {Mgn, ElcFace} )
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
elcUpdate:setIn( {ElcFace, MgnNew} )
elcUpdate:setOut( {ElcFaceNew} )

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
bcElcFaceLeft = createBc(0, "lower", ElcFace, bcCopyFlip)
bcElcFaceRight = createBc(0, "upper", ElcFace, bcCopyFlip)
bcElcFaceTop = createBc(1, "upper", ElcFace, bcCopyFlip)
bcElcFaceBottom = createBc(1, "lower", ElcFace, bcCopyFlip)

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
   runUpdater(bcElcFaceLeft, tCurr, myDt)
   runUpdater(bcElcFaceRight, tCurr, myDt)
   runUpdater(bcElcFaceBottom, tCurr, myDt)
   runUpdater(bcElcFaceTop, tCurr, myDt)

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
   while true do
      ElcFaceDup:copy(ElcFace)
      MgnDup:copy(Mgn)
      if (tCurr+myDt > tEnd) then -- hit tEnd exactly
	 myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = singleStep(tCurr, myDt)

      if ( status == false ) then
	 Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 ElcFace:copy(ElcFaceDup)
	 Mgn:copy(MgnDup)
      else
	 ElcFace:copy(ElcFaceNew)
	 Mgn:copy(MgnNew)

	 tCurr = tCurr + myDt
	 step = step + 1
	 if (tCurr >= tEnd) then
	    break
	 end
      end
   end

   return dtSuggested
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 150e-9

-- compute time-step to move magnetic field to dt/2
dtSuggested = math.min(0.45*dx/light_speed, 0.45*dy/light_speed)
dt2 = 0.5*dtSuggested
-- move the magnetic field by -dt/2
runUpdater(bcElcFaceLeft, tStart, -dt2)
runUpdater(bcElcFaceRight, tStart, -dt2)
runUpdater(bcElcFaceBottom, tStart, -dt2)
runUpdater(bcElcFaceTop, tStart, -dt2)
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
   ElcFace:write( string.format("electricField_%d.h5", frame) )
   Mgn:write( string.format("magneticField_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
