-- Input file to solve transverse magnetic mode Maxwell equations

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

Elc = DataStruct.Field2D {
   -- (electric field is stored on edges)
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
Ex = ElcNew:alias(0, 1)
Ey = ElcNew:alias(1, 2)
Ez = ElcNew:alias(2, 3)

Mgn = DataStruct.Field2D {
   -- (magnetic field is stored on faces)
   onGrid = grid,
   -- [Bx, By, Bz]
   numComponents = 3,
   ghost = {1, 1},
}
Mgn:clear(0.0)

Bx = Mgn:alias(0, 1)
By = Mgn:alias(1, 2)
Bz = Mgn:alias(2, 3)

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
   local m, n = 8, 5
   local a = m*Lucee.Pi/(upperx-lowerx)
   local b = n*Lucee.Pi/(uppery-lowery)
   local xez = x-0.5*dx -- Ez is stored at corner
   local yez = y-0.5*dy
   local Ez = 1.0*math.sin(a*xez)*math.sin(b*yez)
   return 0.0, 0.0, Ez
end
-- initialize electric field
Elc:set(initEz)
ElcNew:copy(Elc)

-- write initial conditions
Elc:write("electricField_0.h5")

-- updater for electric field
elcUpdate = Updater.EdgeFaceCurl2D {
   onGrid = grid,
   -- ElcNew = Elc + c^2*dt*curl(Mgn)
   alpha = light_speed^2,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = 0.45,
}
-- initialize updater
elcUpdate:initialize()
-- set input/out arrays (these do not change so set it once)
elcUpdate:setIn( {Elc, Mgn} )
elcUpdate:setOut( {ElcNew} )

-- updater for magnetic field
mgnUpdate = Updater.FaceEdgeCurl2D {
   onGrid = grid,
   -- MgnNew = Mgn - dt*curl(ElcNew)
   alpha = -1.0,
   -- speed of propagation
   speed = light_speed,
   -- CFL numbers
   cfl = 0.45,
}
-- initialize updater
mgnUpdate:initialize()
-- set input/out arrays (these do not change so set it once)
mgnUpdate:setIn( {Mgn, ElcNew} )
mgnUpdate:setOut( {MgnNew} )

-- electric field BCs on right
bcElcRgt = BoundaryCondition.Const { components = {0, 1, 2}, values = {0, 0, 0} }
bcElcTop = BoundaryCondition.Const { components = {0, 1, 2}, values = {0, 0, 0} }

bcRight = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcElcRgt},
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "upper",
}
-- initialize updater
bcRight:initialize()
-- set input/output arrays (these do not change so set it once)
bcRight:setOut( {ElcNew} )

bcTop = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcElcTop},
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "upper",
}
-- initialize updater
bcTop:initialize()
-- set input/output arrays (these do not change so set it once)
bcTop:setOut( {ElcNew} )

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
      Mgn:applyCopyBc(1, "lower")

      -- update electric field
      elcUpdate:setCurrTime(tCurr)
      elcStatus, elcDtSuggested = elcUpdate:advance(tCurr+myDt)

      -- apply BCs on electric field
      bcRight:setCurrTime(tCurr)
      bcRight:advance(tCurr+myDt)
      bcTop:setCurrTime(tCurr)
      bcTop:advance(tCurr+myDt)

      -- update magnetic field
      mgnUpdate:setCurrTime(tCurr)
      mgnStatus, mgnDtSuggested = mgnUpdate:advance(tCurr+myDt)

      if ( (elcStatus == false) or (mgnStatus == false) ) then
	 dtSuggested = math.min(elcDtSuggested, mgnDtSuggested)
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

dtSuggested = 100.0 -- initial suggested time-step
-- parameters to control time-stepping
tStart = 0.0
tEnd = 150e-9

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
