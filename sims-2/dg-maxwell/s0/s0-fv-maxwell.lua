-- Program to solve Maxwell equations

-- [[Add code to compute field energy]]

-- CFL number
cfl = 0.9

L = 1.0
X = L -- [m]
Y = L -- [m]
kwave = 2
lwave = 0
freq = 2*Lucee.Pi/L*math.sqrt(kwave^2+lwave^2)*Lucee.SpeedOfLight
tperiod = 2*Lucee.Pi/freq

print(freq, tperiod)

-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {X, Y},
   cells = {100, 1},
   decomposition = decomp,
   periodicDirs = {0, 1},
}

-- solution
q = DataStruct.Field2D {
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
   local cos = math.cos
   local pi = Lucee.Pi
   local c = Lucee.SpeedOfLight
   local phi = 2*pi/L*(kwave*x + lwave*y)
   local E0 = 1.0
   local Ex, Ey = 0.0, 0.0
   local Ez = E0*cos(phi)
   local Bx = E0/c*cos(phi)*2*pi/L*lwave
   local By = -E0/c*cos(phi)*2*pi/L*kwave
   local Bz = .0
   return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
end

-- apply initial conditions
q:set(init)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- define equation to solve
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = Lucee.SpeedOfLight,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 1.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 1.0,
}

-- updater for Euler equations
maxSlvr = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "no-limiter", 
   cfl = cfl,
   cflm = cfl*1.01
}

-- set input/output arrays (these do not change so set it once)
maxSlvr:setIn( {q} )
maxSlvr:setOut( {qNew} )

-- total energy calculator
emEnergy = DataStruct.DynVector { numComponents = 1 }
emEnergyCalc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  local epsilon0 = Lucee.Epsilon0
		  local mu0 = Lucee.Mu0
		  return 0.5*epsilon0*(ex^2+ey^2+ez^2) + 0.5/mu0*(bx^2+by^2+bz^2)
	       end,
}
emEnergyCalc:setIn( {q} )
emEnergyCalc:setOut( {emEnergy} )

-- boundary condition to apply
function applyBc(fld)
   fld:sync()
end

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
	 applyBc(qNew)

	 -- copy updated solution back
	 q:copy(qNew)
	 emEnergyCalc:setCurrTime(tCurr)
	 emEnergyCalc:advance(tCurr+myDt)

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
tEnd = tperiod

nFrames = 1
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   -- write out data
   q:write( string.format("q_%d.h5", frame), tCurr+tFrame )
   emEnergy:write( string.format("emEnergy_%d.h5", frame) )
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

