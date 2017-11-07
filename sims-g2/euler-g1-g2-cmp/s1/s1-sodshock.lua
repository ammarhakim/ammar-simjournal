-- Program to solve Euler equations

-- pressure
pressure = 1.0 -- [Pa]
-- density
density = 1.0 -- [kg/m^3]
-- adiabatic index
gas_gamma = 1.4
-- nominal speed of sound
sound_speed = math.sqrt(gas_gamma*pressure/density) -- [m/s]
-- CFL number
cfl = 0.9

-- global parameters
lowerx = 0.0
upperx = 1.0
nx = 1024

-- decomposition object to use
decomp = DecompRegionCalc1D.CartProd { cuts = {1} }
-- computational domain
grid = Grid.RectCart1D {
   lower = {lowerx},
   upper = {upperx},
   cells = {nx},
   decomposition = decomp,
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   -- [rho, rho*u, rho*v, rho*w, Er]
   numComponents = 5,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   -- [rho, rho*u, rho*v, rho*w, Er]
   numComponents = 5,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- initial condition to apply
function init(x,y,z)
   local sloc = 0.5 -- location of shock
   local rho, pr = density, pressure -- density, pressure on right
   if (x<sloc) then
      rho = 3*density
      pr = 3*pressure
   end
   local er = pr/(gas_gamma-1)
   return rho, 0.0, 0.0, 0.0, er
end
-- apply initial conditions
q:set(init)
-- sync ghost cell values
q:sync()
-- copy initial condition
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- define equation to solve
eulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gas_gamma,
}

-- updater for Euler equations
eulerSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = eulerEqn,
   -- one of no-limiter, minmod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "monotonized-centered",
   cfl = cfl,
   cflm = 1.01*cfl,
}

-- set input/output arrays (these do not change so set it once)
eulerSlvr:setIn( {q} )
eulerSlvr:setOut( {qNew} )

myDt = 100.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 0.1

tCurr = tStart
step = 1

-- main loop
while true do
   -- copy qNew in case we need to take this step again
   qNewDup:copy(qNew)

   -- if needed adjust dt to hit tEnd exactly
   if (tCurr+myDt > tEnd) then
      myDt = tEnd-tCurr
   end

   Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

   -- set current time
   eulerSlvr:setCurrTime(tCurr)
   -- advance solution
   status, dtSuggested = eulerSlvr:advance(tCurr+myDt)

   if (dtSuggested < myDt) then
      -- time-step too large
      Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      myDt = dtSuggested
      qNew:copy(qNewDup)
   else
      -- apply copy BCs on lower and upper edges
      qNew:applyCopyBc(0, "lower")
      qNew:applyCopyBc(0, "upper")

      -- sync values in ghost cells
      qNew:sync()

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
-- write final solution
q:write("q_1.h5")

Lucee.logInfo(string.format("WavePropagation1D updater took: %g sec", eulerSlvr:totalAdvanceTime()))
