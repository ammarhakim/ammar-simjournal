-- Program to solve Euler equations

-- gas adiabatic index
gasGamma = 1.4

-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {100},
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   -- [rho, rho*u, Er]
   numComponents = 3,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   -- [rho, rho*u, Er]
   numComponents = 3,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

-- initial condition to apply
function init(x,y,z)
   xs = 0.5
   rhol, ul, pl = 1.4, 0.1, 1.0
   rhor, ur, pr = 1.0, 0.1, 1.0
   if (x<xs) then
      rho, u, press = rhol, ul, pl
   else
      rho, u, press = rhor, ur, pr
   end
   return rho, rho*u, press/(gasGamma-1) + 0.5*rho*u*u
end

-- apply initial conditions
q:set(init)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- CFL number
mycfl = 0.45
-- updater for Euler equations
eulerSlvr = Updater.MusclHancock1D {
   onGrid = grid,
   gas_gamma = gasGamma,
   cfl = mycfl,
   -- one of: average, minmod, superbee, zero, epsilon
   limiter = "epsilon",
}

-- set input/output arrays (these do not change so set it once)
eulerSlvr:setIn( {q} )
eulerSlvr:setOut( {qNew} )

myDt = 100.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 2.0

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

   print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

   -- set current time
   eulerSlvr:setCurrTime(tCurr)
   -- advance solution
   status, dtSuggested = eulerSlvr:advance(tCurr+myDt)

   if (dtSuggested < myDt) then
      -- time-step too large
      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      myDt = dtSuggested
      qNew:copy(qNewDup)
   else
      -- apply copy BCs on lower and upper edges
      qNew:applyCopyBc(0, "lower")
      qNew:applyCopyBc(0, "upper")

      -- check if a nan occured
      if (qNew:hasNan()) then
	 print (string.format("** Nan occured at %g! Writing out data just before nan", tCurr))
	 q:write("q_pre_nan.h5")
	 break
      end

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