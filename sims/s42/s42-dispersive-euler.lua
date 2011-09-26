-- Program to dispersive Euler equations

-- simulation parameters
cfl = 0.9

-- global parameters
charge = 10.0
mass = 1.0
gasGamma = 2.0
Bz = 1.0
lambda = charge/mass

-- computational domain
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {200},
}

-- solution (We store 11 components as this allows use of the
-- LorentzForce to compute source term. The field quanties are not
-- evolved.)
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 11,
   ghost = {2, 2},
}

-- updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 11,
   ghost = {2, 2},
}

-- create duplicate copy in case we need to take step again
qNewDup = qNew:duplicate()

q:clear(0.0) -- zap everything out

-- create aliases to various sub-system
fluid = q:alias(0, 5)
mgnFld = q:alias(8, 11)
fluidNew = qNew:alias(0, 5)

-- initial conditions to apply
function initFluid(x,y,z)
   local rho0, p0 = 1.0, 1.0
   local U0 = 1e-8 -- velocity perturbation
   local nModes = 9
   local cs0 = math.sqrt(gasGamma*p0/rho0)
   local wc = lambda*Bz

   -- sum over modes to use exact solution to initialize problem
   local u, rho, v, p = 0.0, 0.0, 0.0, 0.0
   for n = 0, nModes do
      local kn = 2*Lucee.Pi*(2*n+1)
      local wn = math.sqrt(kn^2*cs0^2 + wc^2)
      local u1 = - U0/(2*n+1)*math.sin(kn*x)
      u = u + u1
      rho = rho - kn*rho0/wn*u1
      v = v - lambda*Bz/wn*U0/(2*n+1)*math.cos(kn*x)
      p = p - gasGamma*kn*p0/wn*u1
   end
   rho = rho0 + rho
   p = p0 + p
   return rho, rho*u, rho*v, 0.0, p/(gasGamma-1) + 0.5*rho*(u^2 + v^2)
end
fluid:set(initFluid)

-- "magnetic" field (this is set once and never evolved)
function initField(x,y,z)
   return 0.0, 0.0, Bz
end
mgnFld:set(initField)

-- copy initial conditions over
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- define various equations to solve
eulerEqn = HyperEquation.Euler {
   -- gas adiabatic constant
   gasGamma = gasGamma,
}

-- updater for electron equations
eulerSlvr = Updater.WavePropagation1D {
   onGrid = grid,
   equation = eulerEqn,
   -- one of no-limiter, min-mod, superbee, van-leer, monotonized-centered, beam-warming
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
}
-- initialize updater
eulerSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
eulerSlvr:setIn( {fluid} )
eulerSlvr:setOut( {fluidNew} )

-- Lorentz force on fluid (this works as the source uses the "magnetic field"
-- quanties stored in the qNew field to compute the source terms)
force = PointSource.LorentzForce {
   -- takes density, momentum and EM fields
   inpComponents = {0, 1, 2, 3, 5, 6, 7, 8, 9, 10},
   -- sets momentum and energy source
   outComponents = {1, 2, 3, 4},

   -- species charge and mass
   charge = charge,
   mass = mass,
}

-- updater to solve ODEs for source-term splitting scheme
sourceSlvr = Updater.GridOdePointIntegrator1D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {force},
}
-- initialize updater
sourceSlvr:initialize()
-- set input/output arrays (these do not change so set it once)
sourceSlvr:setOut( {qNew} )

-- function to take one time-step
function solveSystem(tCurr, t)
   -- advance fluids
   eulerSlvr:setCurrTime(tCurr)
   status, dtSuggested = eulerSlvr:advance(t)

   -- update source terms
   if (status) then
      sourceSlvr:setCurrTime(tCurr)
      sourceSlvr:advance(t)
   end

   return status, dtSuggested
end

myDt = 100.0 -- initial time-step to use (this will be discarded and adjusted to CFL value)
-- parameters to control time-stepping
tStart = 0.0
tEnd = 3.0

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
   -- advance fluids
   status, dtSuggested = solveSystem(tCurr, tCurr+myDt)

   if (dtSuggested < myDt) then
      -- time-step too large
      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      myDt = dtSuggested
      qNew:copy(qNewDup)
   else
      -- apply periodic BCs
      qNew:applyPeriodicBc(0)

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