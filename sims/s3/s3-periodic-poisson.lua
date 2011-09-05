-- Program to test Poisson solver in 2D on a periodic domain

-- basic parameters
L = 10.0
Nx = 128
Ny = 128

-- computational domain
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {L, L},
   cells = {Nx, Ny},
}

-- solution
sol = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}

-- source term
src = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}

-- function to set source
function srcFunc(x,y,z)
   local x1, y1 = 3.5, 5.0
   local x2, y2 = 6.5, 5.0
   local r1 = (x-x1)*(x-x1) + (y-y1)*(y-y1)
   local r2 = (x-x2)*(x-x2) + (y-y2)*(y-y2)
   return math.exp(-r1/0.8) + math.exp(-r2/0.8)
end

src:set(srcFunc)
-- write out source
src:write("src.h5")

-- create Poisson solver
poissonSlvr = Updater.PeriodicPoisson2D { onGrid = grid, }
-- initialize updater
poissonSlvr:initialize()
-- set in/out arrays
poissonSlvr:setIn( {src} )
poissonSlvr:setOut( {sol} )

-- solve equation and write out solution
poissonSlvr:advance(0.0) -- time does not matter
sol:write("sol.h5")
