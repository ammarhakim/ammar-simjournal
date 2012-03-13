-- Program to test Poisson solver in 2D on a periodic domain

-- basic parameters
L = 2.0
Nx = 128
Ny = 128

-- computational domain
grid = Grid.RectCart2D {
   lower = {-L/2, -L/2},
   upper = {L/2, L/2},
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
   return math.exp(-10*(x*x + y*y))
end

src:set(srcFunc)
-- write out source
src:write("src.h5")

-- create Poisson solver
poissonSlvr = Updater.PeriodicPoisson2D { onGrid = grid, }
-- set in/out arrays
poissonSlvr:setIn( {src} )
poissonSlvr:setOut( {sol} )

-- solve equation and write out solution
poissonSlvr:advance(0.0) -- time does not matter
sol:write("sol.h5")

-- now compute central-difference of solution
-- for central-difference of solution
solCD = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 1,
   ghost = {2, 2},
}

-- create CD updater
cdCalc = Updater.RectSecondOrderCentralDiff2D { onGrid = grid, }
-- set in/out arrays
cdCalc:setIn( {sol} )
cdCalc:setOut( {solCD} )

-- compute
cdCalc:advance(0.0) -- time does not matter
solCD:write("solCD.h5")
