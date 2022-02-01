-- Input file to solve 1D Poisson equations using FEM

-- polynomial order
polyOrder = 1

-- Determine number of global nodes per cell for use in creating
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numNodesPerCell = 1
elseif (polyOrder == 2) then
   numNodesPerCell = 2
elseif (polyOrder == 3) then
   numNodesPerCell = 3
end

grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {32},
}

-- source term
src = DataStruct.Field1D {
   onGrid = grid,
   location = "vertex", -- this will not work in general for polyOrder > 1
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numNodesPerCell,
   ghost = {1, 1},
}

-- function to initialize source
function initSrc(x,y,z)
   return 1-2*x^2
end
-- initialize source
src:set(initSrc)

-- write it to disk
src:write("src.h5")

-- field to store potential
phi = DataStruct.Field1D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer of right to get nodes
}
-- clear out contents
phi:clear(0.0)

-- create FEM nodal basis
lobattoBasis = NodalFiniteElement1D.Lobatto {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, 2 or 3. Corresponding
   -- number of nodes are 2, 3, or 4.
   polyOrder = 1,
}

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson1D {
   onGrid = grid,
   -- basis functions to use
   basis = lobattoBasis,
   -- boundary conditions to apply
   bcLeft = { T = "D", V = 0.0 },
   bcRight = { T = "D", V = 0.0 },
}

-- set input output fields
poissonSlvr:setIn( {src} )
poissonSlvr:setOut( {phi} )

-- solve for potential
status, dtSuggested = poissonSlvr:advance(0.0) -- time is irrelevant here

-- output solution
phi:write("phi.h5")