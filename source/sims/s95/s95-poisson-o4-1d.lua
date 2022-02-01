-- Input file to solve 1D Poisson equations using FEM

-- polynomial order
polyOrder = 3

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
   cells = {8},
}

-- create FEM nodal basis
basis = NodalFiniteElement1D.Lobatto {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, 2 or 3. Corresponding
   -- number of nodes are 2, 3, or 4.
   polyOrder = polyOrder,
}

-- source term
src = DataStruct.Field1D {
   onGrid = grid,
   location = "vertex", -- this will not work in general for polyOrder > 1
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numNodesPerCell,
   -- ghost cells
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- updater to initialize source
initSrc = Updater.EvalOnNodes1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local a, b = 2.0, -12.0
		 return 1 + a*x^2 + b*x^4
	      end
}
initSrc:setOut( {src} )
-- initialize source
initSrc:advance(0.0) -- time is irrelevant

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

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson1D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
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