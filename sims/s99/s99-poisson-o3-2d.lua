-- Input file to solve 2D Poisson equations using FEM

-- polynomial order
polyOrder = 2

-- Determine number of global nodes per cell for use in creating
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numNodesPerCell = 1
elseif (polyOrder == 2) then
   numNodesPerCell = 3
end

grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {1.0, 1.0},
   cells = {16, 16},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.Serendipity {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 4 and 8.
   polyOrder = polyOrder,
}

-- source term
src = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex", -- this will not work in general for polyOrder > 1
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- function to initialize source
-- create updater to initialize source
initSrc = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 local a, b = 2, 5
		 local c1, d0 = 0, 0
		 local c0 = a/12 - 1/2
		 local d1 = b/12 - 1/2
		 local t1 = (1-a*x^2)*(-b*y^4/12 + y^2/2 + d0*y + d1)
		 local t2 = (1-b*y^2)*(-a*x^4/12 + x^2/2 + c0*x + c1)
		 return t1+t2
	      end
}
initSrc:setOut( {src} )
-- initialize potential
initSrc:advance(0.0) -- time is irrelevant
-- write it to disk
src:write("src.h5")

-- field to store potential
phi = DataStruct.Field2D {
   onGrid = grid,
   location = "vertex", -- this will not work in general for polyOrder > 1
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- clear out contents
phi:clear(0.0)

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- boundary conditions to apply
   bcLeft = { T = "D", V = 0.0 },
   bcRight = { T = "D", V = 0.0 },
   bcBottom = { T = "N", V = 0.0 },
   bcTop = { T = "D", V = 0.0 },
}

-- set input output fields
poissonSlvr:setIn( {src} )
poissonSlvr:setOut( {phi} )

-- solve for potential (time is irrelevant here)
status, dtSuggested = poissonSlvr:advance(0.0)
-- check if solver converged
if (status == false) then
   Lucee.logError("Poisson solver failed to converge!")
end

-- output solution
phi:write("phi.h5")