-- Test code for moment calculators

log = Lucee.logInfo

polyOrder = 1
VDIM = 3
nMom = VDIM
nPrs = VDIM*(VDIM+1)/2

----------------------------
-- Grids, basis functions --
----------------------------
phaseGrid = Grid.RectCart5D {
   lower = {-1.0, -1.0, -6.0, -6.0, -6.0},
   upper = {1.0, 1.0, 6.0, 6.0, 6.0},
   cells = {32, 4, 8, 8, 8},
   periodicDirs = {0},
}
phaseBasis = NodalFiniteElement5D.SerendipityElement {
   onGrid = phaseGrid,
   polyOrder = polyOrder,
}
confGrid = Grid.RectCart2D {
   lower = {phaseGrid:lower(0), phaseGrid:lower(1)},
   upper = {phaseGrid:upper(0), phaseGrid:upper(1)},
   cells = {phaseGrid:shape(0), phaseGrid:shape(1)},
   periodicDirs = {0},
}
confBasis = NodalFiniteElement2D.SerendipityElement {
   onGrid = confGrid,
   polyOrder = polyOrder,
}

---------------------
-- Data structures --
---------------------

distf = DataStruct.Field5D {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numNodes(),
   ghost = {1, 1},
}
numDensity = DataStruct.Field2D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
momDensity = DataStruct.Field2D {
   onGrid = confGrid,
   numComponents = nMom*confBasis:numNodes(),
   ghost = {1, 1},
}
pressureTensor = DataStruct.Field2D {
   onGrid = confGrid,
   numComponents = nPrs*confBasis:numNodes(),
   ghost = {1, 1},
}
ptclEnergy = DataStruct.Field2D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

--------------
-- Updaters --
--------------

numDensityCalc = Updater.DistFuncMomentCalc2X3V {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = 0,
}
momDensityCalc = Updater.DistFuncMomentCalc2X3V {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = 1,
}		
pressureTensorCalc = Updater.DistFuncMomentCalc2X3V {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = 2,
}
ptclEnergyCalc = Updater.DistFuncMomentCalc2X3V {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = 2,
   scalarPtclEnergy = true,
}

-- initial condition to apply
function maxwellian(x,y,vx,vy,vz)
   local Pi = math.pi   
   local n = 1.0*math.sin(2*Pi*x)
   local ux = 0.1*math.cos(2*Pi*x)
   local uy = 0.2*math.sin(2*Pi*x)
   local uz = 0.1*math.cos(2*Pi*x)

   local Txx = 0.75 + 0.25*math.cos(2*Pi*x)
   local Tyy = 0.75 + 0.25*math.sin(2*Pi*x)
   local Tzz = 0.75 + 0.1*math.sin(2*Pi*x)
   local Txy = 0.5 + 0.1*math.sin(2*Pi*x)
   local Txz = 0.25 + 0.1*math.sin(2*Pi*x)
   local Tyz = 0.125 + 0.1*math.sin(2*Pi*x)

   local cx = vx-ux
   local cy = vy-uy
   local cz = vz-uz

   local detT = Txx*(Tyy*Tzz-Tyz^2)-Txy*(Txy*Tzz-Txz*Tyz)+Txz*(Txy*Tyz-Txz*Tyy)
   local u2 = cx*(cx*(Tyy*Tzz-Tyz^2)+cy*(Txz*Tyz-Txy*Tzz)+cz*(Txy*Tyz-Txz*Tyy))+cy*(cx*(Txz*Tyz-Txy*Tzz)+cy*(Txx*Tzz-Txz^2)+cz*(Txy*Txz-Txx*Tyz))+cz*(cx*(Txy*Tyz-Txz*Tyy)+cy*(Txy*Txz-Txx*Tyz)+cz*(Txx*Tyy-Txy^2))
   u2 = u2/(2*detT)

   return n/(math.pow(2*Pi, VDIM/2)*math.sqrt(detT))*math.exp(-u2)
end


initDistf = Updater.ProjectOnNodalBasis5D {
   onGrid = phaseGrid,
   basis = phaseBasis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   evaluate = function (x,y,vx,vy,vz,t)
      return maxwellian(x,y,vx,vy,vz)
   end
}
initDistf:setOut( {distf} )
initDistf:setCurrTime(0.0)
initDistf:advance(0.0) -- does not matter for initial conditions
distf:write("distf.h5")

---------------
-- Functions --
---------------

function calcNumDensity(curr, dt, distfIn, numDensOut)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:setOut( {numDensOut} )
   numDensityCalc:advance(curr+dt)
end
calcNumDensity(0.0, 0.0, distf, numDensity)
numDensity:write("numDensity.h5")

function calcMomDensity(curr, dt, distfIn, momDensOut)
   momDensityCalc:setCurrTime(curr)
   momDensityCalc:setIn( {distfIn} )
   momDensityCalc:setOut( {momDensOut} )
   momDensityCalc:advance(curr+dt)
end
calcMomDensity(0.0, 0.0, distf, momDensity)
momDensity:write("momDensity.h5")

function calcPressureTensor(curr, dt, distfIn, pressureTensorOut)
   pressureTensorCalc:setCurrTime(curr)
   pressureTensorCalc:setIn( {distfIn} )
   pressureTensorCalc:setOut( {pressureTensorOut} )
   pressureTensorCalc:advance(curr+dt)
end
calcPressureTensor(0.0, 0.0, distf, pressureTensor)
pressureTensor:write("pressureTensor.h5")

function calcPtclEnergy(curr, dt, distfIn, ptclEnergyOut)
   ptclEnergyCalc:setCurrTime(curr)
   ptclEnergyCalc:setIn( {distfIn} )
   ptclEnergyCalc:setOut( {ptclEnergyOut} )
   ptclEnergyCalc:advance(curr+dt)
end
calcPtclEnergy(0.0, 0.0, distf, ptclEnergy)
ptclEnergy:write("ptclEnergy.h5")

log(string.format("Initialization took %g [s]", initDistf:totalAdvanceTime()))
log(string.format("Moment calculation took %g [s]",
		  numDensityCalc:totalAdvanceTime()+
		     momDensityCalc:totalAdvanceTime()+
		     pressureTensorCalc:totalAdvanceTime()+
		     ptclEnergyCalc:totalAdvanceTime()))
