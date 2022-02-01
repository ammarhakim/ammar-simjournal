-- Plasma ------------------------------------------------------------------------
local Basis = require "Basis"
local Constants = require "Lib.Constants"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Time = require "Lib.Time"
local Updater = require "Updater"

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0 = 40*eV 
Ti0 = 40*eV 
B_axis = 0.5 -- [T]
R0 = 0.85  -- [m]
a0  = 0.5 -- [m]
R = R0 + a0
B0 = B_axis*(R0/R) -- [T]
n0 = 7e18 -- [1/m^3]
P_SOL = 8.1e5 -- [W] 
S0 = 5.7691e23
xSource = R - 0.05 -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature
-- derived parameters
vti     = math.sqrt(Ti0/mi)
vte  	= math.sqrt(Te0/me)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s   = c_s/omega_ci

-- box size
Lx = 50*rho_s
Ly = 100*rho_s
Lz = 4 -- [m]

-- source profiles
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local sourceFloor = 0.1
   if math.abs(z) < Lz/4 then
      return 0.90625*S0*math.max(math.exp(-(x-xSource)^2/(2*lambdaSource)^2), sourceFloor)
   else
      return 1e-10
   end
end
sourceTemperature = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < xSource + 3*lambdaSource then
      return 80*eV
   else
      return 30*eV
   end
end

local confGrid = Grid.RectCart {
   lower = {R - Lx/2, -Ly/2, -Lz/2}, -- configuration space lower left
   upper = {R + Lx/2, Ly/2, Lz/2}, -- configuration space upper right
   cells = {24, 48, 16}, -- configuration space cells
}
local confBasis = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = 1 }

local srcDens = DataStruct.Field {
   onGrid = confGrid,
   numComponents = confBasis:numBasis(),
   ghost = {1, 1},
}
local srcTemp = DataStruct.Field {
   onGrid = confGrid,
   numComponents = confBasis:numBasis(),
   ghost = {1, 1},
}

local projectSrcDens = Updater.ProjectOnBasis {
   onGrid = confGrid,
   basis = confBasis,
   evaluate = sourceDensity,
}
projectSrcDens:advance(0.0, 0.0, {}, {srcDens})

local projectSrcTemp = Updater.ProjectOnBasis {
   onGrid = confGrid,
   basis = confBasis,
   evaluate = sourceTemperature,
}
projectSrcTemp:advance(0.0, 0.0, {}, {srcTemp})


srcDens:write("srcDens.bp")
srcTemp:write("srcTemp.bp")
