-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local Constants = require "Lib.Constants"
local Mpi = require "Comm.Mpi"

local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
local seed = rank^2 + 3
math.randomseed(seed)

----------------------
-- HELPER FUNCTIONS --
----------------------
local function log(...)
   if rank == 0 then
      print(string.format(...))
   end
end

----------------
-- PARAMETERS --
----------------
local lightSpeed = Constants.SPEED_OF_LIGHT/100
local mu0 = Constants.MU0
local epsilon0 = 1/mu0/(lightSpeed^2)
local gasGamma = 5/3

local rho0 = 6e6*Constants.PROTON_MASS
local vx0 = 450e3
local vy0 = 0e3
local vz0 = 0e3
local p0 = 6e-12
local Bx0 = 0e-9
local By0 = 0e-9
local Bz0 = -5e-9
local Bnoise = Bz0*1e-3

local Re = 6400e3
local mass_ratio = 1.0 -- 25.0
local pressure_ratio = 1 -- 5.0
local d0_i = 0.1 * Re
local ionMass = Constants.PROTON_MASS
local eleMass = ionMass/mass_ratio
local ionCharge = ionMass / d0_i / math.sqrt(mu0*rho0)
local eleCharge = -ionCharge

local rho0_e = rho0 / (1 + mass_ratio)
local rho0_i = rho0 - rho0_e
local p0_e = p0 / (1 + pressure_ratio)
local p0_i = p0 - p0_e

local xlo, xup, Nx =  -18 * Re, 78 * Re, 512
local ylo, yup, Ny =  -48 * Re, 48 * Re, 512
local lower = {xlo, ylo}
local upper = {xup, yup}
local cells = {Nx, Ny}

local tEnd = 1500
local nFrame = 100
local cfl = 0.5

-- diagnostic parameters
local u0_e = p0_e / (gasGamma-1) + 0.5 * (vx0^2 + vy0^2 + vz0^2) * rho0_e
local u0_i = p0_i / (gasGamma-1) + 0.5 * (vx0^2 + vy0^2 + vz0^2) * rho0_i

local B0 = Bz0
local wc0_e = math.abs(eleCharge * B0 / eleMass)
local vA0_e = B0 / math.sqrt(mu0 * rho0_e)
local vt0_e = math.sqrt(3 * p0_e / rho0_e)
local rt0_e = vt0_e / wc0_e

local wp0_i = lightSpeed / d0_i
local d0_e = d0_i / math.sqrt(mass_ratio)
local wp0_e = lightSpeed / d0_e
--
local Lx = xup - xlo
local Ly = yup - ylo
local dx = Lx / Nx

log("%30s = %g", "lightSpeed", lightSpeed)
log("%30s = %g", "mu0", mu0)
log("%30s = %g", "epsilon0", epsilon0)
log("%30s = %g", "gasGamma", gasGamma)

log("%30s = %g", "eleMass", eleMass)
log("%30s = %g", "eleCharge", eleCharge)
log("%30s = %g", "ionMass", ionMass)
log("%30s = %g", "ionCharge", ionCharge)

log("%30s = %g", "rho0_e", rho0_e)
log("%30s = %g", "p0_e", p0_e)
log("%30s = %g", "u0_e", u0_e)
log("%30s = %g", "rho0_i", rho0_i)
log("%30s = %g", "p0_i", p0_i)
log("%30s = %g", "u0_i", u0_i)
log("%30s = %g", "vx0", vx0)
log("%30s = %g", "vy0", vy0)
log("%30s = %g", "vz0", vz0)

log("%30s = %g", "Bx0", Bx0)
log("%30s = %g", "By0", By0)
log("%30s = %g", "Bz0", Bz0)

log("%30s = %g = 1/%g", "wp0_e", wp0_e, 1 / wp0_e)
log("%30s = %g = 1/%g", "wc0_e", wc0_e, 1 / wc0_e)
log("%30s = %g", "vt0_e=sqrt(2*p0_e/rho0_e)", vt0_e)
log("%30s = %g", "vA0_e", vA0_e)
log("%30s = %g", "d0_e", d0_e)
log("%30s = %g", "rt0_e=vt0_e/wc0_e", rt0_e)


log("%30s = %g", "dx", dx)
log("%30s = %g", "dx/d0_e", dx / d0_e)
log("%30s = %g", "dx/rt0_e", dx / rt0_e)

log("%30s = %g", "Lx", Lx)
log("%30s = %g", "Ly", Ly)
log("%30s = %g", "Nx", Nx)
log("%30s = %g", "Ny", Ny)
log("%30s = %g", "tEnd", tEnd)
log("%30s = %g", "nFrame", nFrame)

-----------------------
-- INITIAL CONDITION --
-----------------------
function calcV(x, y)
   local vx = vx0
   local vy = vy0
   local vz = vz0
   return vx, vy, vz
end

function calcB(x,y,z)
   local Bx = Bx0
   local By = By0
   local Bz = Bz0 + Bnoise*math.sin(2*math.pi*x/Lx)*math.sin(2*math.pi*y/Ly)
   return Bx, By, Bz
end

---------
-- APP --
---------
momentApp = Moments.App {
   logToFile = true,

   lower = lower,
   upper = upper,
   cells = cells,
   tEnd = tEnd,
   nFrame = nFrame,
   timeStepper = "fvDimSplit",
   cfl = cfl,

   -- decomposition for configuration space
   decompCuts = {8, 8}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   ele = Moments.Species {
      charge = eleCharge, mass = eleMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      --forceInv = true,
      limiter = "zero",
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rho_e = rho0_e
         local vx, vy, vz = calcV(x, y)
         local p_e = p0_e

         local rhovx_e = rho_e * vx
         local rhovy_e = rho_e * vy
         local rhovz_e = rho_e * vz
         local u_e = p_e / (gasGamma - 1)
            + 0.5 * (rhovx_e^2 + rhovy_e^2 + rhovz_e^2) / rho_e
         return rho_e, rhovx_e, rhovy_e, rhovz_e, u_e
      end,
      evolve = true,
   },

   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      --forceInv = true,
      limiter = "zero",
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rho_i = rho0_i
         local vx, vy, vz = calcV(x, y)
         local p_i = p0_i

         local rhovx_i = rho_i * vx
         local rhovy_i = rho_i * vy
         local rhovz_i = rho_i * vz
         local u_i = p_i / (gasGamma - 1)
            + 0.5 * (rhovx_i^2 + rhovy_i^2 + rhovz_i^2) / rho_i
         return rho_i, rhovx_i, rhovy_i, rhovz_i, u_i
      end,
      evolve = true,
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      limiter = "zero",
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local vx, vy, vz = calcV(x, y, 0)
         local Bx, By, Bz = calcB(x, y, Bnoise)

         local Ex = -vy*Bz + vz*By
         local Ey = -vz*Bx + vx*Bz
         local Ez = -vx*By + vy*Bx
         return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true,
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"ele", "ion"},
      timeStepper = "direct",
   },
}

momentApp:run()
