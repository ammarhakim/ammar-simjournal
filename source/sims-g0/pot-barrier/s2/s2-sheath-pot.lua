local Plasma = require "Vlasov"

-- jump in phi
deltaPhi = -0.4 -- jump in potential

-- ion parameters
ion_udrift = 0.6 -- inlet drift-speed
ion_vth = 0.5 -- thermal speed

-- electron parameters
elc_udrift = 0.0 -- inlet drift-speed
elc_vth = 0.75 -- thermal speed

function maxwellian(v, n, u, vth)
   return n/math.sqrt(2*math.pi*vth^2)*math.exp(-(v-u)^2/(2*vth^2))
end

local app = Plasma.App { 
   tEnd = 40.0,
   nFrame = 40,
   lower = { 0.0 },
   upper = { 1.0 },
   cells = { 20 },
   basis = "serendipity",
   polyOrder = 2,

   useGPU = false,

   ion = Plasma.Species {
      charge = 1.0,  mass = 1.0,
      lower = {-4.0},
      upper = { 4.0},
      cells = {48},
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian(v, 1.0, ion_udrift, ion_vth)
      end,
      diagnostics = { "M0", "M1i", "M2" },

      bcx = { "bc_copy", "bc_absorb" }
   },

   elc = Plasma.Species {
      charge = -1.0,  mass = 1.0,
      lower = {-4.0},
      upper = { 4.0},
      cells = {48},
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 return maxwellian(v, 1.0, elc_udrift, elc_vth)
      end,
      diagnostics = { "M0", "M1i", "M2" },

      bcx = { "bc_copy", "bc_absorb" }
   },   

   field = Plasma.Field {
      epsilon0 = 1.0,  mu0 = 1.0,
      init = function (t, xn)
	 local x = xn[1]
	 local Ex = 0.0
	 local L = 5/20.0
	 
	 if x>0.5 and x<(0.5+L) then
	    Ex = -deltaPhi/L
	 end
	 
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,

      evolve = false,
      
   },
}
app:run()
