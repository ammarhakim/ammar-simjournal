local ffi = require "ffi"
local Time = require "Lib.Time"
local DataStruct = require "DataStruct"

-- simulation parameters
tEnd = 500*math.pi
w = 0.4567
dt = 0.001

-- Maxwellian
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

-- velocity list
vel = { -4, -2, -1, -0.2, -0.1, 0.1, 0.2, 1, 2, 4 }

local function traceParticles(outPrefix, E0)
   -- initialize particles
   ptcls = { }
   for i = 1, 10 do
      ptcls[i] = DataStruct.DynVector { numComponents = 4 }
      local vx = vel[i]
      local f = maxwellian2D(1.0, vx, 0.0, 0.0, 0.0, 1.0)
      -- x, vx, vy, f(x,vx,vy)
      ptcls[i]:appendData(0.0, { 0.0, vx, 0.0, f })
   end
   -- advance particles
   for i = 1, 10 do
      local tm, xv = ptcls[i]:lastData()
      local xcurr, vxcurr, vycurr, fcurr = xv[1], xv[2], xv[3], xv[4]
      
      local tCurr = 0.0
      while tCurr < tEnd do
	 local tHalf = tCurr+0.5*dt
	 local vxold = vxcurr

	 xcurr = xcurr + dt*vxcurr
	 vxcurr = ((1-dt^2/4)*vxold + dt*E0*math.cos(xcurr-w*tHalf) + dt*vycurr)/(1+dt^2/4)
	 vycurr = vycurr - dt/2*(vxcurr+vxold)

	 -- append data
	 ptcls[i]:appendData(tCurr, { xcurr, vxcurr, vycurr, fcurr })
	 tCurr = tCurr+dt
      end
   end

   -- write out solution
   for i = 1, 10 do
      ptcls[i]:write(string.format("%s_ptcl_%d.bp", outPrefix, i), 0.0)
   end
end

traceParticles("E0_2", 0.2)
traceParticles("E0_4", 0.4)
traceParticles("E0_6", 0.6)
traceParticles("E0_95", 0.95)

