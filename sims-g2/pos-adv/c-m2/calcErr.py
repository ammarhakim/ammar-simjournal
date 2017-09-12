from pylab import *
import postgkyl

d = postgkyl.GData("c-m2-2d-adv-dg_distf_0.bp")
q0 = d.q

d = postgkyl.GData("c-m2-2d-adv-dg_distf_1.bp")
q1 = d.q

dq = q1-q0
err = dq[0]**2 + dq[1]**2 + dq[2]**2 + dq[3]**2

dx = (d.upperBounds[0]-d.lowerBounds[0])/d.numCells[0]
dy = (d.upperBounds[1]-d.lowerBounds[1])/d.numCells[1]
vol = dx*dy

totalErr = dx*dy/4.0*sqrt(sum(err))
print("dx = %g; err = %g" % (dx, totalErr))

