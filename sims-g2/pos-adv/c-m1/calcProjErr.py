from pylab import *
import postgkyl

d = postgkyl.GData("c-m1-2d-adv-dg_distf_0.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, v0 = dg.interpolate(0)

d = postgkyl.GData("c-m1-2d-adv-dg_distf_1.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, v1 = dg.interpolate(0)

dv = v1-v0
err = dv**2

dx = (d.upperBounds[0]-d.lowerBounds[0])/d.numCells[0]
dy = (d.upperBounds[1]-d.lowerBounds[1])/d.numCells[1]
vol = dx*dy

totalErr = sqrt(sum(err))/(d.numCells[0]*d.numCells[1])
print("dx = %g; err = %g" % (dx, totalErr))

