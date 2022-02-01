from pylab import *
import postgkyl

d = postgkyl.GData("c-s3-2d-adv-dg_distf_0.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, v0 = dg.interpolate(0)

d = postgkyl.GData("c-s3-2d-adv-dg_distf_1.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, v1 = dg.interpolate(0)

dv = v1-v0
err = dv**2

dx = X[0][1]-X[0][0]
dy = X[1][1]-X[1][0]

totalErr = sqrt(dx*dy*sum(err))
print("dx = %g; err = %g" % (dx, totalErr))

