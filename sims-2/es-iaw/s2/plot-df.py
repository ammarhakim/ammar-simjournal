from pylab import *
import postgkyl

d = postgkyl.GData("s2-es-iaw_distfIon_0.h5")
dg0 = postgkyl.GInterpNodalSerendipity(d, 2)
X, fv0 = dg0.project(0)

d = postgkyl.GData("s2-es-iaw_distfIon_20.h5")
dg1 = postgkyl.GInterpNodalSerendipity(d, 2)
X, fv1 = dg1.project(0)

df = fv1-fv0

figure(1)
pcolormesh(X[0], X[1], df, cmap='inferno')
axis('tight')
colorbar()

figure(2)
plot(X[1], df[32,:])

figure(3)
plot(X[0], df[:,32])



show()


