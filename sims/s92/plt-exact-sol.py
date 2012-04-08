import pylab
import tables
import math
import numpy

def exactSol(a, b, X):
    c0 = -(1/2.0 + a/12.0 + b/30.0)
    c1 = 0.0
    return X**2/2 + a*X**4/12 + b*X**6/30 + c0*X + c1

fh = tables.openFile("s92-poisson-o3-1d_phi.h5")
q = fh.root.StructGridField
nx, nc = q.shape
Xf = pylab.linspace(0, 1, nx)
qe = q[:,0]

dx = Xf[1]-Xf[0]
Xm = pylab.linspace(0.5*dx, 1-0.5*dx, nx-1)
qm = q[:-1,1]

a = 2.0
b = -12.0

Xhr = pylab.linspace(0, 1, 101)
fhr = exactSol(a, b, Xhr)

# make plot comparing exact to numerical solution
pylab.plot(Xhr, fhr, '-r', Xf, qe, 'ok', Xm, qm, 'ok')

# compute error
fex_e = exactSol(a, b, Xf)
fex_m = exactSol(a, b, Xm)
error = (numpy.abs(fex_e-qe).sum() + numpy.abs(fex_m-qm).sum())/(nx+nx-1);

print "%g %g" % (dx, error)

pylab.show()
