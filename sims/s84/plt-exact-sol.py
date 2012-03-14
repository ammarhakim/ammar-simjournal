import pylab
import tables
import math
import numpy

def exactSol(a, X):
    c0 = 0.0
    c1 = a/12.0 - 1.0/2.0
    return X**2/2 - a*X**4/12 + c0*X + c1

fh = tables.openFile("s84-poisson-1d_phi.h5")
q = fh.root.StructGridField
nx, nc = q.shape
Xf = pylab.linspace(0, 1, nx)
dx = Xf[1]-Xf[0]

Xhr = pylab.linspace(0, 1, 101)
fhr = exactSol(5, Xhr)

# make plot comparing exact to numerical solution
pylab.plot(Xhr, fhr, '-r', Xf, q[:,0], '-ko')

# compute error
fex = exactSol(5, Xf)
error = numpy.abs(fex-q[:,0]).sum()/nx;

print "%g %g" % (dx, error)

pylab.show()
