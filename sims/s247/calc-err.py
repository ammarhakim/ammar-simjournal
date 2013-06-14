from pylab import *
import tables

def exactSol(X, t):
    return exp(-t)*sin(X)

fh = tables.openFile("s247-dg-diffuse_q_1.h5")
q = fh.root.StructGridField
nx, nc = q.shape
dx = 2*pi/nx
Xf = linspace(0, 2*pi-dx, nx)

Xhr = linspace(0, 2*pi, 101)
fhr = exactSol(Xhr, 1.0)

# make plot comparing exact to numerical solution
plot(Xhr, fhr, '-r', Xf, q[:,0], '-k')

# compute error
fex = exactSol(Xf, 1.0)
error = abs(fex-q[:,0]).sum()/nx;

print "%g %g" % (dx, error)

show()
