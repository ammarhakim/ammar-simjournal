import pylab
import tables
import math
import numpy

def exactSol(a, b, X, Y):
    c1, d0 = 0, 0
    c0 = a/12.0 - 1.0/2.0
    d1 = b/12.0 - 1.0/2.0

    XX, YY = pylab.meshgrid(X, Y)
    return (XX**2/2 - a*XX**4/12 + c0*XX + c1)*(YY**2/2 - b*YY**4/12 + d0*YY + d1)

fh = tables.openFile("s99-poisson-o3-2d_phi.h5")
q = fh.root.StructGridField
nx, ny, nc = q.shape

# node 1
Xf = pylab.linspace(0, 1, nx)
Yf = pylab.linspace(0, 1, ny)
q_1 = q[:,:,0]

dx = Xf[1]-Xf[0]
dy = Yf[1]-Yf[0]
hsize = math.sqrt(dx*dy)

# node 5
Xf_5 = pylab.linspace(0.5*dx, 1-0.5*dx, nx-1)
Yf_5 = pylab.linspace(0, 1, ny)
q_5 = q[:-1,:,1]

# node 8
Xf_8 = pylab.linspace(0, 1, nx)
Yf_8 = pylab.linspace(0.5*dy, 1-0.5*dy, ny-1)
q_8 = q[:,:-1,2]

Xhr = pylab.linspace(0, 1, 101)
Yhr = pylab.linspace(0, 1, 101)
fhr = exactSol(2, 5, Xhr, Yhr)
vmin, vmax = fhr.min(), fhr.max()

# make plots
f1 = pylab.figure(1)

pylab.subplot(1,2,1)
cax = pylab.pcolormesh(Xf, Yf, pylab.transpose(q[:,:,0]), vmin=vmin, vmax=vmax)
pylab.axis('image')

pylab.subplot(1,2,2)
cax = pylab.pcolormesh(Xhr, Yhr, fhr)
pylab.axis('image')

# compute error
fex = exactSol(2, 5, Xf, Yf)
error = numpy.abs(fex-pylab.transpose(q_1)).sum()

fex = exactSol(2, 5, Xf_5, Yf_5)
error = error + numpy.abs(fex-pylab.transpose(q_5)).sum()

fex = exactSol(2, 5, Xf_8, Yf_8)
error = error + numpy.abs(fex-pylab.transpose(q_8)).sum()

error = error/(nx*ny + (nx-1)*ny + ny*(nx-1))

print hsize, error

pylab.show()
