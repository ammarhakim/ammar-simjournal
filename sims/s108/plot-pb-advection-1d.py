from pylab import *
import numpy
import tables
import math

def exactSolution(X):
    return exp(-75*(X-0.5)**2)

def projOnFiner(Xc, q):
    Xn = zeros( (2*Xc.shape[0],), float)
    qn = zeros( (2*Xc.shape[0],), float)
    dx = Xc[1]-Xc[0]
    
    Xn[0:2*nx:2] = Xc-0.25*dx
    Xn[1:2*nx:2] = Xc+0.25*dx

    qn[0:2*nx:2] = 3.0/4.0*q[:,2,0] + 1.0/4.0*q[:,2,1]
    qn[1:2*nx:2] = 1.0/4.0*q[:,2,0] + 3.0/4.0*q[:,2,1]

    return Xn, qn

# load data
fh = tables.openFile("s108-pb-advection-1d_chi_1.h5")
q = fh.root.StructGridField

nx = q.shape[0]
dx = 1.0/nx

Xc = linspace(0.5*dx, 1-0.5*dx, nx)
# get numerical solution
Xn, qn = projOnFiner(Xc, q)
# get initial condition on same mesh
fh = tables.openFile("s108-pb-advection-1d_chi_0.h5")
q = fh.root.StructGridField
Xn, qEx = projOnFiner(Xc, q)

# compute error
err = numpy.abs(qEx-qn).sum()/nx
print dx, err

# get exact solution on finer mesh
dx_hr = 1/200.
Xhr = linspace(0.5*dx_hr, 1-0.5*dx_hr, 200)
qEx_hr = exactSolution(Xhr)

# plot
plot(Xhr, qEx_hr, '-k', Xn, qn, '-r')

show()

