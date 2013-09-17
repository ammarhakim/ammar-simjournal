import numpy
from pylab import *
import tables

rc('text', usetex=True)

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*elcMass*n0)
ionCycl = ionCharge*B0/ionMass
start = 0
end  = 60
nFrame = end-start+1
tm = zeros((nFrame,), float)
divb = zeros((nFrame,), float)

def calcDiv(dx, dy, fldx, fldy):
    nx, ny = fldx.shape[0], fldx.shape[1]

    divb = 0.0
    for iy in range(1,ny-1):
        ydiff = (fldy[:,iy+1]-fldy[:,iy-1])/(2*dy)
        divb = divb+sum(ydiff)

    for ix in range(1,nx-1):
        xdiff = (fldx[ix+1.:]-fldx[ix-1,:])/(2*dx)
        divb = divb+sum(xdiff)

    return divb/(nx*ny)

count = 0
for i in range(start, end+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("s285-pulsebox-wave_q_%d.h5" % i)
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]
    YI = ny/2

    X = linspace(-LX/2, LX/2, nx)
    Y = linspace(-LY/2, LY/2, ny)

    dx = X[1]-X[0]
    dy = Y[1]-Y[0]


    tm[count] = fh.root.timeData._v_attrs['vsTime']
    divb[count] = calcDiv(dx, dy, q[:,:,3], q[:,:,4])
    
    count = count+1

plot(tm, divb)
show()
