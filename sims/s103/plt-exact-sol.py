from pylab import *
import tables
import math
import numpy

def fxy(nm, nn, amn, bmn, XX, YY):
    f = zeros(XX.shape, float)
    for m in range(nm):
        for n in range(nn):
            f = f + amn[m][n]*cos(m*XX)*cos(n*YY) + bmn[m][n]*sin(m*XX)*sin(n*YY)
    return f/50

def sxy(nm, nn, amn, bmn, XX, YY):
    f = zeros(XX.shape, float)
    for m in range(nm):
        for n in range(nn):
            f = f + -(m*m+n*n)*(amn[m][n]*cos(m*XX)*cos(n*YY) + bmn[m][n]*sin(m*XX)*sin(n*YY))
    return f/50

amn = [[0,10,0], [10,0,0], [10,0,0]]
bmn = [[0,10,0], [10,0,0], [10,0,0]]

fh = tables.openFile("s103-periodic-poisson-2d_phi.h5")
q = fh.root.StructGridField
nx, ny, nc = q.shape

xup = 2*math.pi

# node 1
Xf = linspace(0, xup, nx)
Yf = linspace(0, xup, ny)
q_1 = q[:,:,0]
q_1 = q_1 - q_1.max()

dx = Xf[1]-Xf[0]
dy = Yf[1]-Yf[0]
hsize = math.sqrt(dx*dy)

XX, YY = meshgrid(Xf, Yf)

Xhr = linspace(0, xup, 101)
Yhr = linspace(0, xup, 101)
XXhr, YYhr = meshgrid(Xhr, Yhr)

fhr = fxy(3, 3, amn, bmn, XXhr, YYhr)
fhr = fhr - fhr.max()
vmin, vmax = fhr.min(), fhr.max()

# make plots
f1 = figure(1)

subplot(1,2,1)
cax = pcolormesh(Xf, Yf, transpose(q_1))
colorbar()
axis('image')

subplot(1,2,2)
cax = pcolormesh(Xhr, Yhr, fhr)
colorbar()
axis('image')


# compute error
fex = fxy(3, 3, amn, bmn, XX, YY)
fex = fex - fex.max()
error = numpy.abs(fex-transpose(q_1)).sum()/(nx*ny)

print hsize, error

# compute source
fh = tables.openFile("s103-periodic-poisson-2d_src.h5")
q = fh.root.StructGridField
src = q[:,:,0]

srcEx = sxy(3, 3, amn, bmn, XX, YY)

show()
