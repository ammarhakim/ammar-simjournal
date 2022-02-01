import matplotlib
matplotlib.use('Agg')
from pylab import *
import tables

me = 1.0/25.0
qe = -1.0

def Jze(qbym, q):
    return qbym*q[:,:,3]

def numElc(m, q):
    return q[:,:,0]/m

def getMeshGrid(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

    return meshgrid(X, Y)

for i in range(0, 18):
    print "Working on %d ..." % i
    fh = tables.openFile("s4-5m-2d-rt_q_%d.h5" % i)
    XX, YY = getMeshGrid(fh.root.StructGrid)

    figure(1)
    q = fh.root.StructGridField
    ne = numElc(me, q)
    pcolormesh(YY, XX, transpose(ne)) #, cmap='hot')
    tm = 5.0*i
    title('T=%g' % tm)
    axis('image')

    savefig('s4-5m-2d-rt-numElc_%05d.png' % i, bbox_inches='tight')#, dpi=300)
    close()
    fh.close()

