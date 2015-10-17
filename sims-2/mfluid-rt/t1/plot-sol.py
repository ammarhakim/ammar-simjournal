import matplotlib
matplotlib.use('Agg')
from pylab import *
import tables

me = 1.0/25.0
qe = -1.0

def Jze(qbym, q):
    return qbym*q[:,:,:,3]

def numElc(m, q):
    return q[:,:,:,0]/m

def getMeshGrid(grid):
    xl, yl, zl = grid._v_attrs.vsLowerBounds
    xu, yu, zu = grid._v_attrs.vsUpperBounds
    nx, ny, nz = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny    
    dz = (zu-zl)/nz
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)
    Z = linspace(zl+0.5*dz, zu-0.5*dz, nz)

    return X, Y, Z

for i in range(0, 28):
    print "Working on %d ..." % i
    fh = tables.openFile("t1-5m-3d-rt_q_%d.h5" % i)
    q = fh.root.StructGridField
    X, Y, Z = getMeshGrid(fh.root.StructGrid)
    tm = fh.root.timeData._v_attrs.vsTime

    f = figure(1)
    ax = subplot(1,2,1)
    im = pcolormesh(X, Y, q[:,:,Z.shape[0]/2,0].transpose())
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    axis('image')

    ax = subplot(1,2,2)
    print(q[:,Y.shape[0]/2,:,0].shape)
    im = ax.pcolormesh(X, Z, q[:,Y.shape[0]/2,:,0].transpose())
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    axis('image')

    suptitle('T=%g' % (tm/10.0))
    savefig('t1-5m-3d-rt_numElc_%05d.png' % i, bbox_inches='tight')#, dpi=300)
    close()
    fh.close()

