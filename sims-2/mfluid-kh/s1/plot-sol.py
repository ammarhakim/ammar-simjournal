from pylab import *
import tables

me = 1.0/25.0
qe = -1.0

def Jze(qbym, q):
    return qbym*q[:,:,3]

def getMeshGrid(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

    return meshgrid(Y, X)

for i in range(0, 101):
    print "Working on %d ..." % i
    fh = tables.openFile("s1-5m-karim-kh_q_%d.h5" % i)
    YY, XX = getMeshGrid(fh.root.StructGrid)

    q = fh.root.StructGridField
    jze = Jze(qe/me, q)
    pcolormesh(YY, XX, abs(jze), cmap='hot')
    tm = 5.0*i
    title('T=%g' % tm)
    axis('image')
    savefig('s1-5m-karim-kh_%05d.png' % i, bbox_inches='tight')
    close()
    fh.close()

