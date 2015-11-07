import matplotlib
matplotlib.use('Agg')
from pylab import *
import tables
import vistools

me = 1.0/100.0
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
    fh = tables.openFile("s9-10m-karim-kh_q_%d.h5" % i)
    YY, XX = getMeshGrid(fh.root.StructGrid)

    figure(1)
    subplot(2,1,1)
    q = fh.root.StructGridField
    jze = Jze(qe/me, q)
    pcolormesh(YY, XX, abs(jze), cmap='hot')
    tm = 5.0*i
    title('T=%g' % tm)
    axis('image')

    subplot(2,1,2)
    Bx = q[:,:,23]
    By = q[:,:,24]
    psi = vistools.calc_psi2d(transpose(Bx), transpose(By))
    contour(YY, XX, transpose(psi), 40, colors='k', linestyles='solid')
    axis('image')

    savefig('s9-10m-karim-kh_%05d.png' % i, bbox_inches='tight', dpi=300)
    close()
    fh.close()

