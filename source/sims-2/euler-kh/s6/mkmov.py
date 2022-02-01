from pylab import *
import tables

def getMeshGrid(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

    return meshgrid(X, Y)

def mkFig(fh, XX, YY, dat, nm):
    tm = fh.root.timeData._v_attrs.vsTime
    Valf = 0.1
    Lx = 4*pi*5.0
    tmAlf = tm
    
    f = figure(1)
    pcolormesh(XX, YY, dat.transpose())
    axis('image')
    colorbar()
    title("T = %.4g" % tmAlf)
    
    savefig(nm)
    close()

for i in range(101):
    print ("Working on %d .." % i)
    fh = tables.openFile("s6-euler-kh_q_%d.h5" % i)
    q = fh.root.StructGridField
    X, Y = getMeshGrid(fh.root.StructGrid)
    mkFig(fh, X, Y, q[:,:,0], 's6-euler-kh_rho_%05d.png' % i)
    fh.close()
