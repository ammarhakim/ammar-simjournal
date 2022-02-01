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
    tmAlf = tm/(Lx/Valf)
    
    f = figure(1)
    pcolormesh(XX, YY, dat.transpose())
    axis('image')
    colorbar()
    title("T = %.4g" % tmAlf)
    
    savefig(nm)
    close()

for i in range(0,51):
    print ("Working on %d .." % i)
    fh = tables.openFile("../s434/s434-is-coal_q_%d.h5" % i)
    q = fh.root.StructGridField
    X, Y = getMeshGrid(fh.root.StructGrid)
    mkFig(fh, X, Y, q[:,:,3], 's434-Jze-%05d.png' % i)
    fh.close()
