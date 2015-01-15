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
    Lx = 4*pi*10.0
    tmAlf = tm/(Lx/Valf)
    
    f = figure(1)
    pcolormesh(XX, YY, dat.transpose())
    axis('image')
    colorbar()
    title("T = %.4g" % tmAlf)
    
    savefig(nm)
    close()

from optparse import OptionParser
# set command line options
parser = OptionParser()
parser.add_option('-p', '--plot', action = 'store',
                  dest = 'fileName',
                  help = 'Hdf5 file to plot')
(options, args) = parser.parse_args()

fn = options.fileName

fh = tables.openFile(fn)
q = fh.root.StructGridField
X, Y = getMeshGrid(fh.root.StructGrid)
mkFig(fh, X, Y, q[:,:,3], "%s.png" % fn[:-3])



