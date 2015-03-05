import pylab
import numpy
from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 18
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
#rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

def colorbar_adj(obj, mode=1, redraw=False, _fig_=None, _ax_=None, aspect=None):
    '''
    Add a colorbar adjacent to obj, with a matching height
    For use of aspect, see http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.set_aspect ; E.g., to fill the rectangle, try "auto"
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if mode == 1:
        _fig_ = obj.figure; _ax_ = obj.axes
    elif mode == 2: # assume obj is in the current figure/axis instance
        _fig_ = plt.gcf(); _ax_ = plt.gca()
    _divider_ = make_axes_locatable(_ax_)
    _cax_ = _divider_.append_axes("right", size="5%", pad=0.05)
    _cbar_ =  _fig_.colorbar(obj, cax=_cax_)
    if aspect != None:
        _ax_.set_aspect(aspect)
    if redraw:
        _fig_.canvas.draw()
    return _cbar_


from pylab import *
import tables

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

def mkFig(fh, XX, YY, dat, nm):
    tm = fh.root.timeData._v_attrs.vsTime
    f = figure(1)
    im = pcolormesh(XX, YY, dat.transpose())
    title("T = %.4g" % tm)
    axis('image')
    colorbar_adj(im)
    
    savefig(nm, bbox_inches='tight')
    close()

for i in range(0,51):
    print ("Working on %d .." % i)
    fh = tables.openFile("s446-euler-rt-3d_q_%d.h5" % i)
    q = fh.root.StructGridField
    X, Y, Z = getMeshGrid(fh.root.StructGrid)
    tm = fh.root.timeData._v_attrs.vsTime

    f = figure(1)
    ax = subplot(1,2,1)
    im = pcolormesh(X, Y, q[:,:,Z.shape[0]/2,0].transpose())
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_xticks([X[0], 0, X[-1]])
    ax.set_xticklabels(['-0.25', '0.0', '0.25'])
    axis('image')
    colorbar_adj(im)

    ax = subplot(1,2,2)
    im = ax.pcolormesh(Z, Y, q[X.shape[0]/2,:,:,0])
    ax.set_xlabel('Z')
    ax.set_ylabel('Y')
    ax.set_xticks([Z[0], 0, Z[-1]])
    ax.set_xticklabels(['-0.25', '0.0', '0.25'])

    axis('image')
    colorbar_adj(im)

    suptitle("T = %.4g" % tm)
    
    savefig('s446-euler-rt-rho_%05d.png' % i, bbox_inches='tight')
    close()
    fh.close()
