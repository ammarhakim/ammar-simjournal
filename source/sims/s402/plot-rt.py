from pylab import *
import tables

import pylab
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# customization for figure
#rcParams['lines.linewidth']            = 2
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

fh = tables.openFile("s402-euler-rt-ds-2d_q_5.h5")
grid = fh.root.StructGrid
nx, ny = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1]
lx = 1.0/6.0
ly = 1.0
dx = lx/nx
dy = ly/ny

def rhoExtend(rho):
    nx, ny = rho.shape[0], rho.shape[1]
    rhoExt = zeros((2*nx, ny), float)
    for j in range(ny):
        for i in range(nx):
            rhoExt[i,j] = rho[i,j]
        for i in range(nx):
            rhoExt[i+nx,j] = rho[nx-i-1,j]
    return rhoExt

q = fh.root.StructGridField
rho = q[:,:,0]
rhoExt = rhoExtend(rho)

X = linspace(0, 2*lx, 2*nx+1)
Y = linspace(0, ly, ny+1)

Xc = linspace(0.5*dx, 2*lx-0.5*dx, 2*nx)
Yc = linspace(0.5*dy, ly-0.5*dy, ny)

figure(1)
gs = gridspec.GridSpec(1, 2, height_ratios=[1, 1])

subplot(1, 2, 1)
im = pcolormesh(X, Y, transpose(rhoExt))
axis('image')
xticks( arange(0.0, 0.4, 0.1) )
colorbar_adj(im)

subplot(1, 2, 2)
contour(Xc, Yc, transpose(rhoExt), [1.5], colors='k', linewidth=2)
yticks( arange(0, 0, 0.2) )
xticks( arange(0, 0, 0.1) )
axis('image')
tight_layout()

savefig('s402-rt.png', dpi=300)#, bbox_inches='tight')
show()
