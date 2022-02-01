from pylab import *
import tables

import pylab
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# customization for figure
#rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 14
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

gasGamma = 1.4

fileNames = ["../s408/s408-riemann-euler-3d_q_2.h5",
             "../s409/s409-riemann-euler-3d_q_2.h5",
             "../s410/s410-riemann-euler-3d_q_2.h5",
             "s411-riemann-euler-rz_q_10.h5"]

count = 1
figure(1)
for f in fileNames:
    fh = tables.openFile(f)
    grid = fh.root.StructGrid
    nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
    lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
    dx = lx/nx
    dz = lz/nz
    Xe = linspace(0.0, lx, nx+1)
    Ze = linspace(0.0, lz, nz+1)

    Xc = linspace(0.5*dx, lx-0.5*dx, nx)
    Zc = linspace(0.5*dz, lz-0.5*dz, nz)

    q = fh.root.StructGridField
    rho = q[:,:,:,0]
    u = q[:,:,:,1]/rho
    v = q[:,:,:,2]/rho
    pr = (q[:,:,:,4] - 0.5*rho*(u**2+v**2))*(gasGamma-1)

    # plot it
    subplot(2,2,count)
    im = pcolormesh(Xe, Ze, transpose(pr[:,0,:]), vmin=0.7, vmax=1.6)
    contour(Xc, Zc, transpose(pr[:,0,:]), 30, colors='k')
    axis('image')

    count = count+1

savefig('euler-spherical-riemann-cmp.png', dpi=300, bbox_inches='tight')

# function to compute lineout along ray with angle theta
def calcLineOut(nv, dx, dy, rho, theta):
    vals = zeros((nv, 1), float)
    cost = cos(theta)
    sint = sin(theta)
    dv = 1.0/nv
    vp = linspace(0.5*dv, 1.5-0.5*dv, nv)
    xp = cost*vp
    yp = sint*vp
    for i in range(nv):
        ix = (int) (xp[i]/dx)
        iy = (int) (yp[i]/dy)
        vals[i] = rho[ix,iy]
    return vals

# scatter plots
fh = tables.openFile("s411-riemann-euler-rz_q_10.h5")
grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx_ax = lx/nx
dy_ax = ly/ny
dz_ax = lz/nz

q = fh.root.StructGridField
rho = q[:,:,:,0]
u = q[:,:,:,1]/rho
v = q[:,:,:,2]/rho
pr = (q[:,:,:,4] - 0.5*rho*(u**2+v**2))*(gasGamma-1)

iz = 0.4/dz_ax

nv = 100
Xax = linspace(0, 1.5, nx)
pr_ax = pr[:,0,iz]

#fileNames = ["../s409/s409-riemann-euler-3d_q_2.h5",
#             "../s410/s410-riemann-euler-3d_q_2.h5"]

fileNames = ["../s408/s408-riemann-euler-3d_q_2.h5",
             "../s409/s409-riemann-euler-3d_q_2.h5",
             "../s410/s410-riemann-euler-3d_q_2.h5",
             "s411-riemann-euler-rz_q_10.h5"]
figure(2)
count = 1
for f in fileNames:
    fh = tables.openFile(f)
    grid = fh.root.StructGrid
    nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
    lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
    dx = lx/nx
    dy = ly/ny
    dz = lz/nz

    q = fh.root.StructGridField
    rho = q[:,:,:,0]
    u = q[:,:,:,1]/rho
    v = q[:,:,:,2]/rho
    pr = (q[:,:,:,4] - 0.5*rho*(u**2+v**2))*(gasGamma-1)

    iz = 0.4/dz

    nv = 40
    Xr = linspace(0, 1.5, nv)
    tRange = linspace(0, pi/4, 20)
    subplot(2,2,count)
    if f != "s411-riemann-euler-rz_q_10.h5":
        for i in range(tRange.shape[0]):
            v = calcLineOut(nv, dx, dy, pr[:,:,iz], tRange[i])
            plot(Xr, v, '-b')
    plot(Xax, pr_ax, '-k', linewidth=2)
    axis('tight')

    count = count+1

savefig('euler-spherical-riemann-lineout.png', dpi=300, bbox_inches='tight')

show()
