from pylab import *
import tables
import numpy

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

def getVars(q):
    gasGamma = 1.4
    rho = q[:,:,:,0]
    u = q[:,:,:,1]/rho
    v = q[:,:,:,2]/rho
    w = q[:,:,:,3]/rho
    Er = q[:,:,:,4]
    ke = 0.5*rho*(u*u+v*v+w*w)
    pr = (gasGamma-1)*(Er - ke)
    cs = sqrt(gasGamma*pr/rho)
    mach = sqrt(u*u+v*v+w*w)/cs
    return rho, pr, mach

# 3D
fh = tables.openFile("s425-super-euler-3d_inOut.h5")
inOut = fh.root.StructGridField
grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx = lx/nx
dy = ly/ny
dz = lz/nz
Xe = linspace(0.0, lx, nx+1)
Ze = linspace(0.0, lz, nz+1)

Xc = linspace(0.5*dx, lx-0.5*dx, nx)
Zc = linspace(0.5*dz, lz-0.5*dz, nz)

fh = tables.openFile("s425-super-euler-3d_q_10.h5")
q = fh.root.StructGridField
rho, pr, mach = getVars(q)

inOutA = inOut[:,0,:]
prMa = numpy.ma.masked_where(inOutA[:,:,0]<0, pr[:,0,:])
rhoMa = numpy.ma.masked_where(inOutA[:,:,0]<0, rho[:,0,:])

figure(1)
contour(Xc, Zc, transpose(prMa), 20, colors='k')
im = pcolormesh(Xc, Zc, transpose(rhoMa))
xlabel('X')
ylabel('Z')
axis('image')
colorbar_adj(im)

savefig('s425-rho-pr-3d.png', dpi=300, bbox_inches='tight')

show()
